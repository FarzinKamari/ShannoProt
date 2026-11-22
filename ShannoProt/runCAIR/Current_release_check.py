#!/usr/bin/env python
# coding: utf-8

# # CEI in the current UniProt reference proteomes -- check of the 2020 proteome values (rev 26)
# * The proteome CEI values of Table S1 were computed from all entries of each organism in UniProt release 2020_02 (runCAIR). This module recomputes CEI for a sample of these organisms from the current UniProt reference proteomes: one canonical sequence per gene (no isoforms), with and without the sequences flagged as fragments, counted as in runCAIR (22 symbols; B, Z and J split between their two residues; X not counted).
# * Sample: up to four organisms per second-rank lineage of Table S1 whose NCBI taxid has a current reference proteome (random, seed 2026), together with the organisms of Table 1.
# * The README file and the FASTA files of the reference proteomes are downloaded from the UniProt FTP server if they are not in the working folder (about 400 MB for the sample).
# ##### Input file:
# >* GC null per proteome.csv (from Proteome_GC_null.py in runGC; the CEI, domain and lineage of each proteome of Table S1)
# ##### Output files:
# >* reference_proteomes/ folder with the README and the downloaded FASTA files
# >* Current release check.csv (one row per organism of the sample) and Current release check summary.json

# In[1]:


import os
import gzip
import json
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats
ftp = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes'
symbols = 'ACDEFGHIKLMNOPQRSTUVWY' # the 22 counted symbols (20 standard amino acids, selenocysteine and pyrrolysine)
table1 = [3702, 6239, 3055, 51511, 7955, 7227, 83333, 9606, 3880, 10090, 367110, 10116, 559292, 284812, 312017, 8364, 4577, 871271] # Table 1


# In[2]:


def cei_of(counts): # CEI (Shannon entropy in base 22) of residue counts
    c = np.array([counts[a] for a in symbols], dtype=float)
    p = c[c > 0] / c.sum()
    return float(-(p * np.log(p)).sum() / np.log(22))


def composition(path): # residue counts of a FASTA file, all sequences and without the fragments
    tot = {a: 0.0 for a in symbols}
    nofrag = {a: 0.0 for a in symbols}
    n = nf = 0
    seq, frag = [], False

    def add(seq, frag):
        s = ''.join(seq)
        for a in symbols:
            k = s.count(a)
            tot[a] += k
            if not frag:
                nofrag[a] += k
        for amb, (x, y) in (('B', ('D', 'N')), ('Z', ('E', 'Q')), ('J', ('I', 'L'))): # ambiguity codes split between two residues
            k = s.count(amb) / 2
            for t in (x, y):
                tot[t] += k
                if not frag:
                    nofrag[t] += k
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if seq:
                    add(seq, frag)
                seq, frag = [], '(Fragment)' in line
                n += 1
                nf += frag
            else:
                seq.append(line.strip())
        if seq:
            add(seq, frag)
    return tot, nofrag, n, nf


def read_readme(path): # the table of reference proteomes in the README file, and the release
    lines = open(path, encoding='utf-8').read().split('\n')
    release = next(l.strip() for l in lines if l.startswith('Release '))
    i0 = next(i for i, l in enumerate(lines) if l.startswith('Proteome_ID'))
    ref = pd.DataFrame([l.split('\t') for l in lines[i0 + 1:] if l.count('\t') >= 7], columns=lines[i0].split('\t'))
    ref['Tax_ID'] = ref.Tax_ID.astype(int)
    return ref, release


# In[3]:


def current_release_check(gc_file='GC null per proteome.csv', folder='reference_proteomes', outfile='Current release check.csv',
                          summary_file='Current release check summary.json'): # refer to runCAIR
    os.makedirs(folder, exist_ok=True)
    readme = os.path.join(folder, 'README')
    if not os.path.exists(readme):
        urllib.request.urlretrieve(ftp + '/README', readme)
    ref, release = read_readme(readme)
    s1 = pd.read_csv(gc_file, low_memory=False)
    s1['taxid'] = s1['Current taxid'].fillna(s1.Organism_ID).astype(int)
    m = s1.merge(ref, left_on='taxid', right_on='Tax_ID') # organisms whose taxid has a current reference proteome
    pick = pd.concat([g.sample(min(len(g), 4), random_state=2026) for _, g in m.groupby('Second_hierarchy')] +
                     [m[m.taxid.isin(table1)]]).drop_duplicates('taxid')
    print('organisms of Table S1 with a current reference proteome of the same taxid:', len(m), '| sample:', len(pick))
    rows = []
    for _, r in pick.iterrows():
        reg = {'archaea': 'Archaea', 'bacteria': 'Bacteria', 'eukaryota': 'Eukaryota'}.get(r.SUPERREGNUM)
        if reg is None:
            continue
        fn = '%s_%d.fasta.gz' % (r.Proteome_ID, r.Tax_ID)
        path = os.path.join(folder, fn)
        if not os.path.exists(path):
            urllib.request.urlretrieve('%s/%s/%s/%s' % (ftp, reg, r.Proteome_ID, fn), path)
        tot, nofrag, n, nf = composition(path)
        rows.append({'taxid': r.taxid, 'Organism': r.Organism, 'Domain': r.Domain, 'Second_hierarchy': r.Second_hierarchy,
                     'Proteome_ID': r.Proteome_ID, 'proteins (current, canonical)': n, 'fragments': nf, 'CEI 2020': r.CEI,
                     'CEI current': cei_of(tot), 'CEI current, no fragments': cei_of(nofrag)})
    res = pd.DataFrame(rows)
    res.to_csv(outfile, index=False) # writing the CSV output
    R = {'UniProt release': release, 'organisms matched': len(m), 'sample': len(res), 'sample by domain': res.Domain.value_counts().to_dict(),
         'lineages in the sample': int(res.Second_hierarchy.nunique())}
    for c in ('CEI current', 'CEI current, no fragments'):
        dlt = res[c] - res['CEI 2020']
        R[c] = {'Spearman with 2020': round(float(stats.spearmanr(res['CEI 2020'], res[c])[0]), 4),
                'Pearson with 2020': round(float(stats.pearsonr(res['CEI 2020'], res[c])[0]), 4),
                'median difference (current - 2020)': float('%.2g' % dlt.median()),
                'median |difference|': float('%.2g' % dlt.abs().median()), '95th percentile |difference|': float('%.2g' % dlt.abs().quantile(0.95)),
                'max |difference|': float('%.2g' % dlt.abs().max())}
    R['identical to 1e-6, by domain'] = {k: [int(((v['CEI current'] - v['CEI 2020']).abs() < 1e-6).sum()), len(v)] for k, v in res.groupby('Domain')}
    R['fragments: median share of canonical sequences'] = float('%.3g' % (res.fragments / res['proteins (current, canonical)']).median())
    R['by domain: median CEI 2020 / current'] = {k: [round(float(v['CEI 2020'].median()), 4), round(float(v['CEI current'].median()), 4)]
                                                 for k, v in res.groupby('Domain')}
    json.dump(R, open(summary_file, 'w'), indent=1, default=str) # writing the JSON output
    for k, v in R.items():
        print(k, ':', v)
    return res
