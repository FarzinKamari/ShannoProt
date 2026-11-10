#!/usr/bin/env python
# coding: utf-8

# # Codon-position GC of one annotated genome per genus -- sensitivity of the GC null
# * The null model of Proteome_GC_null.py assumes the same GC at the three codon positions. Real genes have position-specific GC (GC1, GC2, GC3), which raises the CEI that random translation can reach.
# * For one representative annotated assembly per genus (genus = first word of the organism name, as in the genus-level tests), this module downloads the coding sequences (NCBI "cds_from_genomic" file), counts GC at codon positions 1, 2 and 3, and recalculates the null CEI with position-specific GC.
# * Representative: among the assemblies matched to the proteomes of the genus (GC null per proteome.csv), the annotated one ranked best (RefSeq reference > representative > other; RefSeq before GenBank; most complete level; latest).
# ##### Input file:
# >* GC null per proteome.csv (from Proteome_GC_null.py)
# ##### Output files:
# >* cds/ folder with the downloaded files (not needed after the counts)
# >* Codon position GC per genus.csv

# In[1]:


import os
import gzip
import time
import urllib.request
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from Proteome_GC_null import cei_null


# In[2]:


def genus_of(name): # first word of the organism name, after "Candidatus", "uncultured", "unclassified" (as in Bifurcation tests.R)
    words = str(name).replace('[', '').replace(']', '').replace("'", '').replace('"', '').split()
    while words and words[0].lower() in ('candidatus', 'uncultured', 'unclassified'):
        words = words[1:]
    return words[0] if words else 'NA'


def representatives(gc_file='GC null per proteome.csv'): # one annotated assembly per genus
    g = pd.read_csv(gc_file, dtype={'Annotation provider': str})
    g['Genus'] = g.Organism.apply(genus_of)
    g = g[g['Annotation provider'].notna() & g['FTP path'].notna() & (g['FTP path'] != 'na')] # annotated assemblies only
    g['rank_category'] = g['RefSeq category'].map({'reference genome': 0, 'representative genome': 1}).fillna(2)
    g['rank_source'] = (~g['Assembly accession'].str.startswith('GCF')).astype(int) # RefSeq (GCF) first
    g['rank_level'] = g['Assembly level'].map({'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}).fillna(4)
    g = g.sort_values(['rank_category', 'rank_source', 'rank_level', 'Assembly release date'], ascending=[True, True, True, False])
    return g.drop_duplicates('Genus') # the best assembly of each genus


def download(ftp_path, folder='cds', tries=3): # downloads the cds_from_genomic file of an assembly
    ftp_path = ftp_path.rstrip('/') # some paths of the assembly summaries end with "/"
    url = ftp_path.replace('ftp://', 'https://') + '/' + ftp_path.split('/')[-1] + '_cds_from_genomic.fna.gz'
    outfile = os.path.join(folder, url.split('/')[-1])
    if os.path.exists(outfile):
        return outfile # already downloaded
    for t in range(tries):
        try:
            urllib.request.urlretrieve(url, outfile + '.part')
            os.replace(outfile + '.part', outfile)
            return outfile
        except Exception as e:
            if '404' in str(e):
                return None # the assembly has no CDS file
            time.sleep(5 * (t + 1)) # waiting before trying again
    return None


def position_gc(cds_file): # GC at codon positions 1, 2, 3 over all complete CDS of a genome
    counts = np.zeros((3, 2)) # rows: positions; columns: GC, all
    seq = []
    with gzip.open(cds_file, 'rt') as fileread:
        for line in fileread:
            if line.startswith('>'):
                if seq:
                    counts += codon_counts(''.join(seq))
                seq = []
                skip = ('pseudo=true' in line) or ('partial=' in line) # pseudogenes and partial CDS are skipped
            elif not skip:
                seq.append(line.strip().upper())
        if seq:
            counts += codon_counts(''.join(seq))
    return counts[:, 0] / counts[:, 1], int(counts[0, 1]) # GC1, GC2, GC3 and number of codons


def codon_counts(s): # GC and ACGT counts at each codon position of one CDS
    out = np.zeros((3, 2))
    if len(s) % 3 != 0:
        return out # incomplete reading frame
    for p in range(3):
        sub = s[p::3]
        gc = sub.count('G') + sub.count('C')
        out[p] = [gc, gc + sub.count('A') + sub.count('T')]
    return out


# In[3]:


def codon_position_gc(gc_file='GC null per proteome.csv', folder='cds', outfile='Codon position GC per genus.csv', threads=4): # refer to runGC
    os.makedirs(folder, exist_ok=True)
    rep = representatives(gc_file)
    with ThreadPoolExecutor(threads) as pool: # a few downloads at a time
        files = list(pool.map(lambda p: download(p, folder), rep['FTP path']))
    rows = []
    for (_, r), f in zip(rep.iterrows(), files):
        if f is None:
            continue
        (gc1, gc2, gc3), n_codons = position_gc(f)
        if n_codons == 0:
            continue
        rows.append({'Genus': r.Genus, 'Domain': r.Domain, 'Second_hierarchy': r.Second_hierarchy, 'Organism_ID': r.Organism_ID,
                     'Organism': r.Organism, 'Assembly accession': r['Assembly accession'], 'Genetic code': r['Genetic code'],
                     'CEI': r.CEI, 'Genome GC': r.GC, 'Codons': n_codons, 'GC1': gc1, 'GC2': gc2, 'GC3': gc3,
                     'CDS GC': (gc1 + gc2 + gc3) / 3,
                     'CEI null (genome GC)': r['CEI null'],
                     'CEI null (CDS GC, uniform)': cei_null((gc1 + gc2 + gc3) / 3, r['Genetic code']),
                     'CEI null (GC1, GC2, GC3)': cei_null(None, r['Genetic code'], gc_pos=(gc1, gc2, gc3))})
    out = pd.DataFrame(rows)
    out.to_csv(outfile, index=False) # writing the CSV output
    print(len(rep), 'genera with an annotated assembly;', len(out), 'with coding sequences; written to', outfile)
    return out



# In[4]:


def codon_position_summary(cds_file='Codon position GC per genus.csv', ceiling=None, outfile='Codon position GC summary.csv'): # refer to runGC
    cp = pd.read_csv(cds_file)
    if ceiling is None:
        ceiling = max(cei_null(g / 1000) for g in range(100, 801)) # the highest CEI of the uniform-GC null (standard bacterial code)
    rows = [('genera', len(cp)), ('r genome GC vs CDS GC', cp[['Genome GC', 'CDS GC']].corr().iloc[0, 1]),
            ('mean GC1', cp.GC1.mean()), ('mean GC2', cp.GC2.mean()), ('mean GC3', cp.GC3.mean())]
    nulls = {'uniform genome GC': 'CEI null (genome GC)', 'codon-position GC': 'CEI null (GC1, GC2, GC3)'}
    for name, col in nulls.items():
        rows.append(('R2 of CEI with the null, ' + name + ', all', np.corrcoef(cp.CEI, cp[col])[0, 1] ** 2)) # squared correlation
    for d, x in cp.groupby('Domain'):
        rows.append(('genera, ' + d, len(x)))
        for name, col in nulls.items():
            rows += [('R2 of CEI with the null, ' + name + ', ' + d, np.corrcoef(x.CEI, x[col])[0, 1] ** 2),
                     ('share above the null, ' + name + ', ' + d, (x.CEI > x[col]).mean()),
                     ('median delta CEI, ' + name + ', ' + d, (x.CEI - x[col]).median())]
    euk = cp[(cp.Domain == 'Eukaryota') & (cp.CEI > ceiling)]
    rows += [('eukaryotic genera above the uniform-GC ceiling', len(euk)),
             ('... also above their codon-position null', int((euk.CEI > euk['CEI null (GC1, GC2, GC3)']).sum())),
             ('highest codon-position null', cp['CEI null (GC1, GC2, GC3)'].max()), ('highest CEI of these genera', cp.CEI.max())]
    out = pd.DataFrame(rows, columns=['item', 'value'])
    out.to_csv(outfile, index=False) # writing the CSV output
    print(out.to_string())
    return out
