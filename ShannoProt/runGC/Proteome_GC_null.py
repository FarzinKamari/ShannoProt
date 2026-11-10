#!/usr/bin/env python
# coding: utf-8

# # CEI (formerly CAIR) expected under random translation at the genome GC content
# * This module matches each proteome of Table S1 to the GC content of a genome assembly of the same organism (NCBI), takes the genetic code of the organism from the NCBI taxonomy, and calculates the CEI expected if codons were drawn at random from that GC content (the null model of H2, as described in Methods section).
# * Null model: P(G) = P(C) = GC/2 and P(A) = P(T) = (1 - GC)/2 at each codon position; a codon's probability is the product over its three positions; stop codons of the organism's genetic code are dropped; amino-acid probabilities are renormalized; CEI_null = -sum(p * log22(p)).
# * Matching: the exact taxid (Organism_ID) first, then assemblies of the same species (species_taxid); within each, RefSeq reference > representative > other genomes, RefSeq before GenBank, the most complete assembly level and the latest release. Current assemblies are preferred to historical (suppressed) ones. Taxids merged since 2020 are mapped through merged.dmp.
# ##### Input files (see README):
# >* Complete_proteome_CAIRs.csv (Table S1)
# >* assembly_summary_refseq.txt, assembly_summary_genbank.txt and their historical versions from https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
# >* nodes.dmp and merged.dmp from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# ##### Output files:
# >* GC null per proteome.csv, GC null summary by domain.csv, GC null summary by phylum.csv, GC null curve.csv, GC null maxima.csv
# >* GC null distribution shape.csv (the left skew of the CEI distribution and of the null)

# In[1]:


import numpy as np
import pandas as pd
from itertools import product
from Bio.Data import CodonTable # requires biopython package


# In[2]:


def cei_null(gc, gencode=11, gc_pos=None): # refer to proteome_gc_null
    g = gc_pos if gc_pos is not None else (gc, gc, gc) # GC at codon positions 1, 2 and 3 (uniform unless given)
    table = CodonTable.unambiguous_dna_by_id[int(gencode)] # NCBI genetic code
    aa_prob = {} # probability of each amino acid
    for codon in (''.join(c) for c in product('TCAG', repeat=3)):
        if codon in table.stop_codons:
            continue # stop codons are dropped
        prob = np.prod([g[i] / 2 if base in 'GC' else (1 - g[i]) / 2 for i, base in enumerate(codon)]) # codon probability
        aa = table.forward_table[codon]
        aa_prob[aa] = aa_prob.get(aa, 0) + prob
    p = np.array(list(aa_prob.values()))
    p = p[p > 0] / p.sum() # renormalized over the coded amino acids (0*log0 is taken as 0)
    return float(-np.sum(p * np.log(p)) / np.log(22)) # CEI_null in the same base-22 units as CEI


def read_taxonomy(nodes_file='nodes.dmp', merged_file='merged.dmp'): # parents, ranks and genetic codes of all taxids
    parent, rank, gencode = {}, {}, {}
    with open(nodes_file) as fileread:
        for line in fileread:
            f = line.split('\t|\t') # NCBI dump format
            parent[f[0]] = f[1]
            rank[f[0]] = f[2]
            gencode[f[0]] = f[6] # genetic code id
    merged = {}
    with open(merged_file) as fileread:
        for line in fileread:
            f = line.rstrip('\t|\n').split('\t|\t')
            merged[f[0]] = f[1] # old taxid -> new taxid
    return parent, rank, gencode, merged


def species_of(taxid, parent, rank): # walks up the taxonomy to the species node
    t = taxid
    while t in rank and rank[t] != 'species' and parent[t] != t:
        t = parent[t]
    return t if rank.get(t) == 'species' else None


def read_assemblies(files): # the needed columns of the NCBI assembly summaries
    cols = ['#assembly_accession', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'assembly_level',
            'seq_rel_date', 'gc_percent', 'ftp_path', 'annotation_provider']
    out = []
    for name, source, historical in files:
        a = pd.read_csv(name, sep='\t', skiprows=1, usecols=cols, dtype=str, quoting=3) # the first line is a comment
        a['source'] = source # RefSeq or GenBank
        a['historical'] = historical # suppressed (historical) assemblies
        out.append(a[pd.to_numeric(a.gc_percent, errors='coerce') > 0]) # only assemblies with a GC value (some historical ones have 0, i.e. missing)
    out = pd.concat(out).rename(columns={'#assembly_accession': 'assembly_accession'})
    # ranking of candidate assemblies (lower is better), as described in Methods section
    out['rank_category'] = out.refseq_category.map({'reference genome': 0, 'representative genome': 1}).fillna(2)
    out['rank_source'] = (out.source != 'RefSeq').astype(int)
    out['rank_level'] = out.assembly_level.map({'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}).fillna(4)
    out = out.sort_values(['historical', 'rank_category', 'rank_source', 'rank_level', 'seq_rel_date'],
                          ascending=[True, True, True, True, False]) # best candidate first
    return out


# In[3]:


def proteome_gc_null(cair_file='Complete_proteome_CAIRs.csv', assembly_files=None, nodes_file='nodes.dmp',
                     merged_file='merged.dmp', outfile='GC null per proteome.csv'): # refer to runGC
    if assembly_files is None:
        assembly_files = [('assembly_summary_refseq.txt', 'RefSeq', False), ('assembly_summary_genbank.txt', 'GenBank', False),
                          ('assembly_summary_refseq_historical.txt', 'RefSeq', True),
                          ('assembly_summary_genbank_historical.txt', 'GenBank', True)]
    prot = pd.read_csv(cair_file) # Table S1
    prot['Domain'] = np.where(prot.code <= 56, 'Bacteria', np.where(prot.code <= 76, 'Archaea', 'Eukaryota')) # codes of the tree of life
    parent, rank, gencode, merged = read_taxonomy(nodes_file, merged_file)
    asm = read_assemblies(assembly_files)
    best_taxid = asm.drop_duplicates('taxid') # the best assembly of each taxid
    best_species = asm.drop_duplicates('species_taxid').set_index('species_taxid') # the best assembly of each species
    best_taxid = best_taxid.set_index('taxid')

    rows = []
    for org in prot.Organism_ID.astype(int).astype(str):
        taxid = merged.get(org, org) # current taxid
        sp = species_of(taxid, parent, rank) if taxid in rank else None # species taxid
        if taxid in best_taxid.index:
            hit, match = best_taxid.loc[taxid], 'taxid'
        elif sp is not None and sp in best_species.index:
            hit, match = best_species.loc[sp], 'species'
        else:
            hit, match = None, 'none'
        if hit is not None and hit.historical:
            match += ' (historical assembly)'
        code = gencode.get(taxid, gencode.get(sp, None)) if (taxid in rank or sp) else None # genetic code of the organism
        rows.append({'Organism_ID': int(org),
                     'Current taxid': taxid if taxid in rank else '',
                     'Match type': match,
                     'Assembly accession': '' if hit is None else hit.assembly_accession,
                     'Assembly level': '' if hit is None else hit.assembly_level,
                     'RefSeq category': '' if hit is None else hit.refseq_category,
                     'Assembly release date': '' if hit is None else hit.seq_rel_date,
                     'Annotation provider': '' if hit is None else hit.annotation_provider,
                     'FTP path': '' if hit is None else hit.ftp_path,
                     'GC': np.nan if hit is None else float(hit.gc_percent) / 100,
                     'Genetic code': code})
    out = prot.merge(pd.DataFrame(rows), on='Organism_ID', how='left') # adding the matches to Table S1
    missing_code = out['Genetic code'].isna() # taxids no longer in the taxonomy: code of the domain
    out.loc[missing_code, 'Genetic code'] = np.where(out.loc[missing_code, 'Domain'] == 'Eukaryota', '1', '11')
    out['Genetic code'] = out['Genetic code'].astype(int)
    nulls = {} # CEI_null is calculated once for each (GC, code) pair
    out['CEI null'] = [nulls.setdefault((gc, c), cei_null(gc, c)) if pd.notna(gc) else np.nan
                       for gc, c in zip(out.GC, out['Genetic code'])]
    out['Delta CEI'] = out.CAIR - out['CEI null'] # observed minus expected CEI
    out.rename(columns={'CAIR': 'CEI'}).to_csv(outfile, index=False) # writing the CSV output
    print(len(out), 'proteomes;', int(out.GC.notna().sum()), 'matched to an assembly GC; written to', outfile)
    return out



# In[4]:


def gc_null_summary(gc_file='GC null per proteome.csv', outfile='GC null summary by domain.csv',
                    phylum_outfile='GC null summary by phylum.csv'): # refer to runGC
    import statsmodels.api as sm
    from scipy import stats
    g = pd.read_csv(gc_file) # output of proteome_gc_null
    ceiling = max(cei_null(gc, 11) for gc in np.arange(0.10, 0.81, 0.0001)) # the highest CEI the uniform-GC null can reach
    rows = []
    for domain, d in [('All', g)] + list(g.groupby('Domain')):
        m = d[d.GC.notna()] # proteomes with a GC value
        fit = sm.OLS(m.CEI, sm.add_constant(pd.DataFrame({'GC': m.GC, 'GC2': m.GC ** 2}))).fit() # CEI ~ GC + GC^2
        rows.append({'Domain': domain,
                     'Proteomes': len(d),
                     'With GC': len(m),
                     'R2 of CEI ~ GC + GC2': fit.rsquared,
                     'Spearman CEI vs CEI null': stats.spearmanr(m.CEI, m['CEI null'])[0],
                     'Above the null (Delta CEI > 0)': int((m['Delta CEI'] > 0).sum()),
                     'Share above the null': float((m['Delta CEI'] > 0).mean()),
                     'Median Delta CEI': float(m['Delta CEI'].median()),
                     'Above the uniform-GC ceiling': int((d.CEI > ceiling).sum()),
                     'Share above the ceiling': float((d.CEI > ceiling).mean()),
                     'Ceiling': ceiling})
    out = pd.DataFrame(rows)
    out.to_csv(outfile, index=False) # writing the CSV output
    m = g[g.GC.notna()]
    phyla = m.groupby(['Domain', 'Second_hierarchy']).agg(Proteomes=('CEI', 'size'), Median_GC=('GC', 'median'),
                                                          Median_CEI=('CEI', 'median'), Median_CEI_null=('CEI null', 'median'),
                                                          Median_Delta_CEI=('Delta CEI', 'median'),
                                                          Share_above_null=('Delta CEI', lambda x: float((x > 0).mean()))).reset_index()
    phyla.to_csv(phylum_outfile, index=False) # writing the CSV output
    gcs = np.round(np.arange(0.10, 0.8001, 0.001), 3) # GC values of the null curve (Figure 3A)
    pd.DataFrame({'GC': gcs, 'CEI null (code 11)': [cei_null(x, 11) for x in gcs],
                  'CEI null (code 4)': [cei_null(x, 4) for x in gcs]}).to_csv('GC null curve.csv', index=False)
    print(out.to_string())
    return out, phyla


# In[5]:


def null_maxima(gencode=11, starts=60, outfile='GC null maxima.csv'): # the highest CEI each null model can reach (Methods, Section 3.2)
    from scipy.optimize import minimize
    rng = np.random.default_rng(2026)
    def free_cei(x): # x: 12 logits -> probabilities of T, C, A, G at codon positions 1-3 (no constraint between bases)
        p = np.exp(x.reshape(3, 4))
        p = p / p.sum(axis=1, keepdims=True)
        table = CodonTable.unambiguous_dna_by_id[int(gencode)]
        aa_prob = {}
        for codon in (''.join(c) for c in product('TCAG', repeat=3)):
            if codon in table.stop_codons:
                continue
            prob = np.prod([p[i, 'TCAG'.index(base)] for i, base in enumerate(codon)])
            aa = table.forward_table[codon]
            aa_prob[aa] = aa_prob.get(aa, 0) + prob
        q = np.array(list(aa_prob.values()))
        q = q[q > 0] / q.sum()
        return float(-np.sum(q * np.log(q)) / np.log(22))
    uniform = max((cei_null(g / 1000, gencode), g / 1000) for g in range(1001)) # same GC at all codon positions
    by_position = minimize(lambda x: -cei_null(None, gencode, gc_pos=tuple(np.clip(x, 0, 1))), [0.45, 0.4, 0.9], method='Nelder-Mead',
                           options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 4000}) # GC of each codon position (G = C, A = T)
    free = min((minimize(lambda x: -free_cei(x), rng.normal(0, 1, 12), method='Nelder-Mead',
                         options={'maxiter': 20000, 'xatol': 1e-7, 'fatol': 1e-10}) for _ in range(starts)), key=lambda r: r.fun)
    out = pd.DataFrame({'Null model': ['uniform GC', 'GC by codon position', 'free base composition by codon position'],
                        'Highest CEI': [uniform[0], -by_position.fun, -free.fun],
                        'At': ['GC %.3f' % uniform[1], 'GC1-GC3 %s' % np.round(np.clip(by_position.x, 0, 1), 3).tolist(), '']})
    out.to_csv(outfile, index=False) # writing the CSV output
    print(out.to_string())
    return out


# In[6]:


def cei_distribution_shape(gc_file='GC null per proteome.csv', outfile='GC null distribution shape.csv'): # refer to runGC
    # is the left skew of the proteome CEI distribution (Section 3.1) expected from GC content alone?
    from scipy import stats
    g = pd.read_csv(gc_file) # output of proteome_gc_null
    m = g[g.GC.notna()] # proteomes with a GC value
    low = m[m.CEI <= m.CEI.quantile(0.01)] # the 1% least even proteomes
    extreme = lambda x: float(((x.GC < 0.30) | (x.GC > 0.65)).mean()) # share of AT-rich (GC < 30%) or GC-rich (GC > 65%) genomes
    out = pd.DataFrame({'Quantity': ['Skewness of CEI (all proteomes)', 'Skewness of CEI (proteomes with GC)', 'Skewness of the CEI null',
                                     'Proteomes with GC', 'Least even 1% (proteomes)', 'Least even 1%: share with GC < 30% or > 65%',
                                     'All proteomes with GC: share with GC < 30% or > 65%'],
                        'Value': [stats.skew(g.CEI), stats.skew(m.CEI), stats.skew(m['CEI null']),
                                  len(m), len(low), extreme(low), extreme(m)]})
    out.to_csv(outfile, index=False) # writing the CSV output
    print(out.to_string())
    return out
