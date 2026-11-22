#!/usr/bin/env python
# coding: utf-8

# # Protein length and other properties of the genes -- covariates of the length association (rev 26)
# * Is the association of protein length with rare-disease involvement explained by the number of coding exons, expression breadth, tissue specificity, expression level, essentiality or loss-of-function constraint, and how does length relate to the size of the mutational target of a gene?
# * Exon count, expected numbers of variants under the mutation model of gnomAD (missense, loss-of-function) and loss-of-function constraint (LOEUF, upper bound of the observed/expected ratio): gnomAD v2.1.1 constraint table (canonical transcripts).
# * Expression breadth (number of tissues with a median TPM of at least 1), tissue specificity tau (on log2(TPM + 1)) and expression level (highest median TPM): GTEx v8 median TPM of 54 tissues.
# * Essentiality: core essential genes CEGv2 (Hart et al. 2017).
# * Genes are matched to the reviewed human proteins through the Ensembl gene and GeneID cross-references of UniProt release 2020_02; where an entry has several Ensembl genes, the gnomAD transcript whose coding length best matches the protein length is used.
# ##### Input files:
# >* Human analysis data.csv (from Human_indicator_statistics.py) and Human proteome 2020_02.tab (from Human_proteome_2020_02.py)
# >* gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz and CEGv2.txt
# ##### Output file:
# >* Length covariates.json

# In[1]:


import json
import numpy as np
import pandas as pd
from scipy import stats
from Human_indicator_statistics import auc, logistic, odds_ratios


# In[2]:


def gene_covariates(proteome_file='Human proteome 2020_02.tab', data_file='Human analysis data.csv',
                    gnomad_file='gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
                    gtex_file='GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz', essential_file='CEGv2.txt'):
    d = pd.read_csv(data_file, low_memory=False)
    prot = pd.read_csv(proteome_file, sep='\t', usecols=['Entry', 'Status', 'GeneID', 'Ensembl genes'])
    prot = prot[prot.Status == 'reviewed']
    # gnomAD canonical transcripts (protein-coding genes)
    gn = pd.read_csv(gnomad_file, sep='\t', compression='gzip',
                     usecols=['gene', 'gene_id', 'cds_length', 'num_coding_exons', 'exp_mis', 'exp_lof', 'exp_syn', 'oe_lof_upper', 'pLI', 'gene_length'])
    ens = prot[['Entry', 'Ensembl genes']].dropna()
    ens = ens.assign(gene_id=ens['Ensembl genes'].str.split(';')).explode('gene_id')
    ens['gene_id'] = ens.gene_id.str.strip().str.split('.').str[0]
    ens = ens[ens.gene_id.str.startswith('ENSG')].merge(gn, on='gene_id').merge(d[['Entry', 'Protein length']], on='Entry')
    ens['len_diff'] = (ens.cds_length / 3 - 1 - ens['Protein length']).abs() # coding length (without the stop codon) vs protein length
    ens = ens.sort_values(['Entry', 'len_diff']).drop_duplicates('Entry')
    # GTEx: expression breadth, tissue specificity and level
    gt = pd.read_csv(gtex_file, sep='\t', skiprows=2)
    gt['gene_id'] = gt.Name.str.split('.').str[0]
    tpm = gt.drop(columns=['Name', 'Description']).groupby('gene_id').max()
    lx = np.log2(tpm + 1)
    gtx = pd.DataFrame({'tissues TPM>=1': (tpm >= 1).sum(axis=1),
                        'tau': ((1 - lx.div(lx.max(axis=1).replace(0, np.nan), axis=0)).sum(axis=1) / (lx.shape[1] - 1)),
                        'log2 max TPM': lx.max(axis=1)}).reset_index()
    eg = ens[['Entry', 'gene_id']].merge(gtx, on='gene_id', how='left')
    # core essential genes (CEGv2), by GeneID
    ceg = pd.read_csv(essential_file, sep='\t')
    gid = prot[['Entry', 'GeneID']].dropna()
    gid = gid.assign(GeneID=gid.GeneID.astype(str).str.split(';')).explode('GeneID')
    gid['GeneID'] = pd.to_numeric(gid.GeneID.str.strip(), errors='coerce')
    essential = set(gid[gid.GeneID.isin(set(ceg.ENTREZ_ID))].Entry)
    a = d[d.Category.notna()].copy() # the 18,471 proteins of the article
    a = a.merge(ens[['Entry', 'gene_id', 'cds_length', 'num_coding_exons', 'exp_mis', 'exp_lof', 'exp_syn', 'oe_lof_upper', 'pLI', 'gene_length']],
                on='Entry', how='left').merge(eg[['Entry', 'tissues TPM>=1', 'tau', 'log2 max TPM']], on='Entry', how='left')
    a['essential (CEGv2)'] = a.Entry.isin(essential).astype(float)
    a['log exons'] = np.log(a.num_coding_exons)
    a['log exp_mis'] = np.log(a.exp_mis)
    a['log exp_lof'] = np.log(a.exp_lof.clip(lower=0.1))
    a['log gene length'] = np.log(a.gene_length)
    a['LOEUF'] = a.oe_lof_upper
    return a


# In[3]:


def length_covariates(data_file='Human analysis data.csv', proteome_file='Human proteome 2020_02.tab',
                      gnomad_file='gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
                      gtex_file='GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz', essential_file='CEGv2.txt',
                      outfile='Length covariates.json'): # refer to runNEE
    a = gene_covariates(proteome_file, data_file, gnomad_file, gtex_file, essential_file)
    R = {'proteins analysed': len(a), 'matched to gnomAD': int(a.cds_length.notna().sum()), 'matched to GTEx': int(a['tissues TPM>=1'].notna().sum()),
         'core essential (CEGv2)': int(a['essential (CEGv2)'].sum()),
         'check: median |CDS/3 - 1 - protein length|': float((a.cds_length / 3 - 1 - a['Protein length']).abs().median())}
    sp = lambda x, y: float(stats.spearmanr(a[x], a[y], nan_policy='omit')[0])
    R['Spearman with protein length'] = {c: sp('Protein length', c) for c in ['exp_mis', 'exp_lof', 'num_coding_exons', 'gene_length', 'tissues TPM>=1',
                                                                             'tau', 'LOEUF', 'Number of interactions']}
    for c in ['Protein length', 'exp_mis', 'num_coding_exons', 'gene_length', 'tissues TPM>=1', 'tau', 'LOEUF', 'essential (CEGv2)']: # single covariates
        x = a[a[c].notna()]
        v = x[c] if c not in ('LOEUF', 'tau') else -x[c] # lower LOEUF = more constrained; lower tau = broader expression
        R.setdefault('AUC (single covariate, same proteins)', {})[c] = [float(auc(x.Disease, v)), float(auc(x.Disease, x['Protein length'])), len(x)]
    # joint logistic models of the proteins with a gene age and all covariates (with gene-age class), odds ratios per unit
    need = ['log L', 'log(1+degree)', 'CEI_MM', 'log exons', 'tissues TPM>=1', 'tau', 'LOEUF', 'essential (CEGv2)', 'log2 max TPM']
    g = a[a['Age rank'].notna()].dropna(subset=need).copy()
    g['tissues TPM>=1 (per 10)'] = g['tissues TPM>=1'] / 10
    R['n in the joint models'] = [len(g), int(g.Disease.sum())]
    base = ['log L', 'log(1+degree)', 'CEI_MM'] # the covariates of the joint model of the article
    models = {'article covariates': base,
              '+ log exon count': base + ['log exons'],
              '+ expression breadth and specificity': base + ['tissues TPM>=1 (per 10)', 'tau'],
              '+ expression level': base + ['log2 max TPM'],
              '+ essential (CEGv2)': base + ['essential (CEGv2)'],
              '+ LOEUF (LoF constraint)': base + ['LOEUF'],
              '+ all': base + ['log exons', 'tissues TPM>=1 (per 10)', 'tau', 'log2 max TPM', 'essential (CEGv2)', 'LOEUF'],
              'expected missense variants instead of length': ['log exp_mis', 'log(1+degree)', 'CEI_MM']}
    for lab, cols in models.items():
        m = logistic(g, cols)
        R.setdefault('joint logistic', {})[lab] = {'AUC': float(auc(g.Disease, m.predict())), 'OR': odds_ratios(m, cols)} # full precision
    # slope of the log-odds on ln L: 1 would mean a disease probability proportional to length (pure target size, for small probabilities)
    m1 = logistic(a, ['log L'], age=False)
    R['univariable logit slope on ln L'] = [float(m1.params['log L'])] + [float(v) for v in m1.conf_int().loc['log L']]
    mg = logistic(g, ['log L', 'log(1+degree)', 'CEI_MM'])
    R['adjusted logit slope on ln L (article covariates)'] = [float(mg.params['log L'])] + [float(v) for v in mg.conf_int().loc['log L']]
    json.dump(R, open(outfile, 'w'), indent=1, default=str) # writing the JSON output
    for k, v in R.items():
        print(k, ':', v)
    return R
