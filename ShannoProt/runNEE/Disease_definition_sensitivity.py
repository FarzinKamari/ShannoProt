#!/usr/bin/env python
# coding: utf-8

# # Sensitivity of the human analysis to the definition of the disease categories (rev 26)
# * Largest value instead of the sum: the occurrence category of each disease protein from the largest midpoint occurrence among its disorders instead of their sum (Kendall's tau_b, proportional-odds model and the gene-age contrast of Supplementary File S1), and the same analyses restricted to the proteins with one disorder.
# * The 1,879 reviewed proteins linked to Orphanet disorders without a usable worldwide occurrence counted as disease proteins (AUC, logistic models).
# * The reference group without the proteins that carry a UniProt "Involvement in disease" annotation.
# * Point estimates without bootstrap (the question is whether the conclusions change), with the functions of Human_indicator_statistics.py.
# ##### Input files:
# >* Human analysis data.csv (from Human_indicator_statistics.py)
# >* Table S2 of the preprint (one row per protein and disorder; CSV in Supplementary_materials, or the same table as XLSX)
# ##### Output file:
# >* Disease definition sensitivity.json

# In[1]:


import json
import numpy as np
import pandas as pd
from scipy import stats
from Human_indicator_statistics import auc, tau, logistic, ordinal, odds_ratios
bands = {'<1 / 1 000 000': 1e-6, '1-9 / 1 000 000': 5e-6, '1-9 / 100 000': 5e-5, '1-5 / 10 000': 3e-4, '6-9 / 10 000': 7.5e-4,
         '>1 / 1000': 1e-3} # representative (midpoint) values of the Orphanet occurrence bands, as in Methods section
indicators = ['NEE', 'Protein length', 'Number of interactions', 'Age rank', 'CEI', 'CEI_MM']
model_cols = ['log L', 'log(1+degree)', 'CEI_MM'] # covariates of the models of the article (with gene-age class)


# In[2]:


def disorder_values(table_s2='../Supplementary_materials/EMIP_supplementary_files/w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv'):
    pairs = pd.read_excel(table_s2) if table_s2.endswith('.xlsx') else pd.read_csv(table_s2, low_memory=False)
    pairs = pairs[pairs.Status == 'reviewed']
    pairs['value'] = pairs['Occurence value'].map(bands)
    assert pairs.value.notna().all(), pairs.loc[pairs.value.isna(), 'Occurence value'].unique()
    return pairs.groupby('Entry').agg(n_disorders=('OrphaNumber (source: Orphanet)', 'nunique'), occ_sum=('value', 'sum'),
                                      occ_max=('value', 'max')).reset_index() # number of disorders, summed and largest occurrence


def model_ors(m, cols): # odds ratio, 95% interval and p of each covariate (full precision; rounded in Table S6)
    return odds_ratios(m, cols)


def relative_effect(x, y): # Brunner-Munzel relative effect P(X < Y) + 0.5 P(X = Y) and p-value
    rk = stats.rankdata(np.r_[x, y])
    return float((rk[len(x):].mean() - (len(y) + 1) / 2) / len(x)), float(stats.brunnermunzel(x, y).pvalue)


# In[3]:


def disease_definition_sensitivity(data_file='Human analysis data.csv',
                                   table_s2='../Supplementary_materials/EMIP_supplementary_files/w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv',
                                   outfile='Disease definition sensitivity.json'): # refer to runNEE
    d = pd.read_csv(data_file, low_memory=False).merge(disorder_values(table_s2), on='Entry', how='left')
    d['n_disorders'] = d.n_disorders.fillna(0)
    R = {}
    a = d[d.Category.notna()].copy() # the 18,471 proteins of the article
    chk = a[a.Category > 0]
    R['check: sum of the disorders equals the total occurrence'] = bool(np.allclose(chk.occ_sum, chk['Total occurence value'], rtol=1e-9))
    a['Category (max)'] = np.where(a.Category == 0, 0, np.where(a.occ_max <= 1e-6, 1, 2))
    R['disorders per disease protein'] = a[a.Category > 0].n_disorders.value_counts().sort_index().astype(int).to_dict()
    R['reclassified by the maximum'] = pd.crosstab(a.Category, a['Category (max)']).to_dict()

    # Kendall's tau_b with the ordered category: sum (article), largest value, one disorder only
    def taus(x, cat):
        out = {c: float(tau(x[cat], x[c])) for c in indicators if c != 'Age rank'}
        y = x[x['Age rank'].notna()]
        out['Age rank'] = float(tau(y[cat], y['Age rank']))
        return out

    R['tau_b, sum (article)'] = taus(a, 'Category')
    R['tau_b, maximum'] = taus(a, 'Category (max)')
    one = a[(a.Category == 0) | (a.n_disorders == 1)]
    R['tau_b, proteins with one disorder only (and the reference group)'] = taus(one, 'Category')
    R['n with one disorder only'] = {int(k): int(v) for k, v in one.Category.value_counts().sort_index().items()}

    # proportional-odds model of the article (proteins with a gene age)
    g = a[a['Age rank'].notna()].copy()
    R['ordinal ORs, sum (article)'] = model_ors(ordinal(g, model_cols), model_cols)
    R['ordinal ORs, maximum'] = model_ors(ordinal(g.assign(Category=g['Category (max)']), model_cols), model_cols)
    g1 = g[(g.Category == 0) | (g.n_disorders == 1)].copy()
    R['ordinal ORs, one disorder only'] = model_ors(ordinal(g1, model_cols), model_cols)
    for lab, cat in [('sum', 'Category'), ('maximum', 'Category (max)')]: # the contrast that the sum can change: extremely rare vs rare
        dz = a[a[cat] > 0]
        R['medians extremely rare / rare, %s' % lab] = {c: [float(dz[dz[cat] == 1][c].median()), float(dz[dz[cat] == 2][c].median())]
                                                        for c in ['Protein length', 'Number of interactions', 'NEE', 'CEI_MM']}
        R['n extremely rare / rare, %s' % lab] = [int((dz[cat] == 1).sum()), int((dz[cat] == 2).sum())]
        x, y = g[g[cat] == 1]['Age rank'], g[g[cat] == 2]['Age rank'] # gene age of the extremely rare and rare categories (File S1)
        R['gene age, extremely rare vs rare (File S1), %s' % lab] = {'n': [len(x), len(y)], 'relative effect, p': relative_effect(x, y)}

    # binary analyses: the excluded proteins counted as disease proteins, and the reference group without UniProt disease annotation
    d['Disease, incl. excluded'] = np.where(d.Category.isna(), 1, (d.Category > 0).astype(float))
    # the length-degree baseline over all reviewed proteins (some numbers of interactions occur only among the excluded proteins);
    # the 18,471 proteins of the article keep the baseline of the article
    d['Surrogate'] = d.CEI.mean() * d['Protein length'] + d['Number of interactions'].map(d.groupby('Number of interactions')['Network term M'].median())
    cols_b = ['NEE', 'Surrogate', 'Protein length', 'Number of interactions', 'CEI', 'CEI_MM']

    def aucs(x, y):
        out = {c: float(auc(x[y], x[c])) for c in cols_b}
        m = logistic(x.assign(Disease=x[y]), ['log L', 'log(1+degree)'], age=False)
        out['logistic log L + log(1+degree)'] = float(auc(x[y], m.predict()))
        out['n (disease / none)'] = [int(x[y].sum()), int((x[y] == 0).sum())]
        return out

    R['AUC, article (n = 18,471)'] = aucs(a, 'Disease')
    R['AUC, excluded proteins counted as disease (n = 20,350)'] = aucs(d, 'Disease, incl. excluded')
    clean = a[~((a.Category == 0) & (a['UniProt disease'] == 1))]
    R['AUC, reference group without UniProt disease annotation'] = aucs(clean, 'Disease')
    both = d[~((d.Category == 0) & (d['UniProt disease'] == 1))]
    R['AUC, both changes'] = aucs(both, 'Disease, incl. excluded')
    for lab, x, y in [('article', a, 'Disease'), ('excluded counted as disease', d, 'Disease, incl. excluded'),
                      ('reference without UniProt disease', clean, 'Disease')]: # joint logistic model of the article
        xx = x[x['Age rank'].notna()].assign(Disease=lambda t: t[y])
        R['joint logistic ORs, %s' % lab] = model_ors(logistic(xx, model_cols), model_cols)
    json.dump(R, open(outfile, 'w'), indent=1, default=str) # writing the JSON output
    for k, v in R.items():
        print(k, ':', v)
    return R
