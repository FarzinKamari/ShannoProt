#!/usr/bin/env python
# coding: utf-8

# # Quantitative comparison of the human disease indicators -- NEE (formerly EMIP), length, interactions, gene age, CEI
# * This module compares the indicators of rare-disease involvement as described in Methods section:
# >* AUC (any rare disease vs no recorded rare disease) and Kendall's tau_b (ordered categories: none < extremely rare < rare), with 95% percentile intervals from 1,000 bootstrap resamples of proteins
# >* paired bootstrap differences of AUC; the length-degree baseline of NEE ("Surrogate" in the code), built from length and the number of interactions only
# >* logistic and proportional-odds models with log length, log(1 + interactions), gene age, CEI_MM and ln NEE (likelihood-ratio tests), with HuRI degree and PubMed counts (all articles, and focused articles on at most 10 or 100 human genes) as study-bias covariates
# >* gene age of the extremely rare and the rare category (Brunner-Munzel test; Supplementary File S1)
# >* Miller-Madow correction of CEI, degree deciles and the UniProt "Involvement in disease" labels as sensitivity analyses
# >* whether NEE adds anything to the length-degree baseline, and the association of the compositional divergence D_JS
# * Disease categories and gene ages are those of Table S3 of the preprint (Orphanet data of May 2020 and Liebeskind et al. 2016).
# ##### Input files:
# >* Table S3 of the preprint (CSV in Supplementary_materials), NEE decomposition.csv, Study bias covariates.csv, NEE null models.csv
# ##### Output files:
# >* Human analysis data.csv (one row per analysed protein) and Human indicator statistics.json

# In[1]:


import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel
ages = ['Mammalia', 'Vertebrata', 'Eumetazoa', 'Opisthokonta', 'Eukaryota', 'Euk_Archaea', 'Euk+Bac', 'Cellular_organisms'] # youngest to oldest
n_boot = 1000 # bootstrap resamples of proteins
rng = np.random.default_rng(2026) # random numbers of the bootstrap


# In[2]:


def prepare_data(table_s3='../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv',
                 nee_file='NEE decomposition.csv', bias_file='Study bias covariates.csv',
                 outfile='Human analysis data.csv'): # one row per reviewed protein with all indicators
    s3 = pd.read_excel(table_s3) if table_s3.endswith('.xlsx') else pd.read_csv(table_s3, low_memory=False) # Table S3 of the preprint (CSV, or the same table as XLSX)
    s3 = s3[s3['Status'] == 'reviewed'][['Entry', 'Total occurence value', 'Involvement in disease (source: UniProt)', 'Gene age']]
    d = s3.merge(pd.read_csv(nee_file), on='Entry').merge(pd.read_csv(bias_file), on='Entry') # adding the recalculated indicators
    occ = d['Total occurence value'] # NaN: Orphanet link without a usable worldwide occurrence (excluded)
    d['Category'] = np.where(occ.isna(), np.nan, np.where(occ == 0, 0, np.where(occ <= 1e-6, 1, 2))) # 0 none, 1 extremely rare, 2 rare
    d['Disease'] = (d.Category > 0).astype(float).where(d.Category.notna()) # any rare disease vs none
    d['Age rank'] = d['Gene age'].map({a: i for i, a in enumerate(ages)}) # 0 = Mammalia ... 7 = Cellular organisms
    d['UniProt disease'] = d['Involvement in disease (source: UniProt)'].notna().astype(int)
    d['log L'] = np.log(d['Protein length'])
    d['log(1+degree)'] = np.log1p(d['Number of interactions'])
    d['ln NEE'] = np.log(d.NEE)
    d['log L_N'] = np.log(d.L_N)
    d['log(1+HuRI degree)'] = np.log1p(d['HuRI degree'])
    d['log(1+PubMed to 2019)'] = np.log1p(d['PubMed articles to 2019'])
    d['log(1+PubMed all)'] = np.log1p(d['PubMed articles (all)'])
    d['log(1+PubMed focused to 2019)'] = np.log1p(d['PubMed focused articles to 2019']) # articles on at most 10 human genes
    d['log(1+PubMed to 2019, 100 genes)'] = np.log1p(d['PubMed articles to 2019 on at most 100 genes'])
    a = d[d.Category.notna()] # analysed proteins
    d['Surrogate'] = a.CEI.mean() * d['Protein length'] + d['Number of interactions'].map(a.groupby('Number of interactions')['Network term M'].median())
    # (surrogate of NEE from length and interactions only: mean CEI x length + median network term of proteins with the same number of interactions)
    d.to_csv(outfile, index=False) # writing the CSV output
    return d


# In[3]:


def auc(y, x): # AUC of an indicator (higher values = disease)
    return roc_auc_score(y, x)


def tau(c, x): # Kendall's tau_b between an indicator and the ordered categories
    return stats.kendalltau(x, c, variant='b')[0]


def boot(f, n): # 95% percentile interval of a statistic f(index) over bootstrap resamples of proteins
    vals = [f(rng.integers(0, n, n)) for _ in range(n_boot)]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def univariable(a, col): # AUC and tau_b of one indicator with bootstrap intervals
    y, c, x = a.Disease.values, a.Category.values, a[col].values
    return {'n': len(a), 'AUC': auc(y, x), 'AUC CI': boot(lambda i: auc(y[i], x[i]), len(a)),
            'tau_b': tau(c, x), 'tau_b CI': boot(lambda i: tau(c[i], x[i]), len(a))}


def model_scores(a, cols, i=None, age=True): # in-sample predicted probabilities of a logistic model (refitted on resamples)
    s = a if i is None else a.iloc[i]
    X = s[cols].copy()
    if age:
        X = pd.concat([X, pd.get_dummies(s['Gene age'], prefix='age', drop_first=True).astype(float)], axis=1) # gene age as categories
    m = sm.Logit(s.Disease.values, sm.add_constant(X)).fit(disp=0)
    return s.Disease.values, s.Category.values, m.predict(sm.add_constant(X))


def model_univariable(a, cols, age=True): # AUC and tau_b of a model's predictions with bootstrap intervals (models refitted)
    y, c, p = model_scores(a, cols, age=age)
    def f_auc(i):
        y, c, p = model_scores(a, cols, i, age); return auc(y, p)
    def f_tau(i):
        y, c, p = model_scores(a, cols, i, age); return tau(c, p)
    return {'n': len(a), 'AUC': auc(y, p), 'AUC CI': boot(f_auc, len(a)), 'tau_b': tau(c, p), 'tau_b CI': boot(f_tau, len(a))}


def logistic(a, cols, age=True): # binary logistic model
    X = a[cols].copy()
    if age:
        X = pd.concat([X, pd.get_dummies(a['Gene age'], prefix='age', drop_first=True).astype(float)], axis=1)
    return sm.Logit(a.Disease.values, sm.add_constant(X)).fit(disp=0)


def ordinal(a, cols, age=True): # proportional-odds model of the ordered categories
    X = a[cols].copy()
    if age:
        X = pd.concat([X, pd.get_dummies(a['Gene age'], prefix='age', drop_first=True).astype(float)], axis=1)
    return OrderedModel(pd.Categorical(a.Category, ordered=True), X, distr='logit').fit(method='bfgs', disp=0, maxiter=5000)


def odds_ratios(m, cols): # odds ratios with 95% Wald intervals
    ci = np.exp(m.conf_int())
    return {c: [float(np.exp(m.params[c])), float(ci.loc[c].iloc[0]), float(ci.loc[c].iloc[1]), float(m.pvalues[c])] for c in cols}


def lr_test(m0, m1): # likelihood-ratio test of two models that differ by one variable
    lr = 2 * (m1.llf - m0.llf)
    return {'LR': float(lr), 'df': 1, 'p': float(stats.chi2.sf(lr, 1))}



# In[4]:


def human_indicator_statistics(table_s3='../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv',
                               nee_file='NEE decomposition.csv', bias_file='Study bias covariates.csv',
                               nulls_file='NEE null models.csv', data_outfile='Human analysis data.csv',
                               outfile='Human indicator statistics.json'): # refer to runNEE
    d = prepare_data(table_s3, nee_file, bias_file, data_outfile)
    R = {} # results
    a = d[d.Category.notna()].reset_index(drop=True) # analysed proteins (Table 3)
    g = a[a['Age rank'].notna()].reset_index(drop=True) # analysed proteins with a gene age
    excl = d[d.Category.isna()] # Orphanet link without a usable worldwide occurrence

    ## Table 3 and the exclusions
    R['table3'] = {'reviewed': len(d), 'excluded (no usable occurrence)': len(excl),
                   'analysed': len(a), 'groups (none, extremely rare, rare)': a.Category.value_counts().sort_index().astype(int).tolist(),
                   'with gene age': len(g), 'with gene age groups': g.Category.value_counts().sort_index().astype(int).tolist(),
                   'reviewed with gene age': int(d['Age rank'].notna().sum()), 'gene age but excluded': int(excl['Age rank'].notna().sum()),
                   'reference group with UniProt disease': int(a[a.Category == 0]['UniProt disease'].sum()),
                   'reference group size': int((a.Category == 0).sum())}

    ## what NEE is made of: length, connectivity and composition
    sp = lambda x, y: float(stats.spearmanr(a[x], a[y])[0])
    R['spearman'] = {'NEE-L': sp('NEE', 'Protein length'), 'NEE-degree': sp('NEE', 'Number of interactions'),
                     'L-degree': sp('Protein length', 'Number of interactions'), 'M-degree': sp('Network term M', 'Number of interactions'),
                     'NEE-L_N': sp('NEE', 'L_N'), 'NEE-surrogate': sp('NEE', 'Surrogate'), 'NEE-NEE_MM': sp('NEE', 'NEE_MM'),
                     'NEE-NEE on HuRI': sp('NEE', 'NEE on HuRI'), 'degree-HuRI degree': sp('Number of interactions', 'HuRI degree'),
                     'degree-PubMed to 2019': sp('Number of interactions', 'PubMed articles to 2019'),
                     'HuRI degree-PubMed to 2019': sp('HuRI degree', 'PubMed articles to 2019'),
                     'degree-PubMed focused to 2019': sp('Number of interactions', 'PubMed focused articles to 2019'),
                     'median share of own term': float(np.median(a['Own term (L*CEI)'] / a.NEE)),
                     'median D_JS (proteins with partners)': float(a[a['Partners in proteome'] > 0].D_JS.median())}

    ## univariable comparison (Table 4)
    uni = {}
    for col in ['NEE', 'Number of interactions', 'Protein length', 'CEI', 'CEI_MM', 'NEE_MM', 'L_N', 'Surrogate',
                'HuRI degree', 'NEE on HuRI', 'PubMed articles to 2019']:
        uni[col] = univariable(a, col)
    uni['Gene age'] = univariable(g, 'Age rank')
    uni['Model: log L + log(1+degree)'] = model_univariable(a, ['log L', 'log(1+degree)'], age=False)
    uni['Model: full'] = model_univariable(g, ['log L', 'log(1+degree)', 'CEI_MM', 'ln NEE']) # log L, log(1+degree), gene age, CEI_MM, ln NEE
    R['univariable'] = uni
    y, n = a.Disease.values, len(a)
    def d_auc(c1, c2): # paired bootstrap difference of AUC
        x1, x2 = a[c1].values, a[c2].values
        return {'estimate': auc(y, x1) - auc(y, x2), 'CI': boot(lambda i: auc(y[i], x1[i]) - auc(y[i], x2[i]), n)}
    R['paired dAUC'] = {'NEE - length': d_auc('NEE', 'Protein length'), 'NEE - degree': d_auc('NEE', 'Number of interactions'),
                        'NEE - surrogate': d_auc('NEE', 'Surrogate'), 'NEE - L_N': d_auc('NEE', 'L_N')}

    ## multivariable models (proteins with a gene age)
    base = ['log L', 'log(1+degree)']
    def pred(m, cols, age=True): # in-sample predicted probabilities
        X = pd.concat([g[cols]] + ([pd.get_dummies(g['Gene age'], prefix='age', drop_first=True).astype(float)] if age else []), axis=1)
        return m.predict(sm.add_constant(X))
    m_ld, m_lda = logistic(g, base, age=False), logistic(g, base)
    m_ldac, m_full = logistic(g, base + ['CEI_MM']), logistic(g, base + ['CEI_MM', 'ln NEE'])
    R['multivariable'] = {'n': len(g), 'AUC L+degree': auc(g.Disease, pred(m_ld, base, False)), 'AUC +age': auc(g.Disease, pred(m_lda, base)),
                          'AUC +CEI_MM': auc(g.Disease, pred(m_ldac, base + ['CEI_MM'])),
                          'AUC +ln NEE': auc(g.Disease, pred(m_full, base + ['CEI_MM', 'ln NEE'])),
                          'LR ln NEE': lr_test(m_ldac, m_full), 'OR': odds_ratios(m_full, base + ['CEI_MM', 'ln NEE'])}
    o0, o1 = ordinal(g, base + ['CEI_MM']), ordinal(g, base + ['CEI_MM', 'ln NEE'])
    R['ordinal'] = {'n': len(g), 'LR ln NEE': lr_test(o0, o1), 'OR': odds_ratios(o1, base + ['CEI_MM', 'ln NEE'])}
    if (a.NEE_MM > 0).all(): # NEE_MM as a sensitivity analysis
        g['ln NEE_MM'] = np.log(g.NEE_MM)
        mm = logistic(g, base + ['CEI_MM', 'ln NEE_MM'])
        R['multivariable NEE_MM'] = {'LR ln NEE_MM': lr_test(m_ldac, mm), 'OR': odds_ratios(mm, ['ln NEE_MM'])}
    else:
        R['multivariable NEE_MM'] = {'proteins with NEE_MM <= 0': int((a.NEE_MM <= 0).sum())}

    ## study bias: PubMed counts and HuRI
    p = g[g['GeneID mapped']].reset_index(drop=True) # proteins with a GeneID cross-reference
    for label, pub in [('PubMed to 2019', 'log(1+PubMed to 2019)'), ('PubMed all', 'log(1+PubMed all)'),
                       ('PubMed focused to 2019', 'log(1+PubMed focused to 2019)'), ('PubMed to 2019, 100 genes', 'log(1+PubMed to 2019, 100 genes)')]:
        cols = base + ['CEI_MM', pub]
        b0, b1 = logistic(p, cols), logistic(p, cols + ['ln NEE'])
        q0, q1 = ordinal(p, cols), ordinal(p, cols + ['ln NEE'])
        R['study bias ' + label] = {'n': len(p), 'logistic LR ln NEE': lr_test(b0, b1), 'logistic OR': odds_ratios(b1, cols + ['ln NEE']),
                                    'ordinal LR ln NEE': lr_test(q0, q1), 'ordinal OR': odds_ratios(q1, cols + ['ln NEE'])}
    g['ln NEE on HuRI'] = np.log(g['NEE on HuRI'])
    huri = ['log L', 'log(1+HuRI degree)', 'CEI_MM']
    h0, h1 = logistic(g, huri), logistic(g, huri + ['ln NEE on HuRI'])
    hq0, hq1 = ordinal(g, huri), ordinal(g, huri + ['ln NEE on HuRI'])
    hb = logistic(g, ['log L', 'log(1+degree)', 'log(1+HuRI degree)', 'CEI_MM'])
    hh = a[a['HuRI degree'] > 0] # proteins with at least one HuRI partner
    R['HuRI'] = {'reviewed with a HuRI partner': int((d['HuRI degree'] > 0).sum()),
                 'analysed with a HuRI partner': len(hh), 'by group (none, extremely rare, rare)': hh.Category.value_counts().sort_index().astype(int).tolist(),
                 'analysed with an Ensembl gene': int(a['Ensembl gene mapped'].sum()),
                 'AUC HuRI degree (with a HuRI partner)': auc(hh.Disease, hh['HuRI degree']),
                 'AUC NEE on HuRI (with a HuRI partner)': auc(hh.Disease, hh['NEE on HuRI']),
                 'AUC NEE (with a HuRI partner)': auc(hh.Disease, hh.NEE),
                 'AUC length (with a HuRI partner)': auc(hh.Disease, hh['Protein length']),
                 'AUC degree (with a HuRI partner)': auc(hh.Disease, hh['Number of interactions']),
                 'logistic LR ln NEE on HuRI': lr_test(h0, h1), 'logistic OR': odds_ratios(h1, huri + ['ln NEE on HuRI']),
                 'ordinal LR ln NEE on HuRI': lr_test(hq0, hq1), 'ordinal OR': odds_ratios(hq1, huri + ['ln NEE on HuRI']),
                 'logistic OR with both degrees': odds_ratios(hb, ['log L', 'log(1+degree)', 'log(1+HuRI degree)', 'CEI_MM'])}

    ## the network term beyond the number of interactions
    r = a[a['Number of interactions'] > 0].copy()
    r['log M'] = np.log(r['Network term M'].clip(lower = 1e-6))
    fd = sm.OLS(r['log M'], sm.add_constant(pd.DataFrame({'a': r['log(1+degree)'], 'b': r['log(1+degree)'] ** 2}))).fit()
    r['residual'] = r['log M'] - fd.fittedvalues # network term not explained by the number of interactions
    mr = sm.Logit(r.Disease.values, sm.add_constant(r[['log L', 'log(1+degree)', 'residual']])).fit(disp=0)
    R['residual network term'] = {'n': len(r), 'R2 on degree': float(fd.rsquared),
                                  'OR per SD': float(np.exp(mr.params['residual'] * r.residual.std())), 'p': float(mr.pvalues['residual'])}
    a['degree decile'] = pd.qcut(a['Number of interactions'].rank(method='first'), 10, labels=False)
    R['degree deciles'] = [{'degree': '%d-%d' % (x['Number of interactions'].min(), x['Number of interactions'].max()), 'n': len(x),
                            'AUC NEE': auc(x.Disease, x.NEE), 'AUC length': auc(x.Disease, x['Protein length'])}
                           for _, x in a.groupby('degree decile') if x.Disease.nunique() == 2]

    ## UniProt "Involvement in disease" labels (all reviewed proteins)
    R['UniProt labels'] = {'n': len(d), 'positive': int(d['UniProt disease'].sum())}
    d['Surrogate (all reviewed)'] = d.CEI.mean() * d['Protein length'] + d['Number of interactions'].map(d.groupby('Number of interactions')['Network term M'].median())
    # (the surrogate rebuilt over all reviewed proteins, since some numbers of interactions occur only among the excluded ones)
    for c in ['NEE', 'Surrogate (all reviewed)', 'Number of interactions', 'Protein length', 'CEI', 'CEI_MM']:
        R['UniProt labels'][c] = auc(d['UniProt disease'], d[c])

    ## Miller-Madow correction of CEI (exact number of residue types observed)
    sg = lambda x, c: float(stats.spearmanr(x['Age rank'], x[c])[0])
    long = g[g['Protein length'] >= 300]
    R['Miller-Madow'] = {'rho age-CEI': sg(g, 'CEI'), 'rho age-CEI_MM': sg(g, 'CEI_MM'),
                         'rho age-CEI (L>=300)': sg(long, 'CEI'), 'rho age-CEI_MM (L>=300)': sg(long, 'CEI_MM'),
                         'mean correction (L 95-105)': float((a.CEI_MM - a.CEI)[a['Protein length'].between(95, 105)].mean()),
                         'mean correction (L 285-315)': float((a.CEI_MM - a.CEI)[a['Protein length'].between(285, 315)].mean()),
                         'median by age (CEI, CEI_MM, L)': {k: [float(x.CEI.median()), float(x.CEI_MM.median()), float(x['Protein length'].median())]
                                                            for k, x in g.groupby('Gene age')},
                         'median CEI_MM by group (none, extremely rare, rare)': a.groupby('Category').CEI_MM.median().tolist()}

    ## null models of NEE (composition shuffle and degree-preserving rewiring)
    nulls = pd.read_csv(nulls_file)
    R['null models'] = {'observed AUC NEE': auc(y, a.NEE), 'observed tau_b NEE': tau(a.Category.values, a.NEE.values)}
    for kind, x in nulls[nulls['Null model'] != 'observed'].groupby('Null model'): # (the observed row is summarized above)
        R['null models'][kind] = {'replicates': len(x), 'AUC mean': float(x.AUC.mean()), 'AUC SD': float(x.AUC.std()),
                                  'AUC min': float(x.AUC.min()), 'AUC max': float(x.AUC.max()),
                                  'tau_b mean': float(x.tau_b.mean()), 'tau_b SD': float(x.tau_b.std()),
                                  'share of AUC >= observed': float((x.AUC >= R['null models']['observed AUC NEE']).mean())}

    ## study-bias indicators among the proteins that have the information (Section 3.7)
    # (placed last, so that the bootstrap intervals above do not change)
    hh = a[a['HuRI degree'] > 0] # proteins with at least one HuRI partner
    pm = a[a['GeneID mapped']] # proteins with a GeneID cross-reference (PubMed counts)
    R['study bias univariable'] = {'with a HuRI partner': {c: univariable(hh, c) for c in ['HuRI degree', 'NEE on HuRI', 'NEE',
                                                                                          'Protein length', 'Number of interactions']},
                                   'with a GeneID': {c: univariable(pm, c) for c in ['PubMed articles to 2019', 'NEE', 'Protein length',
                                                                                    'Number of interactions']}}
    R['study bias univariable']['with a GeneID']['PubMed focused articles to 2019'] = univariable(pm, 'PubMed focused articles to 2019')
    # (last bootstrap of the module, so that all intervals above stay the same)

    ## gene age of the extremely rare and the rare category (Supplementary File S1; no random numbers)
    xr, ra = g['Age rank'][g.Category == 1], g['Age rank'][g.Category == 2]
    rk = stats.rankdata(np.concatenate([xr, ra])) # mid-ranks of the pooled gene-age classes
    bm = stats.brunnermunzel(xr, ra)
    R['gene age extremely rare vs rare'] = {'n': [len(xr), len(ra)],
                                            'median class': [ages[int(np.median(xr))], ages[int(np.median(ra))]],
                                            'share Euk_Archaea or older': [float((xr >= 5).mean()), float((ra >= 5).mean())],
                                            'relative effect': float((rk[len(xr):].mean() - (len(ra) + 1) / 2) / len(xr)),
                                            'Brunner-Munzel statistic': float(bm.statistic), 'p': float(bm.pvalue)}
    # (relative effect: P(X < Y) + P(X = Y)/2 for the gene-age classes of an extremely rare (X) and a rare (Y) disease protein;
    # below 0.5 when the extremely rare disease proteins are older)

    ## NEE beyond the length-degree baseline, and the compositional divergence D_JS (no random numbers, so that all intervals above stay the same)
    # does NEE add anything to the length-degree baseline and to the length + interactions model?
    a['ln Surrogate'] = np.log(a.Surrogate)
    s0, s1 = logistic(a, ['ln Surrogate'], age=False), logistic(a, ['ln Surrogate', 'ln NEE'], age=False)
    b0, b1 = logistic(a, base, age=False), logistic(a, base + ['ln NEE'], age=False)
    R['NEE beyond length and interactions'] = {'n': len(a),
                                               'baseline: AUC': auc(y, s0.predict()), 'baseline + ln NEE: AUC': auc(y, s1.predict()),
                                               'baseline + ln NEE: LR': lr_test(s0, s1),
                                               'log L + log(1+degree): AUC': auc(y, b0.predict()), '+ ln NEE: AUC': auc(y, b1.predict()),
                                               '+ ln NEE: LR': lr_test(b0, b1), '+ ln NEE: OR': odds_ratios(b1, ['ln NEE'])}
    # the compositional divergence D_JS of the neighborhood (proteins with at least one partner in the proteome)
    w = a[a['Partners in proteome'] > 0].reset_index(drop=True)
    w['D_JS per SD'] = (w.D_JS - w.D_JS.mean()) / w.D_JS.std() # standardized, so that the odds ratio is per standard deviation
    j0, j1 = logistic(w, base, age=False), logistic(w, base + ['D_JS per SD'], age=False)
    wg = w[w['Age rank'].notna()].reset_index(drop=True)
    k0, k1 = logistic(wg, base + ['CEI_MM']), logistic(wg, base + ['CEI_MM', 'D_JS per SD'])
    R['D_JS'] = {'n with partners': len(w), 'median': float(w.D_JS.median()), 'IQR': [float(w.D_JS.quantile(0.25)), float(w.D_JS.quantile(0.75))],
                 'rho D_JS-degree': float(stats.spearmanr(w.D_JS, w['Number of interactions'])[0]),
                 'AUC D_JS': auc(w.Disease, w.D_JS),
                 'L+degree: AUC': auc(w.Disease, j0.predict()), 'L+degree+D_JS: AUC': auc(w.Disease, j1.predict()),
                 'L+degree+D_JS: LR': lr_test(j0, j1), 'L+degree+D_JS: OR per SD': odds_ratios(j1, ['D_JS per SD']),
                 'n with gene age': len(wg), 'joint: AUC': auc(wg.Disease, k0.predict()), 'joint+D_JS: AUC': auc(wg.Disease, k1.predict()),
                 'joint+D_JS: LR': lr_test(k0, k1), 'joint+D_JS: OR per SD': odds_ratios(k1, ['D_JS per SD'])}
    json.dump(R, open(outfile, 'w'), indent=1) # writing the JSON output
    return R
