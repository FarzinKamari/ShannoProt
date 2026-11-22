#!/usr/bin/env python
# coding: utf-8

# # Supplementary Tables S1-S4 and S6 of the article (Table S5: Table_S5.py in runNEE)
# * S1, S2 and S3: taxonomy; S4 and S5: human proteins (Tables S2 and S3 of the preprint, extended); S6: robustness of the human analysis and CEI in the current UniProt release.
# * Table S1: Table S1 of the preprint (CAIR renamed CEI) with the genome GC content, the genetic code and the CEI expected under random translation of each proteome (Proteome_GC_null.py).
# * Tables S2 and S3: Brunner-Munzel tests of CEI (S2) and of the GC-adjusted CEI (S3) at every bifurcation, at proteome, species and genus level ("Bifurcation tests.R").
# * Table S3, added in rev 26: the tests of the CEI expected from the genome GC content alone ("Bifurcation tests.R"), the phylogenetic signal and PGLS on the GTDB trees ("Phylogenetic signal.R") and the proteomes above the highest CEI of the uniform-GC null by genome quality (GTDB_phylogeny_data.py in runGC).
# * Table S6 (rev 26): the human analysis with other definitions of the disease categories (Disease_definition_sensitivity.py), the covariates of protein length (Length_covariates.py) and CEI in the current UniProt release (Current_release_check.py in runCAIR).
# * Table S4: Table S2 of the preprint with renamed headers (CAIR -> CEI, EMIP -> NEE, "Occurence" -> "Occurrence"); values unchanged.
# ##### Input files:
# >* Tables S1 and S2 of the preprint (Complete_proteome_CAIRs.csv and w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv in Supplementary_materials), GC null per proteome.csv
# >* Bifurcation tests (proteome|species|genus).csv, Bifurcation tests summary.csv and the same files "of delta CEI"
# >* rev 26: Bifurcation tests of CEI null (proteome).csv, Phylogenetic signal.csv, Proteomes GTDB.csv, GC null maxima.csv, GC ceiling by genome quality.csv, Disease definition sensitivity.json, Length covariates.json, Current release check.csv and Current release check summary.json
# ##### Output files:
# >* Table_S1.xlsx, Table_S2.xlsx, Table_S3.xlsx, Table_S4.xlsx and Table_S6.xlsx (in the folder given as outdir)

# In[1]:


import os
import json
import numpy as np
import pandas as pd


# In[2]:


def write_sheets(outfile, sheets): # one Excel file, one sheet per table, bold headers and readable widths
    with pd.ExcelWriter(outfile, engine='openpyxl') as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)
            ws = writer.sheets[name]
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 80)
    print('written', outfile)


def table_s1(table_s1='../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv',
              gc_file='GC null per proteome.csv', outfile='Table_S1.xlsx'): # refer to runGC
    s1 = pd.read_excel(table_s1) if table_s1.endswith('.xlsx') else pd.read_csv(table_s1) # Table S1 of the preprint (Complete_proteome_CAIRs.csv, or the same table as XLSX)
    s1 = s1.rename(columns={'CAIR': 'CEI'}) # renamed header, same values
    gc = pd.read_csv(gc_file, usecols=['Organism_ID', 'GC', 'Genetic code', 'CEI null', 'Delta CEI', 'Assembly accession',
                                       'Assembly level', 'RefSeq category', 'Match type'])
    gc['Organism_ID'] = gc['Organism_ID'].astype(int)
    gc = gc.rename(columns={'GC': 'Genome GC (fraction)', 'CEI null': 'CEI expected under random translation',
                            'Delta CEI': 'Delta CEI (CEI - expected)', 'Match type': 'Assembly match'})
    out = s1.merge(gc, on='Organism_ID', how='left') # same rows and order as Table S1 of the preprint
    assert len(out) == len(s1)
    out['Assembly match'] = out['Assembly match'].fillna('no assembly with a GC value')
    write_sheets(outfile, {'Complete proteome CEIs': out})
    return out


def table_s4(table_s2='../Supplementary_materials/EMIP_supplementary_files/w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv',
             outfile='Table_S4.xlsx'): # Table S2 of the preprint is Table S4 of the article
    s2 = pd.read_excel(table_s2) if table_s2.endswith('.xlsx') else pd.read_csv(table_s2, low_memory=False)
    s2 = s2.rename(columns={'CAIR (in pits)': 'CEI', 'EMIP (in pits)': 'NEE (base-22 units)',
                            'Occurence value': 'Occurrence value', 'Occurence type': 'Occurrence type'}) # same values
    write_sheets(outfile, {'Full list of diseases': s2})
    return s2


# In[3]:


def p_value(p): # p-values below the precision of the t distribution are reported as "< 1e-300"
    return '< 1e-300' if p == 0 else float('%.3g' % p)


def test_sheet(tests, unit):
    rows = []
    for _, t in tests.iterrows():
        tested = bool(t.testable)
        rows.append({'Test': t.test, 'Upper clade (n %s)' % unit: int(t.n1), 'Lower clade (n %s)' % unit: int(t.n2),
                     'Tested': 'yes' if tested else 'no (too few %s)' % unit,
                     'RE': round(t.re, 4) if pd.notna(t.re) else None,
                     'RE 95% CI lower (logit)': round(t.logit_low, 4) if pd.notna(t.logit_low) else None,
                     'RE 95% CI upper (logit)': round(t.logit_high, 4) if pd.notna(t.logit_high) else None,
                     'Brunner-Munzel T': round(t.statistic, 3) if tested else None,
                     'df': round(t.df, 1) if tested else None,
                     'p (two-sided)': p_value(t.p) if tested else None,
                     'Adjusted p (BH)': p_value(t.p_adj) if tested else None,
                     'Permutation p': float('%.3g' % t.p_perm) if tested and pd.notna(t.p_perm) else None,
                     'Class': t['class'] if tested else None,
                     'Separated': ('yes' if t.separated else 'no') if tested else None})
    return pd.DataFrame(rows)


def summary_sheet(summary_file):
    labels = {'proteomes': 'Proteomes', 'testable': 'Testable bifurcations (proteome level)', 'significant': 'Significant after BH adjustment',
              'large': '... large effect', 'medium': '... medium effect', 'small': '... small effect', 'negligible': '... negligible effect',
              'separated': 'Separated (significant and at least small)', 'ge_medium': 'Significant with at least a medium effect',
              'nonsignificant': 'Not significant (class)', 'negligible_significant': 'Significant but negligible',
              'not_tested': 'Not tested (n upper vs n lower)', 'small_clade_tests': 'Tests with a clade of fewer than 30 proteomes',
              'permutation_same_conclusion': '... permutation test with the same conclusion at alpha = 0.05',
              'wald_outside_01': 'Wald interval outside [0, 1]'}
    for level in ('species', 'genus'):
        labels.update({level + '_units': '%s (units)' % level.capitalize(),
                       level + '_testable': '%s level: testable bifurcations' % level.capitalize(),
                       level + '_significant': '%s level: significant after BH adjustment' % level.capitalize(),
                       level + '_rho_effect': '%s level: Spearman rho of RE with the proteome level (bifurcations testable at both levels)' % level.capitalize(),
                       level + '_same_direction': '%s level: same direction of the effect' % level.capitalize(),
                       level + '_separated_testable': '%s level: bifurcations separated at proteome level and testable' % level.capitalize(),
                       level + '_still_separated': '%s level: ... still separated' % level.capitalize(),
                       level + '_lost': '%s level: ... no longer separated' % level.capitalize(),
                       level + '_median_abs_effect_proteome': '%s level: median |RE - 0.5| at proteome level (same bifurcations)' % level.capitalize(),
                       level + '_median_abs_effect_level': '%s level: median |RE - 0.5| at %s level' % (level.capitalize(), level)})
    s = pd.read_csv(summary_file, dtype=str, keep_default_na=False)
    s['item'] = s['item'].map(labels)
    s['value'] = [('%.3f' % float(v)) if v.replace('.', '', 1).isdigit() and '.' in v else v for v in s['value']]
    return s.rename(columns={'item': 'Item', 'value': 'Value'})


def comparison_sheet(cei, dcei): # proteome level: bifurcations separated on CEI and on the GC-adjusted CEI
    m = cei.merge(dcei, on='test', suffixes=('', '_d'))
    m = m[m.testable]
    out = pd.DataFrame({'Test': m.test, 'RE (CEI)': m.re.round(4), 'Separated on CEI': np.where(m.separated, 'yes', 'no'),
                        'RE (delta CEI)': m.re_d.round(4), 'Separated on delta CEI': np.where(m.separated_d, 'yes', 'no')})
    return out


# In[4]:


def bifurcation_table(value_name='', outfile='Table_S2.xlsx', read_me=None, extra_sheets=None):
    levels = {'proteome': 'proteomes', 'species': 'species', 'genus': 'genera'}
    sheets = {'Read me': pd.DataFrame({'Table %s' % ('S2' if value_name == '' else 'S3'): read_me})}
    tests = {}
    for level, unit in levels.items():
        tests[level] = pd.read_csv('Bifurcation tests%s (%s).csv' % (value_name, level))
        sheets['%s level' % level.capitalize()] = test_sheet(tests[level], unit)
    sheets['Summary'] = summary_sheet('Bifurcation tests%s summary.csv' % value_name)
    if value_name != '':
        sheets['Comparison with CEI'] = comparison_sheet(pd.read_csv('Bifurcation tests (proteome).csv'), tests['proteome'])
    sheets.update(extra_sheets or {}) # rev 26: expected CEI, phylogenetic signal and genome quality (Table S3)
    write_sheets(outfile, sheets)
    return tests


read_me_common = ['Test: bifurcation label of Figure 2A (the same comparisons as in Figure 2B and Table 2).',
                  'Upper clade / lower clade: the two descendant clades of the bifurcation (upper and lower arm in Figure 2A); n = number of proteomes, species or genera.',
                  'Tested: proteome level, both clades with more than ten proteomes; species and genus level, at least ten species (genera) in each clade. RE and its interval are also shown for untested bifurcations.',
                  'RE: relative effect of the Brunner-Munzel test, P(X < Y) + 0.5 P(X = Y), X = upper clade, Y = lower clade, i.e. the probability that an organism of the upper clade has the lower value; 0.5 = no separation.',
                  '95% CI (logit): range-preserving interval, logit(RE) +/- t(0.975, df) SE(RE) / (RE (1 - RE)), back-transformed.',
                  'p (two-sided): Brunner-Munzel test (R package brunnermunzel 2.0, t approximation with Satterthwaite df); values below 1e-300 are below the precision of the t distribution.',
                  'Adjusted p (BH): Benjamini-Hochberg adjustment over the tested bifurcations of the same level.',
                  'Permutation p: studentized permutation test (Neubert and Brunner 2007; 10,000 permutations, seed 2026) for bifurcations with a clade of fewer than 30 proteomes (proteome level only).',
                  'Class: magnitude of |RE - 0.5| by the Vargha-Delaney thresholds: < 0.06 negligible, < 0.14 small, < 0.21 medium, otherwise large. Separated: significant after adjustment and at least small.',
                  'Species and genus: species = first two words, genus = first word of the organism name (after removing "Candidatus", "uncultured", "unclassified", brackets and quotes); value = median of the proteomes; clade = most frequent clade of its proteomes.']


# In[5]:


def p_small(p): # rev 26: also the p-values below 1e-300 that are not exactly 0 (spreadsheet programs cannot store them)
    return '< 1e-300' if p < 1e-300 else float('%.3g' % p)


def expected_cei_sheet(cei, null): # proteome level: tests of the CEI expected from genome GC, with the RE of CEI (Table S2)
    out = test_sheet(null, 'proteomes')
    c = cei.set_index('test')
    out['RE of CEI (Table S2)'] = out.Test.map(c.re).round(4)
    out['Separated on CEI'] = [('yes' if c.loc[t, 'separated'] else 'no') if c.loc[t, 'testable'] else None for t in out.Test]
    both = (out.Tested == 'yes') & out['Separated on CEI'].notna()
    out['Same direction as CEI'] = np.where(both, np.where((out.RE - 0.5) * (out['RE of CEI (Table S2)'] - 0.5) > 0, 'yes', 'no'), None)
    return out


def phylogenetic_sheets(signal_file='Phylogenetic signal.csv'): # Pagel's lambda of each trait, and PGLS of CEI on GC and GC^2
    w = pd.read_csv(signal_file).pivot(index='Model', columns='Quantity', values='Value')
    traits = {'CEI': 'CEI', 'GC': 'Genome GC', 'CEI_null': 'CEI expected under random translation', 'dCEI': 'Delta CEI'}
    rows = []
    for dom in ('Bacteria', 'Archaea'):
        for t, label in traits.items():
            r = w.loc['%s %s' % (dom, t)]
            rows.append({'Domain': dom, 'Trait': label, 'GTDB species': int(r['n']), "Pagel's lambda": round(r['lambda'], 4),
                         'log-likelihood': round(r['logLik'], 1), 'log-likelihood (lambda = 0)': round(r['logLik_lambda0'], 1),
                         'LR statistic': round(r['LR'], 1), 'p (LR test, 1 df)': p_small(r['p'])})
    lam = pd.DataFrame(rows)
    rows = []
    for dom in ('Bacteria', 'Archaea'):
        r = w.loc['%s PGLS CEI ~ GC + GC2' % dom]
        rows.append({'Domain': dom, 'GTDB species': int(r['n']), "Pagel's lambda": round(r['lambda'], 4), 'Coefficient of GC': round(r['coef_GC'], 4),
                     'p (GC)': p_small(r['p_GC']), 'Coefficient of GC^2': round(r['coef_GC2'], 4), 'p (GC^2)': p_small(r['p_GC2']),
                     'R2 (PGLS)': round(r['R2'], 3), 'Adjusted R2 (PGLS)': round(r['adjR2'], 3), 'R2 (ordinary least squares)': round(r['OLS_R2'], 3)})
    return lam, pd.DataFrame(rows)


def table_s3_extra(cei_file='Bifurcation tests (proteome).csv', null_file='Bifurcation tests of CEI null (proteome).csv',
                   signal_file='Phylogenetic signal.csv', mapping_file='Proteomes GTDB.csv', maxima_file='GC null maxima.csv',
                   ceiling_file='GC ceiling by genome quality.csv'): # the sheets added to Table S3 in rev 26, and their read-me lines
    lam, pgls = phylogenetic_sheets(signal_file)
    sheets = {'Expected CEI': expected_cei_sheet(pd.read_csv(cei_file), pd.read_csv(null_file)), 'Phylogenetic signal': lam, 'PGLS': pgls,
              'Above the GC maximum': pd.read_csv(ceiling_file)}
    g = pd.read_csv(mapping_file, low_memory=False)
    pro = g[g.Domain.isin(['Bacteria', 'Archaea'])]
    match = pro['GTDB match'].fillna('not matched').value_counts()
    match = match.reindex([k for k in ('accession', 'accession (other version)', 'GenBank accession of a RefSeq genome', 'taxid', 'not matched')
                           if k in match.index])
    ceiling = float(pd.read_csv(maxima_file).iloc[0]['Highest CEI'])
    read_me = ['Expected CEI: Brunner-Munzel tests (proteome level, as above) of the CEI expected under random translation at the genome GC content (Table S1), '
               'with the RE and the separation of CEI (Table S2) and whether the two effects point in the same direction: how much of the taxonomic pattern of CEI the genome GC '
               'content alone predicts.',
               "Phylogenetic signal: Pagel's lambda (R package phylolm 2.6.5 with ape 5.8.1) of CEI, genome GC, the expected CEI and delta CEI on the species-level bacterial (bac120) and "
               'archaeal (ar53) phylogenies of GTDB release 232, one value per GTDB species representative (median of its proteomes); likelihood-ratio test against lambda fixed close '
               'to 0 (1 df). Eukaryotes are not in GTDB. Prokaryotic proteomes matched to a GTDB genome by: %s.'
               % '; '.join('%s, %s' % (k, format(int(v), ',')) for k, v in match.items()),
               "PGLS: phylogenetic generalized least squares (Pagel's lambda) of CEI on GC and GC^2, with the R2 of the same model fitted by ordinary least squares.",
               'Above the GC maximum: proteomes with a CEI above the highest CEI of the uniform-GC null (%.4f, genetic code 11), by GC assignment (Match type, from the column Assembly match of Table S1), '
               'assembly level, RefSeq category, genome quality (CheckM2 completeness and contamination of the GTDB genome) and genome origin (NCBI genome category of the GTDB '
               'genome); strictest subset: GC from an assembly of the same taxid, complete or chromosome-level, high quality and not metagenome-assembled; codon-position null: '
               'proteomes above the maximum whose genus has a codon-position null (one annotated genome per genus), and how many of them lie below it.' % ceiling]
    return sheets, read_me


def ci_text(v, per=1): # "odds ratio (lower-upper)"; per: increment of the covariate (odds ratio and interval ** per)
    return '%.2f (%.2f-%.2f)' % tuple(x ** per for x in v[:3])


def table_s6(sensitivity_file='Disease definition sensitivity.json', length_file='Length covariates.json', current_file='Current release check.csv',
             current_summary='Current release check summary.json', outfile='Table_S6.xlsx'): # robustness of the human analysis, CEI in the current release (rev 26)
    s = json.load(open(sensitivity_file, encoding='utf-8'))
    names = {'NEE': 'NEE', 'Surrogate': 'Length-degree baseline', 'Protein length': 'Protein length', 'Number of interactions': 'Number of interactions',
             'Age rank': 'Gene age', 'CEI': 'CEI', 'CEI_MM': 'CEI_MM', 'logistic log L + log(1+degree)': 'Logistic model: log length + log(1 + interactions)',
             'log L': 'log length', 'log(1+degree)': 'log(1 + interactions)'}
    per = {'CEI_MM': 0.1} # CEI_MM varies over a narrow range, so its odds ratios are given per 0.1 rather than per unit (p unchanged)
    names_or = dict(names, CEI_MM='CEI_MM (per 0.1)')
    rows = []
    for k in ['NEE', 'Protein length', 'Number of interactions', 'Age rank', 'CEI', 'CEI_MM']:
        rows.append({'Measure': "Kendall's tau_b", 'Indicator or covariate': names[k], 'Article': round(s['tau_b, sum (article)'][k], 4),
                     'Largest value': round(s['tau_b, maximum'][k], 4),
                     'One disorder only': round(s['tau_b, proteins with one disorder only (and the reference group)'][k], 4)})
    auc_cols = {'Article': 'AUC, article (n = 18,471)', 'Excluded counted as disease': 'AUC, excluded proteins counted as disease (n = 20,350)',
                'Reference without UniProt disease annotation': 'AUC, reference group without UniProt disease annotation', 'Both changes': 'AUC, both changes'}
    for k in ['NEE', 'Surrogate', 'Protein length', 'Number of interactions', 'CEI', 'CEI_MM', 'logistic log L + log(1+degree)']:
        rows.append(dict({'Measure': 'AUC', 'Indicator or covariate': names[k]}, **{c: round(s[key][k], 4) for c, key in auc_cols.items()}))
    rows.append(dict({'Measure': 'Proteins (disease / none)', 'Indicator or covariate': ''},
                     **{c: '%d / %d' % tuple(s[key]['n (disease / none)']) for c, key in auc_cols.items()}))
    for k in ['log L', 'log(1+degree)', 'CEI_MM']:
        rows.append({'Measure': 'Proportional-odds model, odds ratio (95% CI)', 'Indicator or covariate': names_or[k],
                     'Article': ci_text(s['ordinal ORs, sum (article)'][k], per.get(k, 1)),
                     'Largest value': ci_text(s['ordinal ORs, maximum'][k], per.get(k, 1)),
                     'One disorder only': ci_text(s['ordinal ORs, one disorder only'][k], per.get(k, 1))})
    for k in ['log L', 'log(1+degree)', 'CEI_MM']:
        rows.append({'Measure': 'Joint logistic model, odds ratio (95% CI)', 'Indicator or covariate': names_or[k],
                     'Article': ci_text(s['joint logistic ORs, article'][k], per.get(k, 1)),
                     'Excluded counted as disease': ci_text(s['joint logistic ORs, excluded counted as disease'][k], per.get(k, 1)),
                     'Reference without UniProt disease annotation': ci_text(s['joint logistic ORs, reference without UniProt disease'][k],
                                                                             per.get(k, 1))})
    ga = {lab: s['gene age, extremely rare vs rare (File S1), %s' % lab]['relative effect, p'] for lab in ('sum', 'maximum')}
    rows.append({'Measure': 'Gene age, extremely rare vs rare (File S1; proteins with a gene age): relative effect (p)', 'Indicator or covariate': 'Gene age',
                 'Article': '%.3f (p = %.2g)' % tuple(ga['sum']), 'Largest value': '%.3f (p = %.2g)' % tuple(ga['maximum'])})
    rows.append({'Measure': 'Proteins: extremely rare / rare', 'Indicator or covariate': '',
                 'Article': '%d / %d' % tuple(s['n extremely rare / rare, sum']), 'Largest value': '%d / %d' % tuple(s['n extremely rare / rare, maximum'])})
    rows.append({'Measure': 'Proteins moved from rare to extremely rare', 'Indicator or covariate': '', 'Largest value': s['reclassified by the maximum']['1']['2.0']})
    disease = pd.DataFrame(rows)
    L = json.load(open(length_file, encoding='utf-8'))
    labels = {'log exons': 'log number of coding exons', 'tissues TPM>=1 (per 10)': 'expression breadth (per 10 tissues)', 'tau': 'tissue specificity tau',
              'log2 max TPM': 'expression level, log2(highest median TPM + 1)', 'essential (CEGv2)': 'core essential gene (CEGv2)', 'LOEUF': 'LOEUF',
              'log exp_mis': 'log expected missense variants'}
    rows = []
    for model, r in L['joint logistic'].items():
        for k, v in r['OR'].items():
            e = per.get(k, 1)
            rows.append({'Model': model, 'AUC': round(r['AUC'], 4), 'Covariate': names_or.get(k, labels.get(k, k)), 'Odds ratio': round(v[0] ** e, 2),
                         '95% CI lower': round(v[1] ** e, 2), '95% CI upper': round(v[2] ** e, 2), 'p': float('%.2g' % v[3])})
    length = pd.DataFrame(rows)
    covariates = {'Protein length': 'protein length', 'exp_mis': 'expected missense variants', 'num_coding_exons': 'number of coding exons',
                  'gene_length': 'gene length (genomic)', 'tissues TPM>=1': 'expression breadth (tissues with median TPM >= 1)',
                  'tau': 'tissue specificity tau', 'LOEUF': 'LOEUF', 'essential (CEGv2)': 'core essential gene (CEGv2)'}
    rho = L['Spearman with protein length']
    single = pd.DataFrame([{'Covariate': covariates.get(c, c), 'Spearman rho with protein length': round(rho[c], 3) if c in rho else None,
                            'AUC alone': round(v[0], 4), 'AUC of protein length (same proteins)': round(v[1], 4), 'Proteins': v[2]}
                           for c, v in L['AUC (single covariate, same proteins)'].items()])
    cr = pd.read_csv(current_file)
    current = pd.DataFrame({'Organism': cr.Organism, 'NCBI taxid': cr.taxid, 'Domain': cr.Domain, 'Lineage': cr.Second_hierarchy,
                            'Reference proteome': cr.Proteome_ID, 'Canonical sequences': cr['proteins (current, canonical)'], 'Fragments': cr.fragments,
                            'CEI (Table S1, release 2020_02)': cr['CEI 2020'].round(6), 'CEI (current release)': cr['CEI current'].round(6),
                            'CEI (current release, without fragments)': cr['CEI current, no fragments'].round(6),
                            'Difference (current - 2020)': (cr['CEI current'] - cr['CEI 2020']).round(6)})
    C = json.load(open(current_summary, encoding='utf-8'))
    n_joint, n_dis = L['n in the joint models']
    read_me = ['Robustness of the human analysis to the definition of the disease categories, covariates of protein length, and CEI of a sample of organisms in the current '
               'UniProt release.',
               'Disease definition: Article, categories from the sum of the midpoint occurrences of the disorders of each protein (Methods section); Largest value, from the '
               'largest occurrence among its disorders; One disorder only, the rare-disease proteins with one disorder and the reference group; Excluded counted as disease, the '
               '1,879 Orphanet-linked proteins without a usable worldwide occurrence counted as rare-disease proteins; Reference without UniProt disease annotation, the reference '
               'proteins with a UniProt "Involvement in disease" annotation removed; Both changes, the two previous changes together. Point estimates (no bootstrap).',
               "Kendall's tau_b: with the ordered categories (none < extremely rare < rare); AUC: any rare disease vs none; length-degree baseline: the baseline of NEE built "
               'from length and the number of interactions only. Odds ratios per unit of log length and log(1 + interactions) and per 0.1 of CEI_MM, with gene-age '
               'class, proteins with a gene age. Gene age: relative effect of the Brunner-Munzel test, P(X < Y) + 0.5 P(X = Y), X = extremely rare, Y = rare.',
               'Length covariates: joint logistic models of rare-disease involvement (any rare disease vs none) with the covariates of the article (log length, log(1 + '
               'interactions), CEI_MM, gene-age class) and the listed covariates, %s proteins (%s rare-disease proteins) with all covariates; AUC of the model. Covariates: number '
               'of coding exons, expected missense variants and LOEUF of the canonical transcript (gnomAD v2.1.1 constraint table); expression breadth (tissues with a median '
               'TPM of at least 1), tissue specificity tau (on log2(TPM + 1)) and the highest median TPM (GTEx v8); core essential genes CEGv2 (Hart et al. 2017). Genes '
               'matched through the Ensembl gene and GeneID cross-references of UniProt release 2020_02.' % (format(n_joint, ','), format(n_dis, ',')),
               'Length correlates: Spearman correlation with protein length, AUC of the covariate alone (lower LOEUF and lower tau taken as higher values) and AUC of protein '
               'length in the same proteins.',
               'Current release: CEI of a sample of organisms of Table S1 (up to four per second-rank lineage, random with seed 2026, and the organisms of Table 1) recomputed '
               'from the canonical sequences (one per gene) of the UniProt reference proteome of the same NCBI taxid, %s, with and without the sequences flagged as fragments; '
               'counted as in Table S1 (22 symbols; B, Z and J split; X not counted). Spearman correlation with Table S1: %.3f (without fragments, %.3f); 95th percentile of the '
               'absolute difference: %s; largest: %s.'
               % (C['UniProt release'].replace('Release ', 'release ').split(',')[0], C['CEI current']['Spearman with 2020'], C['CEI current, no fragments']['Spearman with 2020'],
                  C['CEI current']['95th percentile |difference|'], C['CEI current']['max |difference|'])]
    write_sheets(outfile, {'Read me': pd.DataFrame({'Table S6': read_me}), 'Disease definition': disease, 'Length covariates': length,
                           'Length correlates': single, 'Current release': current})
    return outfile


# In[6]:


def supplementary_tables(table_s1_file='../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv',
                         table_s2_file='../Supplementary_materials/EMIP_supplementary_files/w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv',
                         outdir='Supplementary tables'):
    os.makedirs(outdir, exist_ok=True)
    table_s1(table_s1_file, outfile=os.path.join(outdir, 'Table_S1.xlsx'))
    table_s4(table_s2_file, outfile=os.path.join(outdir, 'Table_S4.xlsx'))
    bifurcation_table('', os.path.join(outdir, 'Table_S2.xlsx'),
                      ['Brunner-Munzel tests of CEI between the two descendant clades of each bifurcation of the tree of life, '
                       'at proteome, species and genus level.'] + read_me_common)
    extra_sheets, extra_read_me = table_s3_extra() # rev 26
    bifurcation_table(' of delta CEI', os.path.join(outdir, 'Table_S3.xlsx'),
                      ['Brunner-Munzel tests of the GC-adjusted CEI (delta CEI = CEI - CEI expected under random translation at the genome GC content, '
                       'organism\'s genetic code) between the two descendant clades of each bifurcation, at proteome, species and genus level; '
                       '28,899 proteomes with a matched NCBI genome assembly (Table S1).'] + read_me_common +
                      ['Comparison with CEI: proteome-level results of Table S2 and of this table for the tested bifurcations.'] + extra_read_me, extra_sheets)
    table_s6(outfile=os.path.join(outdir, 'Table_S6.xlsx')) # rev 26
    return outdir

