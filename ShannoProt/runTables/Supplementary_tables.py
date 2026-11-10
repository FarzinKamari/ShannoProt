#!/usr/bin/env python
# coding: utf-8

# # Supplementary Tables S1-S4 of the article (Table S5: Table_S5.py in runNEE)
# * S1, S2 and S3: taxonomy; S4 and S5: human proteins (Tables S2 and S3 of the preprint, extended).
# * Table S1: Table S1 of the preprint (CAIR renamed CEI) with the genome GC content, the genetic code and the CEI expected under random translation of each proteome (Proteome_GC_null.py).
# * Tables S2 and S3: Brunner-Munzel tests of CEI (S2) and of the GC-adjusted CEI (S3) at every bifurcation, at proteome, species and genus level ("Bifurcation tests.R").
# * Table S4: Table S2 of the preprint with renamed headers (CAIR -> CEI, EMIP -> NEE, "Occurence" -> "Occurrence"); values unchanged.
# ##### Input files:
# >* Tables S1 and S2 of the preprint (Complete_proteome_CAIRs.csv and w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv in Supplementary_materials), GC null per proteome.csv
# >* Bifurcation tests (proteome|species|genus).csv, Bifurcation tests summary.csv and the same files "of delta CEI"
# ##### Output files:
# >* Table_S1.xlsx, Table_S2.xlsx, Table_S3.xlsx and Table_S4.xlsx (in the folder given as outdir)

# In[1]:


import os
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
              'wald_outside_01': 'Wald interval (submitted method) outside [0, 1]'}
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


def bifurcation_table(value_name='', outfile='Table_S2.xlsx', read_me=None):
    levels = {'proteome': 'proteomes', 'species': 'species', 'genus': 'genera'}
    sheets = {'Read me': pd.DataFrame({'Table %s' % ('S2' if value_name == '' else 'S3'): read_me})}
    tests = {}
    for level, unit in levels.items():
        tests[level] = pd.read_csv('Bifurcation tests%s (%s).csv' % (value_name, level))
        sheets['%s level' % level.capitalize()] = test_sheet(tests[level], unit)
    sheets['Summary'] = summary_sheet('Bifurcation tests%s summary.csv' % value_name)
    if value_name != '':
        sheets['Comparison with CEI'] = comparison_sheet(pd.read_csv('Bifurcation tests (proteome).csv'), tests['proteome'])
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
    bifurcation_table(' of delta CEI', os.path.join(outdir, 'Table_S3.xlsx'),
                      ['Brunner-Munzel tests of the GC-adjusted CEI (delta CEI = CEI - CEI expected under random translation at the genome GC content, '
                       'organism\'s genetic code) between the two descendant clades of each bifurcation, at proteome, species and genus level; '
                       '28,899 proteomes with a matched NCBI genome assembly (Table S1).'] + read_me_common +
                      ['Comparison with CEI: proteome-level results of Table S2 and of this table for the tested bifurcations.'])
    return outdir

