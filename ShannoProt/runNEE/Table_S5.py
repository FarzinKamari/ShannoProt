#!/usr/bin/env python
# coding: utf-8

# # Table S5 -- human proteome with the renamed indicators and the new columns (Table S3 of the preprint, extended)
# * The columns of Table S3 of the preprint are kept with their values; only the headers are renamed (CAIR -> CEI, EMIP -> NEE, "Occurence" -> "Occurrence").
# * New columns (reviewed entries; double precision, see NEE_decomposition.py): disease category, CEI_MM, NEE_MM, L_N, D_JS, HuRI degree, PubMed articles up to 2019 and focused PubMed articles up to 2019 (articles on at most 10 human genes).
# * The CEI and NEE values of 2020 were computed in single precision; they agree with the double-precision values to within 4e-5 (relative), so NEE = L*CEI + L_N*D_JS holds to that precision.
# ##### Input files:
# >* Table S3 of the preprint (CSV in Supplementary_materials), NEE decomposition.csv, Study bias covariates.csv
# ##### Output files:
# >* Table_S5.csv and Table_S5.xlsx

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


def table_s5(table_s3='../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv',
             nee_file='NEE decomposition.csv', bias_file='Study bias covariates.csv',
             outfile='Table_S5.csv', xlsx_outfile='Table_S5.xlsx'): # refer to runNEE
    s3 = pd.read_excel(table_s3) if table_s3.endswith('.xlsx') else pd.read_csv(table_s3, low_memory=False) # Table S3 of the preprint (CSV, or the same table as XLSX); reviewed and unreviewed entries
    s3 = s3.rename(columns={'CAIR (in pits)': 'CEI',
                            'EMIP (in pits)': 'NEE (base-22 units)',
                            'Total occurence value': 'Total occurrence value'}) # renamed headers, same values
    occ = s3['Total occurrence value']
    category = np.where(occ.isna(), 'Excluded (no usable worldwide occurrence)',
                        np.where(occ == 0, 'No recorded rare disease', np.where(occ <= 1e-6, 'Extremely rare', 'Rare')))
    s3['Disease category'] = np.where(s3['Status'] == 'reviewed', category, '') # categories of the analysis (reviewed entries)
    nee = pd.read_csv(nee_file)[['Entry', 'CEI_MM', 'NEE_MM', 'L_N', 'D_JS']].rename(columns={'NEE_MM': 'NEE_MM (base-22 units)',
                                                                                         'L_N': 'L_N (residues)'})
    bias = pd.read_csv(bias_file)[['Entry', 'HuRI degree', 'PubMed articles to 2019', 'PubMed focused articles to 2019']]
    out = s3.merge(nee, on='Entry', how='left').merge(bias, on='Entry', how='left') # new columns (reviewed entries only)
    out.to_csv(outfile, index=False) # writing the CSV output
    if xlsx_outfile:
        out.to_excel(xlsx_outfile, index=False, sheet_name='Full list of human proteins') # writing the Excel output (Supplementary Materials)
    print(len(out), 'rows;', int(out.CEI_MM.notna().sum()), 'reviewed entries with the new columns')
    return out

