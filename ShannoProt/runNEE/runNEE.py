import sys
from Human_proteome_2020_02 import dat_to_proteome_table
from NEE_decomposition import nee_decomposition, check_against_table_s3
from HuRI_PubMed_covariates import study_bias_covariates
from Human_indicator_statistics import prepare_data, human_indicator_statistics
from NEE_nulls import nee_null_models
from Table_S5 import table_s5

# Table S3 of the preprint (Supplementary_materials/EMIP_supplementary_files), or the same table as XLSX
table_s3 = sys.argv[1] if len(sys.argv) > 1 else "../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv"

dat_to_proteome_table("uniprot_sprot.dat.gz", "Human proteome 2020_02.tab")  # requires the Swiss-Prot DAT file of UniProt release 2020_02
nee_decomposition("Human proteome 2020_02.tab", "UniProtNormalizedTabular-default.txt", "NEE decomposition.csv")  # requires the PICKLE 2.5 file
check_against_table_s3("NEE decomposition.csv", table_s3)  # the recalculated NEE (formerly EMIP) and interactions must reproduce Table S3 of the preprint

study_bias_covariates("NEE decomposition.csv", "Human proteome 2020_02.tab", "HuRI.tsv", "gene2pubmed.gz",
                      "Study bias covariates.csv")  # requires the HuRI and gene2pubmed files
prepare_data(table_s3, "NEE decomposition.csv", "Study bias covariates.csv", "Human analysis data.csv")  # the data of the null models
nee_null_models("Human proteome 2020_02.tab", "UniProtNormalizedTabular-default.txt", "Human analysis data.csv",
                "NEE null models.csv")  # 1,000 composition shuffles and 100 rewired networks (about half an hour)
human_indicator_statistics(table_s3, "NEE decomposition.csv", "Study bias covariates.csv", "NEE null models.csv",
                           "Human analysis data.csv", "Human indicator statistics.json")  # Tables 3 and 4, Sections 3.4-3.7 (about 10 min)

table_s5(table_s3, "NEE decomposition.csv", "Study bias covariates.csv", "Table_S5.csv",
         "Table_S5.xlsx")  # Table S5 of the article (input of "DTK test & Errorbars.R")
