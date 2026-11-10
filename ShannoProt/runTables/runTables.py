import os
import shutil
from Supplementary_tables import supplementary_tables

# Tables S1-S4 of the article (run after runNEE.py, runGC.py and "Bifurcation tests.R"): S2 and S3 are the bifurcation tests,
# S4 and S5 the human proteins (Tables S2 and S3 of the preprint, extended)
supplementary_tables("../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv",
                     "../Supplementary_materials/EMIP_supplementary_files/w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv",
                     "Supplementary tables")
if os.path.exists("Table_S5.xlsx"):
    shutil.copy("Table_S5.xlsx", os.path.join("Supplementary tables", "Table_S5.xlsx"))  # Table S5 (from runNEE.py)
