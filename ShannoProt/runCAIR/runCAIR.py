import os
import sys
if "--current" in sys.argv:  # rev 26: only the CEI of a sample of organisms in the current UniProt reference proteomes (about 400 MB of downloads)
    from Current_release_check import current_release_check
    gc_file = "GC null per proteome.csv"  # from runGC; otherwise the copy of the repository
    if not os.path.exists(gc_file):
        gc_file = "../Supplementary_materials/CAIR_supplementary_files/GC null per proteome.csv"
    current_release_check(gc_file, "reference_proteomes", "Current release check.csv", "Current release check summary.json")
    sys.exit()

import import_ipynb
from Calc_Entries_CAIR import fasta_to_entry_cair
from Calc_Proteomes_CAIR import *

fasta_to_entry_cair("uniprot_sprot.fasta", "Entries sprot.csv")  # requires the sprot.fasta input file and the output file name
fasta_to_entry_cair("uniprot_trembl.fasta", "Entries trembl.csv")  # requires the trembl.fasta input file and the output file name

entry_to_species()  # insert entries input file and the output file name

species_cair()  # insert species file and the output file name