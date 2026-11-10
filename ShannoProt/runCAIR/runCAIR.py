import import_ipynb
from Calc_Entries_CAIR import fasta_to_entry_cair
from Calc_Proteomes_CAIR import *

fasta_to_entry_cair("uniprot_sprot.fasta", "Entries sprot.csv")  # requires the sprot.fasta input file and the output file name
fasta_to_entry_cair("uniprot_trembl.fasta", "Entries trembl.csv")  # requires the trembl.fasta input file and the output file name

entry_to_species()  # insert entries input file and the output file name

species_cair()  # insert species file and the output file name