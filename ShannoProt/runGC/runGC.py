import os
import sys
from Proteome_GC_null import proteome_gc_null, gc_null_summary, null_maxima, cei_distribution_shape
from CDS_codon_position_GC import codon_position_gc, codon_position_summary
from GTDB_phylogeny_data import gtdb_matching, gtdb_species_values, gc_ceiling_quality

cair_file = "Complete_proteome_CAIRs.csv"  # Table S1 of the preprint: a copy in the working folder, otherwise the file of the repository
if not os.path.exists(cair_file):
    cair_file = "../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv"
proteome_gc_null(cair_file, None, "nodes.dmp", "merged.dmp",
                 "GC null per proteome.csv")  # requires Table S1, the four NCBI assembly summaries and the NCBI taxonomy (nodes.dmp, merged.dmp)
gc_null_summary("GC null per proteome.csv")  # R2 of CEI on GC, shares above the null by domain and phylum, the null curve of Figure 3
null_maxima(11)  # the highest CEI of the uniform-GC, codon-position GC and free-composition nulls
cei_distribution_shape("GC null per proteome.csv")  # skew of the CEI distribution and of the null

if "--cds" in sys.argv:  # sensitivity to codon-position GC: downloads the coding sequences of one genome per genus (about 11 GB, 1 h)
    codon_position_gc("GC null per proteome.csv", "cds", "Codon position GC per genus.csv")
    codon_position_summary("Codon position GC per genus.csv")

if os.path.exists("bac120_metadata_r232.tsv.gz"):  # rev 26: GTDB genomes of the proteomes (requires the GTDB release 232 metadata files)
    gtdb_matching("GC null per proteome.csv", ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz"), "Proteomes GTDB.csv")
    gtdb_species_values("Proteomes GTDB.csv", "GTDB species values.csv")  # input of "Phylogenetic signal.R"
    gc_ceiling_quality("Proteomes GTDB.csv", "GC null maxima.csv", "Codon position GC per genus.csv",
                       "GC ceiling by genome quality.csv")  # proteomes above the uniform-GC maximum by genome quality
