## 0. Installing necessary packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('svglite')){install.packages('svglite')}
if(!require('igraph')){install.packages('igraph')}
if(!require('tidyverse')){install.packages('tidyverse')}
if(!require('Cairo')){install.packages('Cairo')}

## 0.1. Loading packages
########################
library(data.table)
library(svglite)
library(igraph)
library(tidyverse)
library(Cairo)

###########################
## 1. Preparing the data ##
#############################################################################################
df <- fread("../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases_w%uniprot_wo%orpha.csv", # rev 26: the copy in the repository
            select = c(2, 8))
colnames(df) <- c("ent", "mip")
sizes <- df[, c("ent", "mip")]

# loading the uniprot PPI networks file: the human PPIs of UniProt release 2020_02 (with the PICKLE interactions, the larger
# networks would be all black with edges); UniProt reports no interactions for Q8WXI7 (CA125), so its last three PPIs are from PICKLE
links <- fread("../R_codes/EMIP_analyses/Network visualization/UniProt_PPI_List__without_nonhuman_entries.csv") # rev 26: the copy in the repository
colnames(links) <- c("intA", "intB")
links <- links[links$intA != links$intB, ]
links <- unique(links)

############################################
## 2. Customizing the igraph vertex shape ##
#############################################################################################
# DISCLAIMER: Customization of vertex shape is adopted from the main developer
# of iGraph R libarary, Mr. Gábor Csárdi. See: https://github.com/gaborcsardi
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.lty  <- params("vertex", "lty")
  if (length(vertex.frame.lty) != 1 && !is.null(v)) {
    vertex.frame.lty <- vertex.frame.lty[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width, vertex.frame.lty,
         FUN=function(x, y, bg, fg, size, lwd, lty) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd, lty=lty,
                   circles=size, add=TRUE, inches=FALSE)
         })
}
add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.lty=1,
                                                vertex.frame.width=1))

############################################
## 3. Preparing and drwaing the networks  ##
#############################################################################################
protnet <- function(Label, Entry, Protein_name, Output_file, Title_pos) {
prim_edge <- c()
main_nodes <- c(Entry)
for (i in main_nodes){
  prim_edge = rbind(links[links$intA==i | links$intB==i, ], prim_edge)
}
ents <- c(prim_edge$intA, prim_edge$intB)
ents <- unique(ents)
scnd_edge = c()
for (i in ents){
  scnd_edge = rbind(links[links$intA==i | links$intB==i, ], scnd_edge)
}
scnd_edge <- scnd_edge[, c(1,2)]
sources <- scnd_edge %>%
  distinct(intA) %>%
  rename(ent = intA)
destinations <- scnd_edge %>%
  distinct(intB) %>%
  rename(ent = intB)
nodes <- full_join(sources, destinations, by = "ent")
nodes <- left_join(nodes, sizes, by = "ent")
nodes$mip <- replace_na(nodes$mip, 0)                    
net <- graph.data.frame(scnd_edge, nodes, directed=F)
vrtx <- V(net)
svglite(file = Output_file, width = 1.77, height = 1.77)
par(mar=c(1,1,1,1))
plot(x = net,
     layout=layout.kamada.kawai,
     vertex.shape="fcircle",
     vertex.frame.color="#000000",
     vertex.lty=1,
     vertex.size=0.04*(vrtx$mip)**0.6,
     vertex.frame.width=0.25,
     vertex.label=NA,
     edge.curved=0,
     vertex.color=c("#1587a0"),
     edge.color="#262d3b",
     edge.width=0.12)
title(Protein_name, line = -(7+Title_pos/5), font.main = 1, family="sans", cex.main=0.5)
title(Label, adj=0, line = -0.36, font.main = 2, family="sans", cex.main=0.5)
dev.off()
}

#############################
## 4. Getting output files ##
#############################################################################################
protnet('A', 'Q8WZ42', 'Connectin', 'Fig6A.svg', 1)
protnet('B', 'P05067', 'Amyloid-beta\nprecursor protein', 'Fig6B.svg', 2)
protnet('C', 'P0CG48', 'Polyubiquitin-C', 'Fig6C.svg', 1)
protnet('D', 'Q8WXI7', 'Ovarian carcinoma\nantigen CA125', 'Fig6D.svg', 2) # the drawn network is from the PICKLE database
protnet('E', 'Q9NRI5', 'Disrupted in schizo-\nphrenia 1 protein', 'Fig6E.svg', 2)
protnet('F', 'P04637', 'Cellular tumor\nantigen p53', 'Fig6F.svg', 2)
protnet('G', 'Q09472', 'Histone acetyl-\ntransferase p300', 'Fig6G.svg', 2)
protnet('H', 'P00533', 'Epidermal growth\nfactor receptor', 'Fig6H.svg', 2)
protnet('I', 'P62993', 'Growth factor receptor-\nbound protein 2', 'Fig6I.svg', 2)
protnet('J', 'Q8NF91', 'Nesprin-1', 'Fig6J.svg', 1)
protnet('K', 'P63104', 'Protein kinase C\ninhibitor protein 1', 'Fig6K.svg', 2)
protnet('L', 'P78362', 'SRSF protein kinase 2', 'Fig6L.svg', 2)
protnet('M', 'Q03001', 'Dystonin', 'Fig6M.svg', 1)
protnet('N', 'Q5VST9', 'Obscurin', 'Fig6N.svg', 1)
protnet('O', 'P38398', 'Breast cancer type 1\nsusceptibility protein', 'Fig6O.svg', 2)
protnet('P', 'Q15149', 'Plectin', 'Fig6P.svg', 1)
