## 0. Installing necessary packages ##
#############################################################################################
if(!require('readxl')){install.packages('readxl')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(readxl) # for "read_excel"
library(svglite) # for "svglite"

################################
## 1. Fetching required files ##
#############################################################################################
temp <- "../R_codes/CAIR_comparisons/Wave of Life/WoL_to_E_coli_removed_organism_duplicates.xlsx" # rev 26: the copy in the repository (run from a folder of the repository, e.g. Results)

##############################
## 2. Drawing density plots ##
#############################################################################################
bw <- 0.0003 # bandwidth of density plots
family <- "sans" # font family of the plots
width <- 3.54 # width of the output plots in inches
height <- 1.82 # height of the output plots in inches
pointsize <- 6 # pointsize of the vector image outputs

## 2.1. Wave of life, i.e. all organisms
########################################
all_organisms <- read_excel(temp, sheet = "All organisms")
from <- 0.875
to <- 0.945
wave_life_dens <- density(all_organisms$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3A.svg', pointsize = pointsize, width = width*2, height = height)
plot(wave_life_dens, main = "", col="#cf202d",
     xlab='CAIR of all organisms', xaxt='n',
     ylab='', yaxt='n', family=family)
axis(side=1, at=seq(from = from, to = to, by = 0.005))
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.2. Proteobacteria
######################
proteobacteria <- read_excel(temp, sheet = "Proteobacteria")
from <- 0.890
to <- 0.936
proteo_dens <- density(proteobacteria$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3B.svg', pointsize = pointsize, width = width, height = height)
plot(proteo_dens, main = "", col="#faa41a",
     xlab='CAIR of the Proteobacteria phylum',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.3. Gammaproteobacteria
###########################
gammaproteobacteria <- read_excel(temp, sheet = "Gammaproteobacteria")
from <- 0.900
to <- 0.936
gammaproteo_dens <- density(gammaproteobacteria$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3C.svg', pointsize = pointsize, width = width, height = height)
plot(gammaproteo_dens, main = "", col="#0b7791",
     xlab='CAIR of the Gammaproteobacteria class',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.4. Enterobacterales
########################
enterobacterales <- read_excel(temp, sheet = "Enterobacterales")
from <- 0.910
to <- 0.936
entero_dens <- density(enterobacterales$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3D.svg', pointsize = pointsize, width = width, height = height)
plot(entero_dens, main = "", col="#4b5437",
     xlab='CAIR of the Enterobacterales order',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.5. Enterobacteriaceae
##########################
enterobacteriaceae <- read_excel(temp, sheet = "Enterobacteriaceae")
from <- 0.926
to <- 0.935
enterobacter_dens <- density(enterobacteriaceae$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3E.svg', pointsize = pointsize, width = width, height = height)
plot(enterobacter_dens, main = "", col="#9b4822",
     xlab='CAIR of the Enterobacteriaceae family',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.6. Escherichia
###################
escherichia <- read_excel(temp, sheet = "Escherichia")
from <- 0.930
to <- 0.935
escherichia_dens <- density(escherichia$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3F.svg', pointsize = pointsize, width = width, height = height)
plot(escherichia_dens, main = "", col="#6c5877",
     xlab='CAIR of the Escherichia genus',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 2.7. E. coli
###############
e_coli <- read_excel(temp, sheet = "E. coli")
from <- 0.931
to <- 0.935
e_coli_dens <- density(e_coli$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3G.svg', pointsize = pointsize, width = width, height = height)
plot(e_coli_dens, main = "", col="#52c0a3",
     xlab='CAIR of the E. coli organism',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()
