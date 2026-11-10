## 0. Installing necessary packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('HH')){install.packages('HH')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(data.table)
library(HH)
library(svglite)

###########################
## 1. Preparing the data ##
#############################################################################################
df <- fread("../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv", # rev 26: the copy in the repository
            select = c(1,2,5,6,7,8,10,11,12))
colnames(df) <- c("stat", "ent", "len", "cair", "int", "mip", "occ", "dis", "age")
df <- df[df$stat=='reviewed',]
df$dis[df$dis != ""] <- 'Disease'
df$dis[df$dis == ""] <- 'Non Disease'
df$mip <- floor(10*log(df$mip))
df$stat <- NULL
df$ent <- NULL

## 1.1. Ranking Total Occurences in three categories
####################################################
tot_occ_sbt <- function (x) {
  if (is.na(x)) {
    x <- NA
  }
  else if (x == 0) {
    x <- '3 - No Disease'
  }
  else if (x <= 0.000001) {
    x <- '2 - Extremely rare'
  }
  else if (x > 0.000001) {
    x <- '1 - Rare'
  }
}
df$occ <- sapply(df$occ, tot_occ_sbt)

## 1.2. Renaming 'Gene Ages' to Ranks
#####################################
df$age <- factor(df$age, labels = c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7", "Rank 8"),
                 levels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"), ordered = T)

## 1.3. Ranking 'MIP' in eight categories
#########################################
m <- mean(df$mip)
sd <- sd(df$mip)
if (m-3*sd > min(df$mip) & max(df$mip) > m+3*sd) {
  prcnts <- c(m-3*sd, m-2*sd, m-sd, m, m+sd, m+2*sd, m+3*sd, max(df$mip))
} else {
    prcnts <- as.vector(quantile(df$mip,
                      c(0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1)))
}
mut_inf_sbt <- function (x) {
  if (x <= prcnts[1]){
    x <- 'Rank 1'
  }
  else if (x <= prcnts[2]){
    x <- 'Rank 2'
  }
  else if (x <= prcnts[3]){
    x <- 'Rank 3'
  }
  else if (x <= prcnts[4]){
    x <- 'Rank 4'
  }
  else if (x <= prcnts[5]){
    x <- 'Rank 5'
  }
  else if (x <= prcnts[6]){
    x <- 'Rank 6'
  }
  else if (x <= prcnts[7]){
    x <- 'Rank 7'
  }
  else {x <- 'Rank 8'}
}
df$mip <- sapply(df$mip, mut_inf_sbt)

## 1.4. Ranking 'Number of Interactions' in eight categories
############################################################
m <- mean(df$int)
sd <- sd(df$int)
if (m-3*sd > min(df$int) & max(df$int) > m+3*sd) {
  prcnts <- c(m-3*sd, m-2*sd, m-sd, m, m+sd, m+2*sd, m+3*sd, max(df$int))
} else {
  prcnts <- as.vector(quantile(df$int,
                               c(0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1)))
}
num_int_sbt <- function (x) {
  if (x <= prcnts[1]){
    x <- 'Rank 1'
  }
  else if (x <= prcnts[2]){
    x <- 'Rank 2'
  }
  else if (x <= prcnts[3]){
    x <- 'Rank 3'
  }
  else if (x <= prcnts[4]){
    x <- 'Rank 4'
  }
  else if (x <= prcnts[5]){
    x <- 'Rank 5'
  }
  else if (x <= prcnts[6]){
    x <- 'Rank 6'
  }
  else if (x <= prcnts[7]){
    x <- 'Rank 7'
  }
  else {x <- 'Rank 8'}
}
df$int <- sapply(df$int, num_int_sbt)

## 1.5. Ranking 'CAIR' in eight categories
##########################################
m <- mean(df$cair)
sd <- sd(df$cair)
if (m-3*sd > min(df$cair) & max(df$cair) > m+3*sd) {
  prcnts <- c(m-3*sd, m-2*sd, m-sd, m, m+sd, m+2*sd, m+3*sd, max(df$cair))
} else {
  prcnts <- as.vector(quantile(df$cair,
                               c(0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1)))
}
prot_cair_sbt <- function (x) {
  if (x <= prcnts[1]){
    x <- 'Rank 1'
  }
  else if (x <= prcnts[2]){
    x <- 'Rank 2'
  }
  else if (x <= prcnts[3]){
    x <- 'Rank 3'
  }
  else if (x <= prcnts[4]){
    x <- 'Rank 4'
  }
  else if (x <= prcnts[5]){
    x <- 'Rank 5'
  }
  else if (x <= prcnts[6]){
    x <- 'Rank 6'
  }
  else if (x <= prcnts[7]){
    x <- 'Rank 7'
  }
  else {x <- 'Rank 8'}
}
df$cair <- sapply(df$cair, prot_cair_sbt)

## 1.6. Ranking 'Protein Length' in eight categories
####################################################
m <- mean(df$len)
sd <- sd(df$len)
if (m-3*sd > min(df$len) & max(df$len) > m+3*sd) {
  prcnts <- c(m-3*sd, m-2*sd, m-sd, m, m+sd, m+2*sd, m+3*sd, max(df$len))
} else {
  prcnts <- as.vector(quantile(df$len,
                               c(0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1)))
}
prot_len_sbt <- function (x) {
  if (x <= prcnts[1]){
    x <- 'Rank 1'
  }
  else if (x <= prcnts[2]){
    x <- 'Rank 2'
  }
  else if (x <= prcnts[3]){
    x <- 'Rank 3'
  }
  else if (x <= prcnts[4]){
    x <- 'Rank 4'
  }
  else if (x <= prcnts[5]){
    x <- 'Rank 5'
  }
  else if (x <= prcnts[6]){
    x <- 'Rank 6'
  }
  else if (x <= prcnts[7]){
    x <- 'Rank 7'
  }
  else {x <- 'Rank 8'}
}
df$len <- sapply(df$len, prot_len_sbt)

#############################
## 2. Drawing likert plots ##
#############################################################################################
col <- c("#cd9289", "#d33a41", "#0b7791")
fontlab <- list(
  font = 2,
  cex = 0.5,
  fontfamily = "sans")
fonttext <- list(
  font = 1,
  cex = 0.5,
  fontfamily = "sans")
parset <- list(
  par.xlab.text = fontlab,
  par.ylab.text = fontlab,
  axis.text = fonttext,
  sub.text = fonttext,
  add.text = fonttext)

svglite(file = 'Fig5A.svg', width = 3.4, height = 2.5)
xt <- t(as.data.frame.array(xtabs(~ occ + mip, data = df)))
plot.likert(x = xt, ylab = 'LEMIP', xlab = "Percents", as.percent = T,
            ReferenceZero = 2.5, scales=list(x=list(at=seq(-100,100,10))),
            auto.key = list(between=0.3, between.columns=1, cex = 0.5),
            main = F, rightAxis = F, col = col, par.settings = parset)
dev.off()

svglite(file = 'Fig5B.svg', width = 3.4, height = 2.5)
xt <- t(as.data.frame.array(xtabs(~ occ + int, data = df)))
plot.likert(x = xt, ylab = 'Number of interactions', xlab = "Percents", as.percent = T,
            ReferenceZero = 2.5, scales=list(x=list(at=seq(-100,100,10))),
            auto.key = list(between=0.3, between.columns=1, cex = 0.5),
            main = F, rightAxis = F, col = col, par.settings = parset)
dev.off()

svglite(file = 'Fig5D.svg', width = 3.4, height = 2.5)
xt <- t(as.data.frame.array(xtabs(~ occ + cair, data = df)))
plot.likert(x = xt, ylab = 'CAIR', xlab = "Percents", as.percent = T,
            ReferenceZero = 2.5, scales=list(x=list(at=seq(-100,100,10))),
            auto.key = list(between = 0.3, between.columns = 1, cex = 0.5),
            main = F, rightAxis = F, col = col, par.settings = parset)
dev.off()

svglite(file = 'Fig5E.svg', width = 3.4, height = 2.5)
xt <- t(as.data.frame.array(xtabs(~ occ + len, data = df)))
plot.likert(x = xt, ylab = 'Protein length', xlab = "Percents", as.percent = T,
            ReferenceZero = 2.5, scales=list(x=list(at=seq(-100,100,10))),
            auto.key = list(between = 0.3, between.columns = 1, cex = 0.5),
            main = F, rightAxis = F, col = col, par.settings = parset)
dev.off()

svglite(file = 'Fig5C.svg', width = 3.4, height = 2.5)
df <- df[!is.na(df$age),]
xt <- t(as.data.frame.array(xtabs(~ occ + age, data = df)))
plot.likert(x = xt, ylab = 'Gene age', xlab = "Percents", as.percent = T,
            ReferenceZero = 2.5, scales=list(x=list(at=seq(-100,100,10))),
            auto.key=list(between = 0.3, between.columns = 1, cex = 0.5),
            main = F, rightAxis = F, col = col, par.settings = parset)
dev.off()
