## 0. Installing necessary packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('ggpubr')){install.packages('ggpubr')}
if(!require('ggplot2')){install.packages('ggplot2')}
if(!require('DTK')){install.packages('DTK')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(data.table)
library(ggpubr)
library(ggplot2)
library(DTK)
library(svglite)

###########################
## 1. Preparing the data ##
#############################################################################################
# rev 26: Table S5 of the article (Table_S5.py in ShannoProt/runNEE); column 14 is the Miller-Madow corrected CEI (CEI_MM),
# which replaces CEI (column 6) in panel D; the file of the working folder is used if available, otherwise the copy in the repository
s5_file <- "Table_S5.csv"
if (!file.exists(s5_file)){
  s5_file <- "../Supplementary_materials/EMIP_supplementary_files/Table_S5.csv"
}
df <- fread(s5_file, select = c(1,2,5,14,7,8,10,11,12))
colnames(df) <- c("stat", "ent", "len", "cair", "int", "mip", "occ", "dis", "age")
df <- df[df$stat=='reviewed',]
df$dis[df$dis != ""] <- 'Disease'
df$dis[df$dis == ""] <- 'Non Disease'
df$age <- factor(df$age, labels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"),
                        levels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"), ordered = T)
occ_sbt <- function (x) {
  if (is.na(x)) {
    x <- NA
  }
  else if (x == 0) {
    x <- '3 - No recorded rare disease'
  }
  else if (x <= 0.000001) {
    x <- '2 - Extremely rare'
  }
  else if (x > 0.000001) {
    x <- '1 - Rare'
  }
}
df$occ <- sapply(df$occ, occ_sbt)
df$stat <- NULL
df$ent <- NULL
df <- df[!is.na(df$occ),]
df <- df[order(occ),]
df <- df[!is.na(df$age)]

#######################################
## 2. Drawing Errorplots w/ Gene Age ##
#############################################################################################
err_size <- 0.35
ball_size <- 0.1
fs <- 6
fsl <- 6
font_family <- "sans"
color_palette <- c("#eaa899", "#b99c74", "#a19997", "#52c0a3",
                   "#0b7791", "#6c5877", "#9b4822", "#cf202d")
mip <- ggerrorplot(data = df, x = 'occ', y = 'mip', desc_stat = 'mean_ci',color = "age",
                   error.plot = "errorbar", xlab = "Occurrence category", size =  err_size,
                   ylab = "NEE", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
int <- ggerrorplot(data = df, x = 'occ', y = 'int', desc_stat = 'mean_ci', color = "age",
                   error.plot = "errorbar", xlab = "Occurrence category", size = err_size,
                   ylab = "Number of interactions", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
len <- ggerrorplot(data = df, x = 'occ', y = 'len', desc_stat = 'mean_ci', color = "age",
                   error.plot = "errorbar", xlab = "Occurrence category", size = err_size,
                   ylab = "Protein length", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
cair <- ggerrorplot(data = df, x = 'occ', y = 'cair', desc_stat = 'mean_ci', color = "age",
                    error.plot = "errorbar", xlab = "Occurrence category", size = err_size,
                    ylab = "CEI_MM", palette = color_palette, add = "mean",
                    position = position_dodge(0.5), add.params = list(size = ball_size)) +
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
fig <- ggarrange(mip, int, len, cair, ncol = 2, nrow = 2, vjust = 0,
                 common.legend = T ,legend = "bottom", align = 'v',
                 font.label = list(size = fsl, face = "bold", family = font_family))

p_mip <- ggplot_build(mip)
p_int <- ggplot_build(int)
p_len <- ggplot_build(len)
p_cair<- ggplot_build(cair)

#############################################
## 3. Performing Dunnett-Tukey-Kramer test ##
#############################################################################################
f <- gl.unequal(n=3, k = c(nrow(df[df$occ == "1 - Rare"]),
                           nrow(df[df$occ == "2 - Extremely rare"]),
                           nrow(df[df$occ == "3 - No recorded rare disease"])))
p_iter <- list(0.05, 0.001, 0.0001)
for (p in p_iter){
  mip_dtk <- as.data.frame(DTK.test(x = df$mip, f = f, a = p)[[2]])
  if (mip_dtk[1,2]*mip_dtk[1,3]>0){
    if (p==0.05){
      mip_sig_2_1 <- "*"
    }
    else if (p==0.001){
      mip_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      mip_sig_2_1 <- "***"
    }
  else {
    mip_sig_2_1 <- "ns"
  }
  }
  if (mip_dtk[2,2]*mip_dtk[2,3]>0){
    if (p==0.05){
      mip_sig_3_1 <- "*"
    }
    else if (p==0.001){
      mip_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      mip_sig_3_1 <- "***"
    }
    else {
      mip_sig_3_1 <- "ns"
    }
  }
  if (mip_dtk[3,2]*mip_dtk[3,3]>0){
    if (p==0.05){
      mip_sig_3_2 <- "*"
    }
    else if (p==0.001){
      mip_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      mip_sig_3_2 <- "***"
    }
    else {
      mip_sig_3_2 <- "ns"
    }
  }
}
mip_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No recorded rare disease", "3 - No recorded rare disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(mip_sig_2_1, mip_sig_3_1, mip_sig_3_2),
                      "y_pos" = c(max(p_mip[["data"]][[2]][["ymax"]])*0.985,
                                  max(p_mip[["data"]][[2]][["ymax"]])*1.065, # rev 26: top bracket raised clear of the stars below
                                  max(p_mip[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  int_dtk <- as.data.frame(DTK.test(x = df$int, f = f, a = p)[[2]])
  if (int_dtk[1,2]*int_dtk[1,3]>0){
    if (p==0.05){
      int_sig_2_1 <- "*"
    }
    else if (p==0.001){
      int_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      int_sig_2_1 <- "***"
    }
    else {
      int_sig_2_1 <- "ns"
    }
  }
  if (int_dtk[2,2]*int_dtk[2,3]>0){
    if (p==0.05){
      int_sig_3_1 <- "*"
    }
    else if (p==0.001){
      int_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      int_sig_3_1 <- "***"
    }
    else {
      int_sig_3_1 <- "ns"
    }
  }
  if (int_dtk[3,2]*int_dtk[3,3]>0){
    if (p==0.05){
      int_sig_3_2 <- "*"
    }
    else if (p==0.001){
      int_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      int_sig_3_2 <- "***"
    }
    else {
      int_sig_3_2 <- "ns"
    }
  }
}
int_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No recorded rare disease", "3 - No recorded rare disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(int_sig_2_1, int_sig_3_1, int_sig_3_2),
                      "y_pos" = c(max(p_int[["data"]][[2]][["ymax"]])*0.99,
                                  max(p_int[["data"]][[2]][["ymax"]])*1.08, # rev 26: top bracket raised clear of the stars below
                                  max(p_int[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  len_dtk <- as.data.frame(DTK.test(x = df$len, f = f, a = p)[[2]])
  if (len_dtk[1,2]*len_dtk[1,3]>0){
    if (p==0.05){
      len_sig_2_1 <- "*"
    }
    else if (p==0.001){
      len_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      len_sig_2_1 <- "***"
    }
    else {
      len_sig_2_1 <- "ns"
    }
  }
  if (len_dtk[2,2]*len_dtk[2,3]>0){
    if (p==0.05){
      len_sig_3_1 <- "*"
    }
    else if (p==0.001){
      len_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      len_sig_3_1 <- "***"
    }
    else {
      len_sig_3_1 <- "ns"
    }
  }
  if (len_dtk[3,2]*len_dtk[3,3]>0){
    if (p==0.05){
      len_sig_3_2 <- "*"
    }
    else if (p==0.001){
      len_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      len_sig_3_2 <- "***"
    }
    else {
      len_sig_3_2 <- "ns"
    }
  }
}
len_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No recorded rare disease", "3 - No recorded rare disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(len_sig_2_1, len_sig_3_1, len_sig_3_2),
                      "y_pos" = c(max(p_len[["data"]][[2]][["ymax"]])*0.99,
                                  max(p_len[["data"]][[2]][["ymax"]])*1.065, # rev 26: top bracket raised clear of the stars below
                                  max(p_len[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  cair_dtk <- as.data.frame(DTK.test(x = df$cair, f = f, a = p)[[2]])
  if (cair_dtk[1,2]*cair_dtk[1,3]>0){
    if (p==0.05){
      cair_sig_2_1 <- "*"
    }
    else if (p==0.001){
      cair_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      cair_sig_2_1 <- "***"
    }
    else {
      cair_sig_2_1 <- "ns"
    }
  }
  if (cair_dtk[2,2]*cair_dtk[2,3]>0){
    if (p==0.05){
      cair_sig_3_1 <- "*"
    }
    else if (p==0.001){
      cair_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      cair_sig_3_1 <- "***"
    }
    else {
      cair_sig_3_1 <- "ns"
    }
  }
  if (cair_dtk[3,2]*cair_dtk[3,3]>0){
    if (p==0.05){
      cair_sig_3_2 <- "*"
    }
    else if (p==0.001){
      cair_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      cair_sig_3_2 <- "***"
    }
    else {
      cair_sig_3_2 <- "ns"
    }
  }
}
cair_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No recorded rare disease", "3 - No recorded rare disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(cair_sig_2_1, cair_sig_3_1, cair_sig_3_2),
                      "y_pos" = c(max(p_cair[["data"]][[2]][["ymax"]])*0.998,
                                  max(p_cair[["data"]][[2]][["ymax"]])*1.0065, # rev 26: top bracket raised clear of the stars below
                                  max(p_cair[["data"]][[2]][["ymax"]])*0.996))

## 3.1. Axis limits (rev 26: each vertical axis starts and ends at a labeled tick)
##########################################################################################
y_limits <- function (p, sig, step) {
  lo <- min(p[["data"]][[2]][["ymin"]]) # lowest end of the error bars
  hi <- max(sig$y_pos)
  hi <- hi + 0.06 * (hi - lo) # room for the significance stars above the highest bracket
  round(c(step * floor(lo / step), step * ceiling(hi / step)), 3)
}
y_scale <- function (lim, step) {
  list(scale_y_continuous(breaks = round(seq(lim[1], lim[2], by = step), 3), expand = c(0, 0)), # no padding on y only
       coord_cartesian(ylim = lim)) # zoom only: no protein is dropped before the means are computed
}
lim_mip <- y_limits(p_mip, mip_sig, 500)
lim_int <- y_limits(p_int, int_sig, 25)
lim_len <- y_limits(p_len, len_sig, 400)
lim_cair <- y_limits(p_cair, cair_sig, 0.025)
mip <- mip + y_scale(lim_mip, 500)
int <- int + y_scale(lim_int, 25)
len <- len + y_scale(lim_len, 400)
cair <- cair + y_scale(lim_cair, 0.025)
fig <- ggarrange(mip, int, len, cair, ncol = 2, nrow = 2, vjust = 0,
                 common.legend = T ,legend = "bottom", align = 'v',
                 font.label = list(size = fsl, face = "bold", family = font_family))

################################
## 4. Drawing overall gglines ##
#############################################################################################
l_mip <- ggline(data = df, x = 'occ', y = 'mip', add = 'mean_ci', error.plot = "errorbar",
              # rev 26: same axis as the colored panel (y_scale below)
              ylab = "NEE", xlab = "Occurrence category", plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = mip_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001) +
  y_scale(lim_mip, 500)

l_int <- ggline(data = df, x = 'occ', y = 'int', add = 'mean_ci', error.plot = "errorbar",
              # rev 26: same axis as the colored panel (y_scale below)
              ylab = "Number of interactions", xlab = "Occurrence category",
              plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = int_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001) +
  y_scale(lim_int, 25)

l_len <- ggline(data = df, x = 'occ', y = 'len', add = 'mean_ci', error.plot = "errorbar",
              # rev 26: same axis as the colored panel (y_scale below)
              ylab = "Protein length", xlab = "Occurrence category",
              plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = len_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001) +
  y_scale(lim_len, 400)

l_cair <- ggline(data = df, x = 'occ', y = 'cair', add = 'mean_ci', error.plot = "errorbar",
               # rev 26: same axis as the colored panel (y_scale below)
               ylab = "CEI_MM", xlab = "Occurrence category", plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = cair_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.004) +
  y_scale(lim_cair, 0.025)

overlay_fig <- ggarrange(l_mip, l_int, l_len, l_cair, labels = c("a", "b", "c", "d"),
                         font.label = list(size = fsl, face = "bold", family = font_family),
                         ncol = 2, nrow = 2)

#############################
## 5. Getting output files ##
#############################################################################################
svglite(file = "Fig5_ColouredAges.svg", width = 7.08, height = 4.75) # rev 26: height of the panels in the published figure (was 4)
fig
dev.off()

svglite(file = "Fig5_BlackMain.svg", width = 7.08, height = 4.75, bg = "transparent") # to be overlaid on "Fig5_ColouredAges.svg"
overlay_fig
dev.off()
