## Bifurcation tests (rev 26) ##
#############################################################################################
# Brunner-Munzel tests of CEI (formerly CAIR) for the bifurcations of the tree of life, i.e. the
# same comparisons as "Violin Plots.R", tabulated at three levels:
#   proteome level: one value per proteome (Table S1)
#   species level:  median CEI per species (pseudoreplication check)
#   genus level:    median CEI per genus
# Added in rev 26 (as described in Methods section):
#   * range-preserving (logit) 95% confidence intervals of the relative effect
#   * Benjamini-Hochberg adjustment over the testable bifurcations
#   * Vargha-Delaney magnitude classes of the relative effect
#   * studentized permutation tests (Neubert & Brunner 2007) when a clade has fewer than 30 proteomes
# The relative effect (RE) is P(X<Y)+0.5*P(X=Y), X = upper clade (first group), Y = lower clade,
# i.e. the probability that an organism of the upper clade has the lower CEI.
# Input:  Complete_proteome_CAIRs.csv (Table S1); optionally "GC null per proteome.csv" (section 4)
# Output: Bifurcation tests (proteome).csv, (species).csv, (genus).csv, Bifurcation tests summary.csv
#         and, with the GC file, the same tables for the GC-adjusted CEI (Table S3)

## 0. Installing necessary packages ##
#############################################################################################
if(!require('brunnermunzel')){install.packages('brunnermunzel')}
if(!require('data.table')){install.packages('data.table')}

## 0.1. Loading packages
########################
library(brunnermunzel) # for "brunnermunzel.test"
library(data.table) # for "fread"

###########################
## 1. Preparing the data ##
#############################################################################################
# Table S1 in the repository (the script is run from a folder of the repository, e.g. Results)
cair_file <- "../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv"
proteomes <- fread(cair_file, select = c(1,2,3,5))
colnames(proteomes) <- c("org", "cair", "tax", "code")
n_perm <- 10000 # number of permutations of the studentized permutation test
small_clade <- 30 # clades smaller than this are also tested by permutation
set.seed(2026) # to make the permutation p-values reproducible
df <- proteomes # the data used by "bmtest" (switched between levels in section 3)
permute <- TRUE # whether "bmtest" runs the permutation test

## 1.1. One value per species and per genus
###########################################
# species = first two words of the organism name, genus = first word, after removing
# "Candidatus", "uncultured", "unclassified", brackets and quotes; value = median CEI;
# phylum code = the most frequent code of the organisms of the species (or genus)
name_words <- function (x) {
  x <- gsub("\\[|\\]|'|\"", "", x)
  words <- strsplit(trimws(x), "\\s+")[[1]]
  while (length(words) > 0 && tolower(words[1]) %in% c("candidatus", "uncultured", "unclassified")) {
    words <- words[-1]
  }
  words
}
words <- lapply(proteomes$tax, name_words)
proteomes$genus <- sapply(words, function (w) if (length(w) > 0) w[1] else "NA")
proteomes$species <- sapply(words, function (w) if (length(w) > 1) paste(w[1], w[2]) else if (length(w) == 1) w[1] else "NA")
mode_code <- function (x) as.numeric(names(which.max(table(x)))) # the smallest code wins ties
species <- proteomes[, .(cair = median(cair), code = mode_code(code)), by = species]
genera <- proteomes[, .(cair = median(cair), code = mode_code(code)), by = genus]

#######################################################
## 2. Creating a function for the Brunner-Munzel test ##
#############################################################################################
vd_class <- function (re) {
  # Vargha & Delaney (2000) thresholds on |RE - 0.5|
  d <- abs(re - 0.5)
  if (is.na(d)) {
    NA
  } else if (d < 0.06) {
    'negligible'
  } else if (d < 0.14) {
    'small'
  } else if (d < 0.21) {
    'medium'
  } else {
    'large'
  }
}

bm_permutation <- function (x, y) {
  # studentized permutation test: the pooled values are permuted and the test statistic recalculated
  t_obs <- brunnermunzel.test(x = x, y = y)$statistic
  z <- c(x, y)
  t_perm <- replicate(n_perm, {
    zz <- sample(z)
    brunnermunzel.test(x = zz[seq_along(x)], y = zz[-seq_along(x)])$statistic
  })
  (sum(abs(t_perm) >= abs(t_obs) - 1e-12) + 1)/(n_perm + 1)
}

bmtest <- function(test_num, grp1, grp2){
  x <- df$cair[grp1] # upper clade
  y <- df$cair[grp2] # lower clade
  out <- data.frame(test = test_num, n1 = length(x), n2 = length(y), re = NA, ci_low = NA, ci_high = NA,
                    logit_low = NA, logit_high = NA, df = NA, statistic = NA, p = NA, p_perm = NA)
  if (min(length(x), length(y)) >= 3){
    BMtest <- brunnermunzel.test(x = x, y = y)
    re <- as.numeric(BMtest$estimate)
    tq <- qt(0.975, BMtest$parameter) # the t quantile used by brunnermunzel.test
    se <- (BMtest$conf.int[2] - re)/tq # standard error of the relative effect
    # logit-scale interval, back-transformed, so that the limits stay inside [0, 1]
    logit_se <- se/(re*(1 - re))
    out[, c("re", "ci_low", "ci_high", "df", "statistic", "p")] <- c(re, BMtest$conf.int, BMtest$parameter,
                                                                    BMtest$statistic, BMtest$p.value)
    out$logit_low <- plogis(qlogis(re) - tq*logit_se)
    out$logit_high <- plogis(qlogis(re) + tq*logit_se)
    if (permute & min(length(x), length(y)) < small_clade){
      out$p_perm <- bm_permutation(x, y)
    }
  }
  out
}

###############################################
## 3. Running the tests of all bifurcations ##
#############################################################################################
# WARNING: the bifurcations are those of "Violin Plots.R"; the tests with ten or fewer
# proteomes in a clade are hashed there and are kept here as well, so that they can be
# counted at species and genus level. Three definitions are corrected (T54, T73, T75; see
# the comments), and T64 (8 vs 67 proteomes) is tabulated but not tested at proteome level.
run_tests <- function(){
  rbind(
    bmtest('T1', df$code<=56, df$code>=57),
    bmtest('T2', df$code<=55, df$code==56),
    bmtest('T3', df$code<=51, df$code<=55 & df$code>=52),
    bmtest('T4', df$code<=42, df$code<=51 & df$code>=43),
    bmtest('T5', df$code<=34, df$code<=42 & df$code>=35),
    bmtest('T6', df$code<=18, df$code<=34 & df$code>=19),
    bmtest('T7', df$code<=10, df$code<=18 & df$code>=11),
    bmtest('T8', df$code<=6, df$code<=10 & df$code>=7),
    bmtest('T9', df$code<=4, df$code<=6 & df$code>=5),
    bmtest('T10', df$code==1, df$code<=4 & df$code>=2),
  # bmtest('T11', df$code<=3 & df$code>=2, df$code==4),
  # bmtest('T12', df$code==2, df$code==3),
  # bmtest('T13', df$code==5, df$code==6),
    bmtest('T14', df$code==7, df$code<=10 & df$code>=8),
    bmtest('T15', df$code==8, df$code<=10 & df$code>=9),
    bmtest('T16', df$code==9, df$code==10),
    bmtest('T17', df$code==11, df$code<=18 & df$code>=12),
    bmtest('T18', df$code==12, df$code<=18 & df$code>=13),
    bmtest('T19', df$code==13, df$code<=18 & df$code>=14),
    bmtest('T20', df$code==14, df$code<=18 & df$code>=15),
    bmtest('T21', df$code<=16 & df$code>=15, df$code<=18 & df$code>=17),
  # bmtest('T22', df$code==15, df$code==16),
  # bmtest('T23', df$code==17, df$code==18),
  # bmtest('T24', df$code==19, df$code<=34 & df$code>=20),
  # bmtest('T25', df$code==20, df$code<=34 & df$code>=21),
    bmtest('T26', df$code<=22 & df$code>=21, df$code<=34 & df$code>=23),
  # bmtest('T27', df$code==21, df$code==22),
    bmtest('T28', df$code<=32 & df$code>=23, df$code<=34 & df$code>=33),
    bmtest('T29', df$code==23, df$code<=32 & df$code>=24),
  # bmtest('T30', df$code==24, df$code<=32 & df$code>=25),
  # bmtest('T31', df$code==25, df$code<=32 & df$code>=26),
    bmtest('T32', df$code<=29 & df$code>=26, df$code<=32 & df$code>=30),
    bmtest('T33', df$code<=27 & df$code>=26, df$code<=29 & df$code>=28),
  # bmtest('T34', df$code==26, df$code==27),
  # bmtest('T35', df$code==28, df$code==29),
    bmtest('T36', df$code==30, df$code<=32 & df$code>=31),
  # bmtest('T37', df$code==31, df$code==32),
    bmtest('T38', df$code==33, df$code==34),
    bmtest('T39', df$code<=36 & df$code>=35, df$code<=42 & df$code>=37),
    bmtest('T40', df$code==35, df$code==36),
    bmtest('T41', df$code==37, df$code<=42 & df$code>=38),
    bmtest('T42', df$code<=39 & df$code>=38, df$code<=42 & df$code>=40),
  # bmtest('T43', df$code==38, df$code==39),
    bmtest('T44', df$code==40, df$code<=42 & df$code>=41),
  # bmtest('T45', df$code==41, df$code==42),
    bmtest('T46', df$code<=49 & df$code>=43, df$code<=51 & df$code>=50),
    bmtest('T47', df$code<=48 & df$code>=43, df$code==49),
  # bmtest('T48', df$code<=47 & df$code>=43, df$code==48),
    bmtest('T49', df$code<=45 & df$code>=43, df$code<=47 & df$code>=46),
  # bmtest('T50', df$code<=44 & df$code>=43, df$code==45),
    bmtest('T51', df$code==43, df$code==44),
  # bmtest('T52', df$code==46, df$code==47),
  # bmtest('T53', df$code==50, df$code==51),
    bmtest('T54', df$code<=53 & df$code>=52, df$code<=55 & df$code>=54), # corrected in rev 26: "df$code>=42" compared codes 42-53 instead of codes 52-53
    bmtest('T55', df$code==52, df$code==53),
    bmtest('T56', df$code==54, df$code==55),
    bmtest('T57', df$code<=66 & df$code>=57, df$code>=67),
  # bmtest('T58', df$code==57, df$code<=66 & df$code>=58),
    bmtest('T59', df$code<=60 & df$code>=58, df$code<=66 & df$code>=61),
    bmtest('T60', df$code==58, df$code<=60 & df$code>=59),
  # bmtest('T61', df$code==59, df$code==60),
    bmtest('T62', df$code==61, df$code<=66 & df$code>=62),
  # bmtest('T63', df$code==62, df$code<=66 & df$code>=63),
    bmtest('T64', df$code<=64 & df$code>=63, df$code<=66 & df$code>=65), # 8 vs 67 proteomes: tabulated but not tested
  # bmtest('T65', df$code==63, df$code==64),
    bmtest('T66', df$code==65, df$code==66),
    bmtest('T67', df$code==67, df$code>=68),
    bmtest('T68', df$code<=75 & df$code>=68, df$code<=92 & df$code>=76),
    bmtest('T69', df$code<=70 & df$code>=68, df$code<=75 & df$code>=71),
    bmtest('T70', df$code==68, df$code<=70 & df$code>=69),
  # bmtest('T71', df$code==69, df$code==70),
  # bmtest('T72', df$code==71, df$code<=75 & df$code>=72),
    bmtest('T73', df$code<=73 & df$code>=72, df$code<=75 & df$code>=74), # corrected in rev 26: "df$code>=62" instead of "df$code>=72"; 4 vs 150 proteomes: not tested
  # bmtest('T74', df$code==72, df$code==73),
    bmtest('T75', df$code==74, df$code==75), # corrected in rev 26: "df$code>=75" compared with all lineages above code 74 instead of code 75
    bmtest('T76', df$code==76, df$code>=77),
    bmtest('T77', df$code==77, df$code>=78),
  # bmtest('T78', df$code==78, df$code>=79),
    bmtest('T79', df$code<=87 & df$code>=79, df$code<=92 & df$code>=88),
    bmtest('T80', df$code==79, df$code<=87 & df$code>=80),
    bmtest('T81', df$code<=86 & df$code>=80, df$code==87),
  # bmtest('T82', df$code==80, df$code<=86 & df$code>=81),
    bmtest('T83', df$code<=82 & df$code>=81, df$code<=86 & df$code>=83),
  # bmtest('T84', df$code==81, df$code==82),
  # bmtest('T85', df$code==83, df$code<=86 & df$code>=84),
  # bmtest('T86', df$code==84, df$code<=86 & df$code>=85),
  # bmtest('T87', df$code==85, df$code==86),
    bmtest('T88', df$code<=90 & df$code>=88, df$code<=92 & df$code>=91),
  # bmtest('T89', df$code<=89 & df$code>=88, df$code==90),
  # bmtest('T90', df$code==88, df$code==89),
  # bmtest('T91', df$code==91, df$code==92),
    NULL)
}

## 3.1. Testing at proteome, species and genus level
####################################################
test_levels <- function(value_name){
  for (level in c("proteome", "species", "genus")){
    df <<- switch(level, proteome = proteomes, species = species, genus = genera)
    permute <<- (level == "proteome") # permutation tests only at proteome level
    tests <- run_tests()
    if (level == "proteome"){
      tests$testable <- tests$n1 > 10 & tests$n2 > 10 # ten or fewer proteomes in a clade: not tested
    } else {
      tests$testable <- tests$n1 >= 10 & tests$n2 >= 10 # at least ten species (or genera) on each side
    }
    tests$p_adj <- NA
    tests$p_adj[tests$testable] <- p.adjust(tests$p[tests$testable], method = "BH")
    tests$class <- sapply(tests$re, vd_class)
    tests$separated <- tests$testable & tests$p_adj < 0.05 & tests$class != 'negligible'
    tests$level <- level
    fwrite(tests, paste0("Bifurcation tests", value_name, " (", level, ").csv"))
    assign(paste0("tests_", level), tests, envir = .GlobalEnv)
  }
}

## 3.2. Summary of the tests
############################
summarize_tests <- function(value_name){
  tp <- tests_proteome[tests_proteome$testable, ]
  sig <- tp[tp$p_adj < 0.05, ]
  small <- tp[!is.na(tp$p_perm), ]
  test_list <- function(x) paste(x, collapse = " ")
  smry <- data.frame(item = c("proteomes", "testable", "significant", "large", "medium", "small", "negligible",
                              "separated", "ge_medium", "nonsignificant", "negligible_significant", "not_tested",
                              "small_clade_tests", "permutation_same_conclusion",
                              "wald_outside_01"),
                     value = c(nrow(proteomes), nrow(tp), nrow(sig), sum(sig$class == 'large'), sum(sig$class == 'medium'),
                               sum(sig$class == 'small'), sum(sig$class == 'negligible'), sum(tp$separated),
                               sum(sig$class %in% c('medium', 'large')),
                               test_list(paste0(tp$test[tp$p_adj >= 0.05], ":", tp$class[tp$p_adj >= 0.05])),
                               test_list(sig$test[sig$class == 'negligible']),
                               test_list(paste0(tests_proteome$test[!tests_proteome$testable], ":", tests_proteome$n1[!tests_proteome$testable],
                                                "vs", tests_proteome$n2[!tests_proteome$testable])),
                               nrow(small), sum((small$p < 0.05) == (small$p_perm < 0.05)),
                               test_list(tp$test[tp$ci_low < 0 | tp$ci_high > 1])))
  for (level in c("species", "genus")){
    tl <- get(paste0("tests_", level))
    both <- merge(tp, tl[tl$testable, ], by = "test", suffixes = c("", "_l")) # testable at both levels
    sep_testable <- both[both$separated, ]
    smry <- rbind(smry, data.frame(item = paste0(level, "_", c("units", "testable", "significant", "rho_effect", "same_direction",
                                                             "separated_testable", "still_separated", "lost",
                                                             "median_abs_effect_proteome", "median_abs_effect_level")),
                                   value = c(nrow(get(ifelse(level == "species", "species", "genera"))), sum(tl$testable),
                                             sum(tl$p_adj < 0.05, na.rm = T),
                                             cor(both$re, both$re_l, method = "spearman"),
                                             paste0(sum(sign(both$re - 0.5) == sign(both$re_l - 0.5)), "/", nrow(both)),
                                             nrow(sep_testable), sum(sep_testable$separated_l),
                                             test_list(sep_testable$test[!sep_testable$separated_l]),
                                             median(abs(both$re - 0.5)), median(abs(both$re_l - 0.5)))))
  }
  fwrite(smry, paste0("Bifurcation tests", value_name, " summary.csv"))
  print(smry)
}

test_levels("")
summarize_tests("")

#####################################################
## 4. Tests on the GC-adjusted CEI (Table S3)       ##
#############################################################################################
# delta CEI = CEI - CEI expected under random translation at the genome GC content
# (from "GC null per proteome.csv", Proteome_GC_null.py); only the proteomes with a GC value
gc_file <- "GC null per proteome.csv"
if (!file.exists(gc_file)){
  gc_file <- "../Supplementary_materials/CAIR_supplementary_files/GC null per proteome.csv" # the copy in the repository
}
if (file.exists(gc_file)){
  gc <- fread(gc_file, select = c("Organism_ID", "Delta CEI"))
  colnames(gc) <- c("org", "dcei")
  gc <- gc[!is.na(gc$dcei), ] # proteomes without a matched genome assembly are left out (not counted as missing values)
  proteomes <- merge(proteomes, gc, by = "org")
  proteomes$cair <- proteomes$dcei # the tests below use delta CEI in place of CEI
  species <- proteomes[, .(cair = median(cair), code = mode_code(code)), by = species]
  genera <- proteomes[, .(cair = median(cair), code = mode_code(code)), by = genus]
  test_levels(" of delta CEI")
  summarize_tests(" of delta CEI")

  ## 4.1. Agreement with the tests of CEI (proteome level)
  cei_tests <- fread("Bifurcation tests (proteome).csv")
  both <- merge(cei_tests[cei_tests$testable, ], tests_proteome[tests_proteome$testable, ], by = "test", suffixes = c("", "_d"))
  both <- both[order(as.numeric(sub("T", "", both$test))), ] # T1, T2, ... (not T1, T10, ...)
  comparison <- data.frame(item = c("tested", "rho_effect_vs_cei", "same_direction_vs_cei", "separated_on_cei", "separated_on_both",
                                    "lost_after_gc_adjustment", "gained_after_gc_adjustment"),
                           value = c(nrow(both), cor(both$re, both$re_d, method = "spearman"),
                                     paste0(sum(sign(both$re - 0.5) == sign(both$re_d - 0.5)), "/", nrow(both)),
                                     sum(both$separated), sum(both$separated & both$separated_d),
                                     paste(both$test[both$separated & !both$separated_d], collapse = " "),
                                     paste(both$test[!both$separated & both$separated_d], collapse = " ")))
  fwrite(comparison, "Bifurcation tests of delta CEI vs CEI.csv")
  print(comparison)
}
