## Phylogenetic signal (rev 26) ##
#############################################################################################
# Phylogenetic signal (Pagel's lambda) of CEI (formerly CAIR), genome GC, the CEI expected under
# random translation at the genome GC content (CEI null) and delta CEI, and phylogenetic
# generalized least squares (PGLS) of CEI on GC and GC^2, on the bacterial (bac120) and archaeal
# (ar53) trees of the Genome Taxonomy Database (GTDB release 232).
# One value per GTDB species representative (median of its proteomes; GTDB_phylogeny_data.py in
# runGC); the trees are pruned to the species with proteomes. Lambda is tested against 0 (no
# phylogenetic structure) with a likelihood-ratio test, the model under the null having lambda
# fixed close to 0.
# Input:  GTDB species values.csv; bac120_r232.tree and ar53_r232.tree (GTDB release 232)
# Output: Phylogenetic signal.csv (one row per model and quantity)

## 0. Installing necessary packages ##
#############################################################################################
if(!require('ape')){install.packages('ape')}
if(!require('phylolm')){install.packages('phylolm')}

## 0.1. Loading packages
########################
library(ape) # for "read.tree" and "keep.tip"
library(phylolm) # for "phylolm" (Pagel's lambda model)

###########################
## 1. Preparing the data ##
#############################################################################################
species_file <- "GTDB species values.csv"
if (!file.exists(species_file)){
  species_file <- "../Supplementary_materials/CAIR_supplementary_files/GTDB species values.csv" # the copy in the repository
}
tips <- read.csv(species_file, check.names = FALSE)
tips <- data.frame(rep = tips[["GTDB species representative"]], domain = tips$Domain, CEI = tips$CEI, GC = tips$GC,
                   CEI_null = tips[["CEI null"]], dCEI = tips[["Delta CEI"]])
tree_files <- c(Bacteria = "bac120_r232.tree", Archaea = "ar53_r232.tree")
results <- list()

###############################################
## 2. Pagel's lambda and PGLS for each domain ##
#############################################################################################
for (dom in c("Archaea", "Bacteria")) {
  tr <- read.tree(tree_files[[dom]])
  d <- tips[tips$domain == dom & tips$rep %in% tr$tip.label, ]
  tr <- keep.tip(tr, d$rep)
  tr$edge.length[tr$edge.length <= 0] <- 1e-6 # a few zero-length branches
  rownames(d) <- d$rep
  cat(dom, ": ", nrow(d), " species, ", Ntip(tr), " tips\n", sep = "")

  ## 2.1. Phylogenetic signal of each trait
  for (v in c("CEI", "GC", "CEI_null", "dCEI")) {
    f <- as.formula(paste(v, "~ 1"))
    m1 <- phylolm(f, data = d, phy = tr, model = "lambda")
    m0 <- phylolm(f, data = d, phy = tr, model = "lambda", lower.bound = 1e-7, upper.bound = 2e-7, starting.value = 1.5e-7) # lambda ~ 0
    lr <- 2 * (m1$logLik - m0$logLik)
    results[[paste(dom, v)]] <- c(n = nrow(d), lambda = m1$optpar, logLik = m1$logLik, logLik_lambda0 = m0$logLik, LR = lr,
                                  p = pchisq(lr, 1, lower.tail = FALSE))
    cat(sprintf("  %-8s lambda = %.3f  LR = %.1f  p = %.3g\n", v, m1$optpar, lr, pchisq(lr, 1, lower.tail = FALSE)))
  }

  ## 2.2. PGLS of CEI on GC and GC^2 (and the same model without phylogeny)
  ols <- lm(CEI ~ GC + I(GC^2), data = d)
  pg <- phylolm(CEI ~ GC + I(GC^2), data = d, phy = tr, model = "lambda")
  s <- summary(pg)
  results[[paste(dom, "PGLS CEI ~ GC + GC2")]] <- c(n = nrow(d), lambda = pg$optpar, R2 = pg$r.squared, adjR2 = pg$adj.r.squared,
                                                     coef_GC = unname(coef(pg)["GC"]), p_GC = s$coefficients["GC", "p.value"],
                                                     coef_GC2 = unname(coef(pg)["I(GC^2)"]), p_GC2 = s$coefficients["I(GC^2)", "p.value"],
                                                     OLS_R2 = summary(ols)$r.squared)
  cat(sprintf("  PGLS CEI ~ GC + GC^2: lambda = %.3f, R2 = %.3f (OLS R2 = %.3f); GC p = %.3g, GC^2 p = %.3g\n", pg$optpar, pg$r.squared,
              summary(ols)$r.squared, s$coefficients["GC", "p.value"], s$coefficients["I(GC^2)", "p.value"]))
}

#################################
## 3. Writing the results table ##
#############################################################################################
res <- do.call(rbind, lapply(names(results), function(k) data.frame(Model = k, Quantity = names(results[[k]]), Value = unname(results[[k]]))))
write.csv(res, "Phylogenetic signal.csv", row.names = FALSE)
cat("written Phylogenetic signal.csv\n")
