## CEI vs genome GC (rev 26, Figure 3) ##
#############################################################################################
# A: CEI (formerly CAIR) of each proteome against the GC content of its genome assembly, with the
#    CEI expected under random translation at that GC content (the null model of H2)
# B: GC-adjusted CEI (delta CEI = observed - expected) of the phyla with at least 50 proteomes
# Input:  GC null per proteome.csv and GC null curve.csv (Proteome_GC_null.py, or the copies in Supplementary_materials)
# Output: Fig3.svg
# Drawn with ggpubr like Figure 5 (6-pt sans, bold axis titles),
#   horizontal boxes in B so that the phylum names can be read, and a boxed legend below the panels
#   as in Figures 2, 4 and 5

## 0. Installing necessary packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('ggpubr')){install.packages('ggpubr')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(data.table) # for "fread"
library(ggpubr) # for "ggarrange", "theme_pubr" and "font" (loads ggplot2)
library(svglite) # for "svglite"

###########################
## 1. Preparing the data ##
#############################################################################################
gc_file <- "GC null per proteome.csv" # from runGC, otherwise the copy in the repository
if (!file.exists(gc_file)){
  gc_file <- "../Supplementary_materials/CAIR_supplementary_files/GC null per proteome.csv"
}
gc <- fread(gc_file)
gc <- gc[!is.na(gc$GC), ] # proteomes matched to an assembly
gc$Domain <- factor(gc$Domain, levels = c("Bacteria", "Archaea", "Eukaryota"))
gc <- gc[order(gc$Domain), ] # the smaller domains are drawn on top
curve_file <- "GC null curve.csv" # from runGC, otherwise the copy in the repository
if (!file.exists(curve_file)){
  curve_file <- "../Supplementary_materials/CAIR_supplementary_files/GC null curve.csv"
}
null_curve <- fread(curve_file)
ceiling_cei <- max(null_curve$`CEI null (code 11)`) # the highest CEI of the uniform-GC null
min_proteomes <- 50 # phyla shown in panel B
phyla <- gc[, .(n = .N, median_dcei = median(`Delta CEI`), Domain = Domain[1]), by = Second_hierarchy]
phyla <- phyla[!is.na(phyla$Second_hierarchy) & phyla$Second_hierarchy != "" & phyla$n >= min_proteomes, ] # proteomes without a second-rank lineage are not a phylum
phyla <- phyla[order(phyla$median_dcei), ] # lowest median at the bottom of panel B
short_name <- function (x) sub("^[Cc]andidate division ", "Cand. div. ", sub("^unclassified ", "uncl. ", sub("^Candidatus ", "Ca. ", x))) # shorter axis labels
gc_b <- gc[gc$Second_hierarchy %in% phyla$Second_hierarchy, ]
gc_b$phylum <- factor(short_name(gc_b$Second_hierarchy), levels = short_name(phyla$Second_hierarchy))

#####################################
## 2. Drawing the scatter and boxes ##
#############################################################################################
fs <- 6 # font size (as Figure 5)
fsl <- 6 # font size of the panel labels
font_family <- "sans"
line_size <- 0.25 # thin lines (mm) of the axes, boxes and reference lines
domain_col <- c(Bacteria = "#0f4d82", Archaea = "#4b5437", Eukaryota = "#9b4822") # colours of the domain names in Figure 2A
point_alpha <- c(Bacteria = 0.25, Archaea = 0.6, Eukaryota = 0.6) # the bacterial cloud is lighter, so the other domains stay visible
null_col <- "black" # the null model
ref_col <- "#a19997" # dashed reference lines

theme_fig3 <- theme_pubr(base_size = fs, base_family = font_family) +
  theme(axis.line = element_line(linewidth = line_size), axis.ticks = element_line(linewidth = line_size),
        axis.ticks.length = unit(1, "mm"), legend.position = "none", plot.margin = margin(8, 6, 2, 4))

## 2.1. CEI against genome GC
#############################
A <- ggplot(gc, aes(x = 100 * GC, y = CEI, colour = Domain, alpha = Domain)) +
  geom_point(size = 0.35, stroke = 0, shape = 16) +
  geom_hline(yintercept = ceiling_cei, linetype = "dashed", colour = ref_col, linewidth = line_size) +
  geom_line(data = null_curve, aes(x = 100 * GC, y = `CEI null (code 11)`), inherit.aes = FALSE, colour = null_col, linewidth = 0.45) +
  scale_colour_manual(values = domain_col) +
  scale_alpha_manual(values = point_alpha) +
  scale_x_continuous(breaks = seq(10, 80, by = 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.80, 0.95, by = 0.05), expand = c(0, 0)) +
  coord_cartesian(xlim = c(10, 80), ylim = c(0.80, 0.95)) + # axes end at labeled ticks (zoom only, so no data are dropped)
  labs(x = "Genome GC content (%)", y = "CEI") +
  theme_fig3 +
  font("xy", size = fs, family = font_family, face = 'bold') +
  font("xy.text", size = fs, family = font_family)

## 2.2. Delta CEI of the phyla
##############################
box_fill <- setNames(apply(col2rgb(domain_col) / 255 * 0.6 + 0.4, 2, function (v) rgb(v[1], v[2], v[3])), names(domain_col))
# (the domain colours at 60% on white, as opaque fills, so that the zero line stays behind the boxes)
B <- ggplot(gc_b, aes(x = `Delta CEI`, y = phylum, colour = Domain, fill = Domain)) +
  geom_vline(xintercept = 0, colour = ref_col, linewidth = 0.3) +
  # (delta CEI = 0, no deviation from the null: a solid grey line behind the boxes; the dashed line of the legend is panel A only)
  geom_boxplot(outlier.shape = NA, coef = 1.5, width = 0.62, linewidth = line_size, staplewidth = 0.4) +
  scale_colour_manual(values = domain_col) +
  scale_fill_manual(values = box_fill) +
  scale_x_continuous(breaks = seq(-0.04, 0.10, by = 0.02), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.04, 0.10)) + # beyond every whisker (zoom only: the boxes are computed from all proteomes)
  labs(x = "ΔCEI (observed − expected)", y = NULL) +
  theme_fig3 +
  font("xy", size = fs, family = font_family, face = 'bold') +
  font("xy.text", size = fs, family = font_family)

## 2.3. Legend box (as in Figures 2, 4 and 5)
#############################################
svglite(file = 'Fig3.svg', width = 7.08, height = 4.1) # 180 mm wide, like the other figures
par(ps = fs, family = font_family) # font metrics for spacing the legend items
legend_w <- 7.08 - 8 / 72 # width of the legend box (in), inside the plot margins
sw <- function (s, bold = FALSE) strwidth(s, units = "inches", font = ifelse(bold, 2, 1)) / legend_w # text width as a fraction of the box
items <- data.frame(kind = c("title", "point", "point", "point", "title", "line", "dashed"),
                    label = c("Domain:", names(domain_col), "Random translation:", "expected CEI at the genome GC content",
                              sprintf("highest expected CEI (%.3f)", ceiling_cei)),
                    col = c("black", domain_col, "black", null_col, ref_col), stringsAsFactors = FALSE)
key_w <- c(title = 0, point = 0.010, line = 0.030, dashed = 0.030) # width of the keys
x <- 0.015; items$x_key <- NA; items$x_text <- NA
for (i in seq_len(nrow(items))) { # items placed one after the other: key, label, space
  items$x_key[i] <- x
  items$x_text[i] <- x + key_w[items$kind[i]] + ifelse(items$kind[i] == "title", 0, 0.006)
  x <- items$x_text[i] + sw(items$label[i], items$kind[i] == "title") + ifelse(items$kind[i] == "title", 0.012, 0.03)
}
x0 <- (1 - (x - 0.03 + 0.015)) / 2 # the row of items is centred in the box
items$x_key <- items$x_key + x0; items$x_text <- items$x_text + x0
L <- ggplot() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = font_family) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), plot.margin = margin(2, 4, 4, 4)) +
  annotate("text", x = items$x_text, y = 0.5, label = items$label, fontface = ifelse(items$kind == "title", "bold", "plain"),
           size = fs / .pt, hjust = 0, family = font_family) +
  annotate("point", x = items$x_key[items$kind == "point"] + 0.005, y = 0.5, colour = items$col[items$kind == "point"], size = 1.6, shape = 16) +
  annotate("segment", x = items$x_key[items$kind == "line"], xend = items$x_key[items$kind == "line"] + key_w["line"], y = 0.5, yend = 0.5,
           colour = null_col, linewidth = 0.45) +
  annotate("segment", x = items$x_key[items$kind == "dashed"], xend = items$x_key[items$kind == "dashed"] + key_w["dashed"], y = 0.5, yend = 0.5,
           colour = ref_col, linewidth = line_size, linetype = "dashed")

fig <- ggarrange(ggarrange(A, B, ncol = 2, widths = c(1, 1.08), labels = c("A", "B"), hjust = 0, vjust = 1.2,
                           font.label = list(size = fsl, face = "bold", family = font_family)),
                 L, nrow = 2, heights = c(1, 0.085))
print(fig)
dev.off()
