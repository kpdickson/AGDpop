# ==============================================================================
# AGDpop Analysis Pipeline (R)
# Description: Population genetics, co-evolution, and network analysis 
#              of Paramoeba host and Perkinsela symbiont AmpSeq.
# ==============================================================================

# --- DEPENDENCIES ---
library(vcfR)
library(adegenet)
library(poppr)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(paco)
library(igraph)
library(ggraph)
library(SNPRelate)
library(geosphere)
library(gridExtra)
library(ape)
library(StAMPP)

# ==============================================================================
# --- SETUP & USER INPUTS ---
# ==============================================================================
RUN_PREFIX <- "AGD" 
OUT_DIR <- "AGDout/"

para_vcf_path <- "AGD_VCFs/para-nuc.vcf.gz"
perk_vcf_path <- "AGD_VCFs/perk-nuc.vcf.gz"
metadata_path <- "AGDmetadata.csv"

# Check input files exist before proceeding
if (!file.exists(para_vcf_path)) stop("Host VCF not found at: ", para_vcf_path)
if (!file.exists(perk_vcf_path)) stop("Symbiont VCF not found at: ", perk_vcf_path)
if (!file.exists(metadata_path)) stop("Metadata not found at: ", metadata_path)

dir.create(OUT_DIR, showWarnings = FALSE)

# ==============================================================================
# --- DATA IMPORT & FORMATTING ---
# ==============================================================================

# 1. Import & Assign Populations
import_and_clean <- function(vcf_obj) {
  gl <- vcfR2genlight(vcf_obj)
  pop(gl) <- ifelse(grepl("^N", indNames(gl)), "Norway", "Scotland")
  ploidy(gl) <- 2
  return(gl)
}

para_gl <- import_and_clean(read.vcfR(para_vcf_path, verbose = FALSE))
perk_gl <- import_and_clean(read.vcfR(perk_vcf_path, verbose = FALSE))

# 2. Convert to Genclone Objects
build_genclone <- function(gl) {
  mat <- as.matrix(gl)
  char_mat <- matrix(NA_character_, nrow = nrow(mat), ncol = ncol(mat))
  
  char_mat[mat == 0] <- "1/1"
  char_mat[mat == 1] <- "1/2"
  char_mat[mat == 2] <- "2/2"
  rownames(char_mat) <- rownames(mat)
  colnames(char_mat) <- colnames(mat)
  
  gc <- as.genclone(df2genind(char_mat, sep = "/", ploidy = 2, type = "codom"))
  pop(gc) <- pop(gl)
  return(gc)
}

para_gc <- build_genclone(para_gl)
mlg.filter(para_gc, distance = bitwise.dist) <- 0.03

perk_gc <- build_genclone(perk_gl)
mlg.filter(perk_gc, distance = bitwise.dist) <- 0.03

# 3. Import Metadata
meta <- read.csv(metadata_path, stringsAsFactors = FALSE) %>% 
  dplyr::rename(Sample_Type = Sample, Sample = ID, Lon = Long) %>% 
  dplyr::filter(Sample != "")

# ==============================================================================
# --- PART 1: STATS & SUMMARY TABLE ---
# ==============================================================================
cat("\n--- Running Core Stats ---\n")
stats_out <- list()

# A. Index of Association (Ia)
set.seed(123)
ia_para <- poppr::ia(para_gc, sample = 999, quiet = TRUE, valuereturn = TRUE)
set.seed(123)
ia_perk <- poppr::ia(perk_gc, sample = 999, quiet = TRUE, valuereturn = TRUE)

stats_out$Host_Ia_p <- ia_para$index["p.rD"]
stats_out$Host_rbarD <- ia_para$index["rbarD"]
stats_out$Symb_Ia_p <- ia_perk$index["p.rD"]
stats_out$Symb_rbarD <- ia_perk$index["rbarD"]

# B. Pairwise Fst
set.seed(123)
fst_para <- stamppFst(stamppConvert(para_gl, type = "genlight"), nboots = 1000, percent = 95)
set.seed(123)
fst_perk <- stamppFst(stamppConvert(perk_gl, type = "genlight"), nboots = 1000, percent = 95)

stats_out$Host_Fst <- fst_para$Fsts[2,1]
stats_out$Host_Fst_p <- fst_para$Pvalues[2,1]
stats_out$Symb_Fst <- fst_perk$Fsts[2,1]
stats_out$Symb_Fst_p <- fst_perk$Pvalues[2,1]

# C. Co-evolution Data Prep (Matching common individuals)
common_inds <- intersect(indNames(para_gl), indNames(perk_gl))
para_sub <- para_gl[match(common_inds, indNames(para_gl)), ]
perk_sub <- perk_gl[match(common_inds, indNames(perk_gl)), ]

dist_para <- bitwise.dist(para_sub)
dist_perk <- bitwise.dist(perk_sub)

# D. PACo & ParaFit
assoc_mat <- diag(length(common_inds))
rownames(assoc_mat) <- common_inds
colnames(assoc_mat) <- common_inds

paco_data <- add_pcoord(prepare_paco_data(H = as.matrix(dist_para), P = as.matrix(dist_perk), HP = assoc_mat), correction = "cailliez")

set.seed(123)
paco_res <- PACo(paco_data, nperm = 1000, seed = 123, method = "r0", symmetric = FALSE, shuffled = TRUE)
stats_out$PACo_m2 <- paco_res$gof$ss
stats_out$PACo_p <- paco_res$gof$p

set.seed(123)
parafit_res <- ape::parafit(dist_para, dist_perk, assoc_mat, nperm = 999, test.links = TRUE, correction = "cailliez")
stats_out$ParaFit_Global <- parafit_res$ParaFitGlobal
stats_out$ParaFit_p <- parafit_res$p.global

# E. Export Summary Table
summary_df <- data.frame(
  Metric = names(stats_out),
  Value = as.numeric(stats_out)
) %>% 
  mutate(Value = round(Value, 4))

summary_file <- paste0(OUT_DIR, RUN_PREFIX, "_Stat-Summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Statistical analyses complete. Saved to:", summary_file, "\n")

# ==============================================================================
# --- PART 2: FIGURES ---
# ==============================================================================
cat("\n--- Generating Figures ---\n")

# --- DAPC ---
# optim.a.score(dapc(para_gl, pop(para_gl), n.pca = 20, n.da = 1))
# optim.a.score(dapc(perk_gl, pop(perk_gl), n.pca = 20, n.da = 1))
set.seed(123)
dapc_para <- dapc(para_gl, pop(para_gl), n.pca = 1, n.da = 1) 
set.seed(123)
dapc_perk <- dapc(perk_gl, pop(perk_gl), n.pca = 5, n.da = 1) 

pdf(paste0(OUT_DIR, RUN_PREFIX, "_DAPC.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), mgp = c(2, 0.7, 0))
scatter(dapc_para, col = c("#E41A1C", "#377EB8"), bg = "white", solid = 0.7, main = "Paramoeba Host", scree.da = FALSE, legend = TRUE)
scatter(dapc_perk, col = c("#E41A1C", "#377EB8"), bg = "white", solid = 0.7, main = "Perkinsela Symbiont", scree.da = FALSE, legend = TRUE)
invisible(dev.off())

# --- CLUSTER ASSIGNMENT (ADMIXTURE BARPLOT) ---
build_admixture_df <- function(dapc_obj, pop_labels, dataset_name) {
  df <- as.data.frame(dapc_obj$posterior)
  df$Sample <- rownames(df)
  df$Original_Pop <- pop_labels
  
  df_long <- tidyr::pivot_longer(df, cols = c("Norway", "Scotland"), names_to = "Assigned_Cluster", values_to = "Probability")
  df_long$Dataset <- dataset_name
  return(df_long)
}

admix_combined <- dplyr::bind_rows(
  build_admixture_df(dapc_para, pop(para_gl), "Paramoeba Host"), 
  build_admixture_df(dapc_perk, pop(perk_gl), "Perkinsela Symbiont")
)

admix_combined$Original_Pop <- dplyr::recode(admix_combined$Original_Pop, "Norway" = "Norwegian Isolates", "Scotland" = "Scottish Isolates")
admix_combined$Sample <- factor(admix_combined$Sample, levels = stringr::str_sort(unique(admix_combined$Sample), numeric = TRUE))

p_admix <- ggplot(admix_combined, aes(x = Sample, y = Probability)) + 
  geom_col(aes(fill = Assigned_Cluster), width = 1) +
  facet_grid(Dataset ~ Original_Pop, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = c("Norway" = "#E06666", "Scotland" = "#6FA8DC")) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) + 
  theme_minimal() +
  labs(title = "Cluster Membership Probability", x = NULL, y = "Membership Probability", fill = "Genetic Cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), 
        strip.text = element_text(face = "bold", size = 9), 
        strip.placement = "outside")

ggsave(filename = paste0(OUT_DIR, RUN_PREFIX, "_Admix.pdf"), plot = p_admix, width = 10, height = 4, bg = "white")

# --- ISOLATION BY DISTANCE (IBD) ---
run_ibd_plot <- function(gl, metadata, title) {
  shared <- intersect(indNames(gl), metadata$Sample)
  gl_sub <- gl[match(shared, indNames(gl)), ]
  meta_sub <- metadata[match(shared, metadata$Sample), ]
  
  gen_dist <- bitwise.dist(gl_sub)
  geo_dist_km <- as.dist(distm(as.matrix(meta_sub[, c("Lon", "Lat")]), fun = distGeo) / 1000) 
  
  set.seed(123)
  ibd_mantel <- vegan::mantel(gen_dist, geo_dist_km, permutations = 999)
  
  df_plot <- data.frame(Geo_Distance_km = as.numeric(geo_dist_km), Genetic_Distance = as.numeric(gen_dist))
  
  p <- ggplot(df_plot, aes(x = Geo_Distance_km, y = Genetic_Distance)) + 
    geom_point(alpha = 0.2, color = "#377EB8") +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") + 
    theme_bw() + 
    labs(title = paste("Isolation by Distance:", title), 
         subtitle = paste("Mantel r =", round(ibd_mantel$statistic, 3), "| p-value =", round(ibd_mantel$signif, 3)))
  
  ggsave(filename = paste0(OUT_DIR, RUN_PREFIX, "_IBD_", gsub(" ", "_", title), ".pdf"), plot = p, width = 6, height = 3, bg = "white")
}

run_ibd_plot(para_gl, meta, "Paramoeba Host")
run_ibd_plot(perk_gl, meta, "Perkinsela Symbiont")

# --- MINIMUM SPANNING NETWORKS (MSN) ---
pdf(paste0(OUT_DIR, RUN_PREFIX, "_MSN.pdf"), width = 10, height = 7)

# Host MSN
set.seed(42) 
para_msn <- poppr.msn(para_gc, bitwise.dist(para_gc), showplot = FALSE, include.ties = TRUE)
clean_labels_para <- gsub("MLG\\.", "", V(para_msn$graph)$name)
para_widths <- (1 / E(para_msn$graph)$weight)
para_widths <- (para_widths / max(para_widths)) * 3 

plot_poppr_msn(para_gc, para_msn, palette = c("#E41A1C", "#377EB8"), 
               main = paste("Paramoeba Host\n(", mlg(para_gc), "MLLs)"), 
               gscale = TRUE, pop.leg = FALSE, mlg = FALSE, 
               wscale = FALSE, edge.width = para_widths, gadj = 10, 
               vertex.frame.color = NA, vertex.label = clean_labels_para, 
               vertex.label.cex = 0.8, vertex.label.color = "black", vertex.label.dist = 0)

# Symbiont MSN
set.seed(42)
perk_msn <- poppr.msn(perk_gc, bitwise.dist(perk_gc), showplot = FALSE, include.ties = TRUE)
clean_labels_perk <- gsub("MLG\\.", "", V(perk_msn$graph)$name)
perk_widths <- (1 / E(perk_msn$graph)$weight)
perk_widths <- (perk_widths / max(perk_widths)) * 3 

plot_poppr_msn(perk_gc, perk_msn, palette = c("#E41A1C", "#377EB8"), 
               main = paste("Perkinsela Symbiont\n(", mlg(perk_gc), "MLLs)"), 
               gscale = TRUE, pop.leg = FALSE, mlg = FALSE, 
               wscale = FALSE, edge.width = perk_widths, gadj = 10, 
               vertex.frame.color = NA, vertex.label = clean_labels_perk,
               vertex.label.cex = 0.8, vertex.label.color = "black", vertex.label.dist = 0)
invisible(dev.off())

# --- PAIRWISE LINKAGE DISEQUILIBRIUM (LD) ---
run_pairwise_ld <- function(vcf_path, gds_filename, org_name, plot_color) {
  gds_path <- paste0(OUT_DIR, gds_filename)
  snpgdsVCF2GDS(vcf_path, gds_path, method = "biallelic.only", verbose = FALSE)
  genofile <- snpgdsOpen(gds_path)
  
  ld_matrix <- snpgdsLDMat(genofile, method = "corr", slide = -1, verbose = FALSE)
  r2_values <- ld_matrix$LD^2
  r2_values <- r2_values[upper.tri(r2_values)] 
  
  total_pairs <- length(r2_values[!is.na(r2_values)])
  percent_in_LD <- (sum(r2_values > 0.2, na.rm = TRUE) / total_pairs) * 100
  
  plot(density(r2_values, na.rm = TRUE), main = org_name, 
       xlab = expression(paste("Pairwise Correlation (", r^2, ")")), 
       ylab = "Density", lwd = 2, col = plot_color, xlim = c(-0.05, 1.1)) 
  abline(v = 0.2, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = paste(round(percent_in_LD, 2), "% of loci > 0.2"), col = "red", lty = 2, bty = "n", cex = 0.9)
  
  snpgdsClose(genofile)
  file.remove(gds_path)
}

pdf(paste0(OUT_DIR, RUN_PREFIX, "_LD.pdf"), width = 8, height = 3)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))
run_pairwise_ld(para_vcf_path, "para_temp.gds", "Paramoeba Host", "#377EB8")
run_pairwise_ld(perk_vcf_path, "perk_temp.gds", "Perkinsela Symbiont", "#4DAF4A")
invisible(dev.off())

# --- BIPARTITE NETWORK ---
mlg_df <- data.frame(
  Sample = common_inds,
  Host_MLG = paste0("Host_", mlg.vector(para_gc[common_inds, ])),
  Symbiont_MLG = paste0("Symb_", mlg.vector(perk_gc[common_inds, ])),
  Country = as.character(pop(para_sub)) 
)

edges_weighted <- mlg_df %>% 
  group_by(Host_MLG, Symbiont_MLG, Country) %>% 
  summarise(Weight = n(), .groups = "drop") %>% 
  select(from = Host_MLG, to = Symbiont_MLG, Weight, Country) 

g <- graph_from_data_frame(edges_weighted, directed = FALSE)
V(g)$type <- grepl("Host", V(g)$name) 

lo <- create_layout(g, layout = 'bipartite')
original_x <- lo$x
lo$x <- ifelse(lo$type, 0, 1)
lo$y <- original_x
lo$clean_label <- gsub("Host_|Symb_", "", lo$name)

p_network <- ggraph(lo) + 
  geom_edge_link(aes(edge_width = Weight, edge_colour = Country), alpha = 0.6) +
  scale_edge_width_continuous(range = c(0.4, 2.0)) + 
  geom_node_point(aes(fill = type), shape = 21, size = 3.5, color = "black") +
  geom_node_text(aes(label = clean_label, hjust = ifelse(type, 1.5, -0.5)), size = 3.5, fontface="bold") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
  theme_void() +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) + 
  scale_edge_color_manual(values = c("Norway" = "#E41A1C", "Scotland" = "#377EB8")) +
  theme(legend.position = "bottom", plot.margin = margin(20, 40, 20, 40)) + 
  guides(edge_width = "none")

ggsave(filename = paste0(OUT_DIR, RUN_PREFIX, "_BipartNet.pdf"), plot = p_network, width = 6, height = 8, bg = "white")


# ==============================================================================
# --- SUPPLEMENTARY FIGURES ---
# ==============================================================================

# --- INDEX OF ASSOCIATION ---
p_ia_para <- plot(ia_para, index = "rbarD") + 
  ggtitle(paste("Population: Total\nN:", nInd(para_gl), "\nData: Paramoeba Host\nPermutations: 999"))
p_ia_perk <- plot(ia_perk, index = "rbarD") + 
  ggtitle(paste("Population: Total\nN:", nInd(perk_gl), "\nData: Perkinsela Symbiont\nPermutations: 999"))

pdf(paste0(OUT_DIR, RUN_PREFIX, "_Supp_Ia.pdf"), width = 12, height = 6)
grid.arrange(p_ia_para, p_ia_perk, ncol = 2)
invisible(dev.off())

# --- PACO BELL CURVE ---
pdf(paste0(OUT_DIR, RUN_PREFIX, "_Supp_PACobell.pdf"), width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(density(paco_res$shuffled), main = "Global PACo Significance", 
     xlab = expression(paste("Procrustes Sum of Squared Residuals (", m^2, ")")), 
     ylab = "Frequency", col = "darkgrey", lwd = 2)
polygon(density(paco_res$shuffled), col = adjustcolor("grey", alpha.f = 0.3), border = "darkgrey")
abline(v = paco_res$gof$ss, col = "#E41A1C", lwd = 3, lty = 1)

p_text <- ifelse(paco_res$gof$p == 0, "P < 0.001", paste("P =", round(paco_res$gof$p, 3)))
legend("topright", legend = c("Null Distribution", paste("Observed Error (m2 =", round(paco_res$gof$ss, 2), ")\n", p_text)), 
       col = c("darkgrey", "#E41A1C"), lwd = c(2, 3), bty = "n", cex = 0.9)
invisible(dev.off())

# --- PACO RESIDUALS ---
res <- residuals_paco(paco_res$proc)
res_sorted <- sort(res, decreasing = TRUE)
raw_sample_ids <- sapply(strsplit(names(res_sorted), "-"), `[`, 1)
matched_mlgs <- mlg_df[match(raw_sample_ids, mlg_df$Sample), ]
names(res_sorted) <- paste0(raw_sample_ids, " (", gsub("Host_", "H", matched_mlgs$Host_MLG), "-", gsub("Symb_", "S", matched_mlgs$Symbiont_MLG), ")")

threshold <- quantile(res_sorted, 0.75)

pdf(paste0(OUT_DIR, RUN_PREFIX, "_PACo_Residuals.pdf"), width = 10, height = 7)
par(mar = c(8, 5, 4, 2) + 0.1) 
barplot(res_sorted, las = 2, cex.names = 0.65, col = ifelse(res_sorted > threshold, "#E69F00", "#999999"), 
        border = NA, ylim = c(0, max(res_sorted) * 1.05), ylab = "Procrustes Residual")
abline(h = threshold, lty = 2, col = "black", lwd = 2)
par(mar = c(5, 4, 4, 2) + 0.1) 
invisible(dev.off())

cat("\n--- Pipeline Complete ---\n")