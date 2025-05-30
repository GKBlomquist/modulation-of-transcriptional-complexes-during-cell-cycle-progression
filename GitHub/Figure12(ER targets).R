## RNA-seq vulcano plot (figure 12a) ## ----
# load data
Fulv_data <- read.delim("files/FULV_ER_target_genes.txt") %>%  # C1 = FULV/VEH, C2 = E2/VEH
  mutate(
    group_FULV = case_when( # subset genes according to p-value and log2fc
      log2FC.c1 > 0.5 & padj.c1 < 0.05 ~ "Induced",
      log2FC.c1 < -0.5 & padj.c1 < 0.05 ~ "Repressed",
      TRUE ~ "Insignificant"
    ))

# select top10 most significant (negative log2fc)
highlight_FULV <- Fulv_data %>% 
  filter(group_FULV == "Repressed") %>% 
  slice_min(padj.c1, n = 10) %>% 
  pull(Geneid)

# visualize
ggplot(Fulv_data, aes(x = log2FC.c1, y = -log10(padj.c1), color = group_FULV))+
  geom_point(size = 0.5)+
  geom_text_repel(data = Fulv_data %>% filter(Geneid %in% highlight_FULV), aes(label = Geneid), color = "#2A3D45", size = 3.5)+
  labs(x = "log2(FULV/VEH)", y = "-log10(p.adjusted)")+
  scale_color_manual(values = rev(c("#4E86AC", "lightgrey", "#BB3037")))+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

# grab ER targets
Fulv_genes <- Fulv_data %>% 
  filter(group_FULV == "Repressed") %>% 
  pull(Geneid) %>% 
  base::intersect(rownames(IPAX_seurat_sce))

## GAM smoothing (figure 12c) ## ----
# grab counts for all ER target genes
counts <- IPAX_seurat_sce@assays@data$counts[Fulv_genes,] # The IPAX_seurat_sce object is established in the "Figure9(cell cycle complexes).R" script along with the pseudotime scores. Make sure to run that script before this one

# grab pseudotime values for all cells
pseudotime <- slingPseudotime(IPAX_seurat_sce, na = FALSE)
cellWeights <- slingCurveWeights(IPAX_seurat_sce)

library(BiocParallel)
BPPARAM <- SnowParam(workers = 6, type = "SOCK") # optimize performance (works best for MAC, may need tweaking to run better on windows)

# Fit NB-GAM model to the selected ER targets
sce_model_ER_targets <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                               nknots = 7, verbose = TRUE, parallel=TRUE, BPPARAM = BPPARAM)

# Perform statistical test to check whether gene expression is constant across pseudotime within a lineage
assoRes_ER_targets <- associationTest(sce_model_ER_targets, contrastType = "consecutive", inverse = "eigen", lineages = F)

# Extract genes along lineage
genes_traj_ER_targets = assoRes_ER_targets %>%
  mutate(padjust = p.adjust(pvalue, method = "BH")) %>% 
  filter(padjust < 0.05) %>%
  arrange(desc(waldStat)) %>%
  top_n(waldStat, n = 300)

genes_of_interest = rownames(genes_traj_ER_targets)

#Get smoothers estimated counts by tradeSeq along a grid
Smooth_R1 <- predictSmooth(sce_model_ER_targets, gene = genes_of_interest, tidy=F)

## Calculate the row Z-score
zHM_meanSub <- (Smooth_R1 - rowMeans(Smooth_R1))
zHM_final <- (zHM_meanSub)/(matrixStats::rowSds(as.matrix(Smooth_R1)))

#Saturating it a bit. So all genes with a scaled rlog value greater than 2 should just be viewed as 2.
zHM_final[zHM_final>2] <- 2
# Everything lower than -2 should just be called -2
zHM_final[zHM_final< -2] <- -2

# Estimate linear order of smooth data
Smooth <- zHM_final[get_order(seriate(zHM_final, method="PCA_angle")),]


# Build heatmap
HM <- Heatmap(Smooth, 
              cluster_columns=F, 
              cluster_rows=T, 
              col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(1000), 
              show_column_names = F,
              show_row_names = T,
              km = 4,
              column_title = "Pseudo Time",
              column_title_side = "bottom")
HM <- draw(HM)

## GSEA analysis (figure 12d) ## ----
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(DOSE)

# extract gene clusters
rcl.list <- row_order(HM)

# import hallmarks
Hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

# perform enrichment analysis for each gene cluster
Hallmark_results <- list()
for (i in seq_along(rcl.list)) {
  Hallmark_results[[i]] <- enricher(gene = rownames(Smooth)[rcl.list[[i]]], 
                                    pvalueCutoff = 0.05, 
                                    pAdjustMethod = "BH",
                                    TERM2GENE = Hallmark_sets %>% dplyr::select(gs_name, gene_symbol))
}

# Combine the data frames and add an 'ID' column
combined_Hallmark_results <- bind_rows(lapply(seq_along(Hallmark_results), function(i) {
  df <- as.data.frame(Hallmark_results[[i]])
  df$ID <- clusters[i]  # Add the ID column denoting the list element
  return(df)
}))

combined_Hallmark_heatmap <- combined_Hallmark_results %>% 
  mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>% 
  dplyr::select(`-log10(p.adjust)`, Description, ID) %>% 
  pivot_wider(names_from = ID, values_from = `-log10(p.adjust)`) %>% 
  column_to_rownames("Description")

# remove NAs
combined_Hallmark_heatmap[is.na(combined_Hallmark_heatmap)] <- 0

pheatmap::pheatmap(combined_Hallmark_heatmap, color = colorRampPalette(c("white", "#2D7AB5","#2D7AB5"))(200), angle_col = 0, cluster_rows = F)


## ChIP average plots (figure 12e) ## ----
cluster_colors <- c("#b2df8a","#33a02c", "#a6cee3", "#1f78b4")

# ER
ER_Cluster1_hist <- read.delim("files/ChIP/ER/ER_Cluster1_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.ER_Cluster1_peaks.bed.hg38..d.ER_TagDir..hist.50..size.5000.,
                "Count" = ER_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 1")

ER_Cluster2_hist <- read.delim("files/ChIP/ER/ER_Cluster2_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.ER_Cluster2_peaks.bed.hg38..d.ER_TagDir..hist.50..size.5000.,
                "Count" = ER_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 2")

ER_Cluster3_hist <- read.delim("files/ChIP/ER/ER_Cluster3_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.ER_Cluster3_peaks.bed.hg38..d.ER_TagDir..hist.50..size.5000.,
                "Count" = ER_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 3")

ER_Cluster4_hist <- read.delim("files/ChIP/ER/ER_Cluster4_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.ER_Cluster4_peaks.bed.hg38..d.ER_TagDir..hist.50..size.5000.,
                "Count" = ER_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 4")

ER_hist_merged <- ER_Cluster1_hist %>% 
  rbind(ER_Cluster2_hist) %>% 
  rbind(ER_Cluster3_hist) %>% 
  rbind(ER_Cluster4_hist) %>% 
  add_column(Bait = "ER")

# MYC
MYC_Cluster1_hist <- read.delim("files/ChIP/MYC/MYC_Cluster1_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.MYC_Cluster1_peaks.bed.hg38..d.MYC_VEH_TagDir..hist.50..size.5000.,
                "Count" = MYC_VEH_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 1")

MYC_Cluster2_hist <- read.delim("files/ChIP/MYC/MYC_Cluster2_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.MYC_Cluster2_peaks.bed.hg38..d.MYC_VEH_TagDir..hist.50..size.5000.,
                "Count" = MYC_VEH_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 2")

MYC_Cluster3_hist <- read.delim("files/ChIP/MYC/MYC_Cluster3_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.MYC_Cluster3_peaks.bed.hg38..d.MYC_VEH_TagDir..hist.50..size.5000.,
                "Count" = MYC_VEH_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 3")

MYC_Cluster4_hist <- read.delim("files/ChIP/MYC/MYC_Cluster4_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.MYC_Cluster4_peaks.bed.hg38..d.MYC_VEH_TagDir..hist.50..size.5000.,
                "Count" = MYC_VEH_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 4")

MYC_hist_merged <- MYC_Cluster1_hist %>% 
  rbind(MYC_Cluster2_hist) %>% 
  rbind(MYC_Cluster3_hist) %>% 
  rbind(MYC_Cluster4_hist) %>% 
  add_column(Bait = "MYC")

# BRG1
BRG1_Cluster1_hist <- read.delim("files/ChIP/BRG1/BRG1_ER_Cluster1_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.BRG1_ER_Cluster1_peaks.bed.hg38..d.BRG1_TagDir..hist.50..size.5000.,
                "Count" = BRG1_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 1")

BRG1_Cluster2_hist <- read.delim("files/ChIP/BRG1/BRG1_ER_Cluster2_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.BRG1_ER_Cluster2_peaks.bed.hg38..d.BRG1_TagDir..hist.50..size.5000.,
                "Count" = BRG1_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 2")

BRG1_Cluster3_hist <- read.delim("files/ChIP/BRG1/BRG1_ER_Cluster3_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.BRG1_ER_Cluster3_peaks.bed.hg38..d.BRG1_TagDir..hist.50..size.5000.,
                "Count" = BRG1_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 3")

BRG1_Cluster4_hist <- read.delim("files/ChIP/BRG1/BRG1_ER_Cluster4_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.BRG1_ER_Cluster4_peaks.bed.hg38..d.BRG1_TagDir..hist.50..size.5000.,
                "Count" = BRG1_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 4")

BRG1_hist_merged <- BRG1_Cluster1_hist %>% 
  rbind(BRG1_Cluster2_hist) %>% 
  rbind(BRG1_Cluster3_hist) %>% 
  rbind(BRG1_Cluster4_hist) %>% 
  add_column(Bait = "BRG1")

# P300
P300_Cluster1_hist <- read.delim("files/ChIP/P300/P300_ER_Cluster1_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.P300_ER_Cluster1_peaks.bed.hg38..d.P300_TagDir..hist.50..size.5000.,
                "Count" = P300_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 1")

P300_Cluster2_hist <- read.delim("files/ChIP/P300/P300_ER_Cluster2_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.P300_ER_Cluster2_peaks.bed.hg38..d.P300_TagDir..hist.50..size.5000.,
                "Count" = P300_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 2")

P300_Cluster3_hist <- read.delim("files/ChIP/P300/P300_ER_Cluster3_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.P300_ER_Cluster3_peaks.bed.hg38..d.P300_TagDir..hist.50..size.5000.,
                "Count" = P300_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 3")

P300_Cluster4_hist <- read.delim("files/ChIP/P300/P300_ER_Cluster4_hist.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.P300_ER_Cluster4_peaks.bed.hg38..d.P300_TagDir..hist.50..size.5000.,
                "Count" = P300_TagDir.Coverage) %>% 
  add_column(Cluster = "Cluster 4")

P300_hist_merged <- P300_Cluster1_hist %>% 
  rbind(P300_Cluster2_hist) %>% 
  rbind(P300_Cluster3_hist) %>% 
  rbind(P300_Cluster4_hist) %>% 
  add_column(Bait = "P300")

## Merge data frames
ChIP_data_merged <- ER_hist_merged %>% 
  rbind(MYC_hist_merged) %>% 
  rbind(BRG1_hist_merged) %>% 
  rbind(P300_hist_merged)

# visualize
ggplot(ChIP_data_merged, aes(x = DFC, y = Count, color = Cluster))+
  facet_wrap(.~factor(Bait, levels = c("ER", "MYC", "BRG1", "P300")), nrow = 1, scales = "free")+
  geom_line()+
  labs(x = NULL, y = "Tag Count", color = NULL)+
  scale_color_manual(values = cluster_colors)+
  scale_x_continuous(labels = c("-2kb", "-1kb", "Center", "1kb", "2kb"), breaks = c(-2000, -1000, 0, 1000, 2000))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))
## MYC tag distribution (figure 12f) ## ----
MYC_Cluster1_annotated <- read.delim("files/ChIP/MYC/MYC_Cluster1_annotated.txt") %>% 
  dplyr::select(Annotation, MYC_VEH_TagDir.Tag.Count.in.given.bp..58211597.5.Total..normalization.factor...0.17..effective.total...10000000.) %>% 
  add_column(Cluster = "Cluster 1")

MYC_Cluster2_annotated <- read.delim("files/ChIP/MYC/MYC_Cluster2_annotated.txt") %>% 
  dplyr::select(Annotation, MYC_VEH_TagDir.Tag.Count.in.given.bp..58211597.5.Total..normalization.factor...0.17..effective.total...10000000.) %>% 
  add_column(Cluster = "Cluster 2")

MYC_Cluster3_annotated <- read.delim("files/ChIP/MYC/MYC_Cluster3_annotated.txt") %>% 
  dplyr::select(Annotation, MYC_VEH_TagDir.Tag.Count.in.given.bp..58211597.5.Total..normalization.factor...0.17..effective.total...10000000.) %>% 
  add_column(Cluster = "Cluster 3")

MYC_Cluster4_annotated <- read.delim("files/ChIP/MYC/MYC_Cluster4_annotated.txt") %>% 
  dplyr::select(Annotation, MYC_VEH_TagDir.Tag.Count.in.given.bp..58211597.5.Total..normalization.factor...0.17..effective.total...10000000.) %>% 
  add_column(Cluster = "Cluster 4")

MYC_annotated_merged <- MYC_Cluster1_annotated %>% 
  rbind(MYC_Cluster2_annotated) %>% 
  rbind(MYC_Cluster3_annotated) %>% 
  rbind(MYC_Cluster4_annotated) %>% 
  mutate(Annotation = gsub(" \\(.*", "", Annotation))

# visualize
ggplot(MYC_annotated_merged, aes(x = Cluster, y = MYC_VEH_TagDir.Tag.Count.in.given.bp..58211597.5.Total..normalization.factor...0.17..effective.total...10000000., fill = factor(Annotation, levels = c("promoter-TSS", "5' UTR", "exon", "intron", "3' UTR", "TTS", "non-coding", "Intergenic", "Other"))))+
  geom_col(position = "fill")+
  labs(x = NULL, y = "Tag Count", fill = NULL)+
  scale_fill_manual(values = brewer.pal(9, "RdBu"))+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9))

## DRIP-seq average plots (figure 12g) ## ----
Cluster1_hist <- read.delim("files/DRIP/Cluster1_hist.txt") %>%
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.Cluster1_DRIP_peaks.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Veh" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  add_column(Cluster = 1)

Cluster2_hist <- read.delim("files/DRIP/Cluster2_hist.txt") %>%
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.Cluster2_DRIP_peaks.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Veh" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  add_column(Cluster = 2)

Cluster3_hist <- read.delim("files/DRIP/Cluster3_hist.txt") %>%
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.Cluster3_DRIP_peaks.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Veh" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  add_column(Cluster = 3)

Cluster4_hist <- read.delim("files/DRIP/Cluster4_hist.txt") %>%
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.Cluster4_DRIP_peaks.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Veh" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  add_column(Cluster = 4)

# combine results
Cluster_hist_merged <- rbind(Cluster1_hist, Cluster2_hist) %>% 
  rbind(Cluster3_hist) %>% 
  rbind(Cluster4_hist) %>% 
  pivot_longer(Veh:E2, names_to = "Condition", values_to = "Count") %>% 
  mutate(Cluster = paste("Cluster", Cluster))

# visualize
ggplot(Cluster_hist_merged, aes(x = DFC, y = Count, color = Condition))+
  facet_wrap(.~Cluster, nrow = 1)+
  geom_line(show.legend = F)+
  scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2kb", "-1kb", "Center", "1kb", "2kb"))+
  scale_y_continuous(breaks = seq(0, 12, 3))+
  labs(x = NULL, y = "Tag Count")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major.y = element_line(linewidth = 0.2, linetype = "dashed", color = "darkgrey"))

## Gene cluster E2 expression (figure 12h) ## ----
# Extract cluster genes
Cluster_index <- row_order(HM)
cluster_matrices <- lapply(names(Cluster_index), function(cluster_name) {
  submatrix <- Smooth[Cluster_index[[cluster_name]], , drop = FALSE]  # Extract rows
  submatrix <- cbind(Cluster = cluster_name, submatrix)  # Add a new column for cluster
  return(submatrix)
})

gene_clusters <- do.call(rbind, cluster_matrices) %>% 
  data.frame() %>% 
  dplyr::select(Cluster) %>% 
  rownames_to_column("Gene")

Cluster1_genes <- gene_clusters %>% 
  filter(Cluster == 1) %>% 
  pull(Gene)

Cluster2_genes <- gene_clusters %>% 
  filter(Cluster == 2) %>% 
  pull(Gene)

Cluster3_genes <- gene_clusters %>% 
  filter(Cluster == 3) %>% 
  pull(Gene)

Cluster4_genes <- gene_clusters %>% 
  filter(Cluster == 4) %>% 
  pull(Gene)

# import RNA-seq data
MCF7_RNA_seq <- read.delim("files/FULV_ER_target_genes.txt") %>% 
  mutate(Cluster = case_when(
    Geneid %in% Cluster1_genes ~ "Cluster 1",
    Geneid %in% Cluster2_genes ~ "Cluster 2",
    Geneid %in% Cluster3_genes ~ "Cluster 3",
    Geneid %in% Cluster4_genes ~ "Cluster 4",
    TRUE ~ "None"
  )) %>% 
  mutate(
    E2_mean = rowMeans(across(matches("E2") & !matches("FULV_E2")))
  )

# visualize
ggplot(MCF7_RNA_seq %>% filter(Cluster != "None"), aes(x = Cluster, y = E2_mean, fill = Cluster))+
  geom_boxplot(outliers = F, show.legend = F)+
  scale_fill_manual(values = cluster_colors)+
  labs(x = NULL, y = "Normalized Count")+
  theme(axis.line = element_line(linewidth = 0.4),
        panel.border = element_rect(linewidth = 0.4))
