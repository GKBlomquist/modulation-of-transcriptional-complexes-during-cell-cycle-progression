## libraries ## ----
library(Seurat)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(peco)
library(SingleCellExperiment)
library(tradeSeq)
library(slingshot)
library(seriation)
library(BiocParallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(scater)

## load data ## ----
IPAX_merged <- loadRDS("files/IPAX_merged_CC.rds")

## DE analysis ## ----
IPAX_MCF7 <- subset(IPAX_merged, subset = Cell_Line == "MCF7") # extract MCF7 cells

DefaultAssay(IPAX_MCF7) <- "PLA"

# perform differential expression analysis on IPAX counts
S_markers <- FindMarkers(IPAX_MCF7, ident.1 = "HC_S", ident.2 = c("HC_G1", "HC_G2M"), group.by = "HC_Phase")
G2M_markers <- FindMarkers(IPAX_MCF7, ident.1 = "HC_G2M", ident.2 = c("HC_S", "HC_G1"), group.by = "HC_Phase")
G1_markers <- FindMarkers(IPAX_MCF7, ident.1 = "HC_G1", ident.2 = c("HC_G2M", "HC_S"), group.by = "HC_Phase")

# select all significant complexes
S_markers_filtered <- S_markers %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(avg_log2FC) %>% 
  rownames_to_column("Complex")

G2M_markers_filtered <- G2M_markers %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(avg_log2FC) %>% 
  rownames_to_column("Complex")

G1_markers_filtered <- G1_markers %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(avg_log2FC) %>% 
  rownames_to_column("Complex")


# Volcano plots (figure 9a)
volcano_data_G1 <- G1_markers %>% 
  mutate(Group = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Downregulated",
    TRUE ~ "Insignificant"
  )) %>% 
  rownames_to_column("Complex")

ggplot(volcano_data_G1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group))+
  geom_point()+
  geom_text_repel(data = volcano_data_G1 %>% 
                    filter(Group != "Insignificant") %>% 
                    slice_min(p_val_adj, n = 15), aes(label = Complex), max.overlaps = 50)+
  scale_color_manual(values = c("#2D7AB5", "grey", "#BB3037"))+
  theme(legend.position = "")+
  labs(x = "log2(G1/S&G2M)", y = "-log10(p.adjusted)", title = "G1")+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

volcano_data_S <- S_markers %>% 
  mutate(Group = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Downregulated",
    TRUE ~ "Insignificant"
  )) %>% 
  rownames_to_column("Complex")

ggplot(volcano_data_S, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group))+
  geom_point()+
  geom_text_repel(data = volcano_data_S %>% 
                    filter(Group != "Insignificant") %>% 
                    slice_min(p_val_adj, n = 15), aes(label = Complex), max.overlaps = 50)+
  scale_color_manual(values = c("#2D7AB5", "grey", "#BB3037"))+
  theme(legend.position = "")+
  labs(x = "log2(S/G1&G2M)", y = "-log10(p.adjusted)", title = "S")+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

volcano_data_G2M <- G2M_markers %>% 
  mutate(Group = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Downregulated",
    TRUE ~ "Insignificant"
  )) %>% 
  rownames_to_column("Complex")

ggplot(volcano_data_G2M, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group))+
  geom_point()+
  geom_text_repel(data = volcano_data_G2M %>% 
                    filter(Group != "Insignificant") %>% 
                    slice_min(p_val_adj, n = 15), aes(label = Complex), max.overlaps = 50)+
  scale_color_manual(values = c("grey", "#BB3037"))+
  theme(legend.position = "")+
  labs(x = "log2(G2M/G1&S)", y = "-log10(p.adjusted)", title = "G2M")+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

## trend filtering ## ----
# Extract MCF7 cells
IPAX_MCF7 <- subset(IPAX_merged, subset = Cell_Line == "MCF7")
DefaultAssay(IPAX_MCF7) <- "PLA"
IPAX_MCF7 <- JoinLayers(IPAX_MCF7)
sce_IPAX_MCF7 <- as.SingleCellExperiment(IPAX_MCF7, assay = "PLA") # convert IPAX assay to sce object

# Normalize data
sce_IPAX_MCF7 <- data_transform_quantile(sce_IPAX_MCF7)

# Remove IgG and FLAG complexes
non_negative_complexes <- rownames(sce_IPAX_MCF7) %>% subset(., subset = !str_detect(., "IgG")) %>% subset(., subset = !str_detect(., "FLAG"))
PLA_yy_input <- assay(sce_IPAX_MCF7, "cpm_quantNormed")[non_negative_complexes,]

# Get theta predictions
PLA_theta_predict <- FetchData(IPAX_MCF7, vars = c("CC_Prediction")) %>% 
  rownames_to_column("Barcode") %>% 
  pull(CC_Prediction, name = Barcode)

# Estimate Trend Filters (!intensive!)
all_complexes_trend <- fit_cyclical_many(Y = PLA_yy_input, theta = PLA_theta_predict)


# Cell cycle variance plot (figure 9b left)
CC_variance <- all_complexes_trend$pve %>% 
  data.frame() %>% 
  rownames_to_column("Complex") %>% 
  arrange(pve) %>% 
  add_column(Ranked_pve = 1:nrow(.))

ggplot(CC_variance, aes(x = Ranked_pve, y = pve))+
  geom_point()+
  geom_text_repel(data = CC_variance %>% slice_max(pve, n = 9), aes(label = Complex), max.overlaps = 30, nudge_x = -40, box.padding = 0.4, nudge_y = 0.0012)

# Visualize trends
range <- seq(0, 2*pi, length.out = 100)

all_trends <- data.frame(
  x = range,
  do.call(cbind, lapply(non_negative_complexes, function(f) all_complexes_trend$cellcycle_function[[f]](range)))
)
colnames(all_trends) <- c("x", non_negative_complexes)

Top9_complexes <- CC_variance %>% slice_max(pve, n = 9) %>% pull(Complex)

# Heatmap (figure 9c)
heatmap_data_trends <- all_trends[,Top9_complexes] %>% 
  base::as.matrix() %>% 
  t()
colnames(heatmap_data) <- seq(0, 2*pi, length.out = 100) %>% as.character()

pheatmap(heatmap_data_trends,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         border_color = "grey",
         labels_col = c("G1", "S", "G2M"),angle_col = 0)

## Trend filtering BT549 
# Extract BT549 cells
IPAX_BT459 <- subset(IPAX_merged, subset = Cell_Line == "BT549")
DefaultAssay(IPAX_BT549) <- "PLA"
IPAX_BT549 <- JoinLayers(IPAX_BT549)
sce_IPAX_BT549 <- as.SingleCellExperiment(IPAX_BT549, assay = "PLA")

# Normalize data
sce_IPAX_BT549 <- data_transform_quantile(sce_IPAX_BT549)


PLA_yy_input_BT549 <- assay(sce_IPAX_BT549, "cpm_quantNormed")[non_negative_complexes,]

# Get theta predictions
PLA_theta_predict_BT549 <- FetchData(IPAX_BT549, vars = c("CC_Prediction")) %>% 
  rownames_to_column("Barcode") %>% 
  pull(CC_Prediction, name = Barcode)

# Estimate Trend Filters (!intensive!)
all_complexes_trend_BT549 <- fit_cyclical_many(Y = PLA_yy_input_BT549, theta = PLA_theta_predict_BT549)

# Visualize CC variance (figure 9b, right)
CC_variance_BT549 <- all_complexes_trend_BT549$pve %>% 
  data.frame() %>% 
  rownames_to_column("Complex") %>% 
  arrange(pve) %>% 
  add_column(Ranked_pve = 1:nrow(.))

ggplot(CC_variance_BT549, aes(x = Ranked_pve, y = pve))+
  geom_point()+
  geom_text_repel(data = CC_variance_BT459 %>% slice_max(pve, n = 9), aes(label = Complex), max.overlaps = 30, nudge_x = -40, box.padding = 0.4, nudge_y = 0.0012)

## tradeseq workflow ## ----
IPAX_seurat = subset(IPAX_merged, subset = Cell_Line == "MCF7") %>% JoinLayers() # extract MCF7 cells

#Convert dataset into a SingleCellExperiment Object
IPAX_seurat_sce = as.SingleCellExperiment(IPAX_seurat, assay = "RNA")
reducedDims(IPAX_seurat_sce)$UMAP <- IPAX_seurat@reductions$RNA_UMAP@cell.embeddings ## attach RNA_UMAP to sce data

# Run Slingshot analysis
IPAX_seurat_sce = slingshot(IPAX_seurat_sce, clusterLabels = "Phase", reducedDim = "RNA_UMAP", dist.method = 'mnn', approx_points = FALSE, start.clus = "G1", end.clus = "G2M", extend = "n")

## Fit data to a NB-GAM model
IPAX_seurat_sce_pla = as.SingleCellExperiment(IPAX_seurat, assay = "PLA")

biological_complexes <- rownames(IPAX_seurat_sce_pla@assays@data$counts) %>% subset(.,!str_detect(., "IgG") & !str_detect(.,"FLAG")) ## remove unwanted complexes to relieve computational load

counts_RNA = IPAX_seurat_sce@assays@data$counts[c("CDKN1A", "MCM4", "CDK1"),] ## extract marker genes
counts_PLA = IPAX_seurat_sce_pla@assays@data$counts ## extract complexes
counts = rbind(counts_RNA, counts_PLA)

pseudotime <- slingPseudotime(IPAX_seurat_sce, na = FALSE)
cellWeights <- slingCurveWeights(IPAX_seurat_sce)

BPPARAM <- BiocParallel::SnowParam(workers = 6, type = "SOCK") # optimize performance (works well on MAC but may need tweaking to work optimally on windows)

#evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, nGenes = 100) # run this function to assess the optimal amount of knots (K). This line is commented out as it is very computationally expensive

sce_model <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                    nknots = 7, verbose = TRUE, parallel=TRUE, BPPARAM = BPPARAM) # Fit GAM model

# Perform statistical test to check whether gene expression is constant across pseudotime
assoRes <- associationTest(sce_model, contrastType = "consecutive", inverse = "eigen", lineages = F)

# Extract genes along lineage
genes_traj = assoRes %>%
  mutate(padjust = p.adjust(pvalue, method = "BH")) %>% # calculate adjusted p-value
  filter(padjust < 0.05) %>%
  arrange(desc(waldStat)) %>%
  top_n(waldStat, n = 300)

genes_of_interest = rownames(genes_traj)

#Get smoothers estimated counts by tradeSeq along a grid
Smooth_R1 <- predictSmooth(sce_model, gene = genes_of_interest, tidy=F)

# save smoothers for later use
write.csv(Smooth_R1, file = "files/complex_splines.csv")

## Calculate the row Z-score
zHM_meanSub <- (Smooth_R1 - rowMeans(Smooth_R1))
zHM_final <- (zHM_meanSub)/(matrixStats::rowSds(as.matrix(Smooth_R1)))

# Values above or below 2/-2 are transformed to saturate the final heatmap
zHM_final[zHM_final>2] <- 2
zHM_final[zHM_final< -2] <- -2

# Estimate linear order of smooth data
Smooth <- zHM_final[get_order(seriate(zHM_final, method="PCA_angle")),]

# Build heatmap (figure 9d)
HM <- Heatmap(Smooth[rownames(Smooth) %in% biological_complexes,], 
              cluster_columns=F, 
              cluster_rows=T, 
              col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(1000), 
              km = 4,
              show_column_names = F
)
HM <- draw(HM)

# build heatmap for phase marker genes
HM_refgenes <- Heatmap(Smooth[c("CDK1", "MCM4", "CDKN1A"),],
                       cluster_columns = F,
                       cluster_rows = T,
                       km = 2,row_gap = unit(0, "cm"),
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(1000), 
                       column_title = "Pseudo Time",
                       column_title_side = "bottom",
                       show_column_names = F)
HM_refgenes <- draw(HM_refgenes)

