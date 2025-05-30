## libraries and functions ## ----
library(Seurat)
library(DoubletFinder)
library(presto)
library(tidyverse)

scale_by_percentile <- function(x) {
  p5 <- quantile(x, 0.05, na.rm = TRUE)
  p95 <- quantile(x, 0.95, na.rm = TRUE)
  x <- (x - p5) / (p95 - p5) # Scale to 0-1
  x[x < 0] <- 0 # Set values below 5th percentile to 0
  x[x > 1] <- 1 # Set values above 95th percentile to 1
  return(x)
}

## LOAD DATA ##----
## iPAX-seq R1
# Load counts
IPAXR1.data <- Read10X("/path/to/iPAXseq_R1_10x")

IPAXR1.PA <- as.sparse(read.delim(file = "/path/to/IPAXR1_CITE_counts", 
                                  header = T, 
                                  row.names = 1))

IPAXR1.PLA <- as.sparse(read.delim(file = "/path/to/IPAXR1_IPAX_counts", 
                                   header = T, 
                                   row.names = 1))

# Identify cells present in each assay
common_cells <- base::intersect(colnames(IPAXR1.PA), colnames(IPAXR1.data)) %>% base::intersect(colnames(IPAXR1.PLA))


IPAXR1.PA <- IPAXR1.PA[, common_cells]
IPAXR1.PLA <- IPAXR1.PLA[, common_cells]

# Scale PA counts using the function defined above
IPAXR1.PA.scaled <- t(apply(IPAXR1.PA, 1, scale_by_percentile))

PA_R1_assay <- CreateAssay5Object(counts = IPAXR1.PA, data = IPAXR1.PA.scaled)
PLA_R1_assay <- CreateAssay5Object(counts = IPAXR1.PLA)

IPAXR1 <- CreateSeuratObject(IPAXR1.data)
IPAXR1 <- subset(IPAXR1, cells = common_cells)

IPAXR1[["PA"]] <- PA_R1_assay
IPAXR1[["PLA"]] <- PLA_R1_assay

IPAXR1[["percent.mt"]] <- PercentageFeatureSet(IPAXR1, pattern = "^MT-")


## iPAX-seq R2 (repeat procedure)
IPAXR2.data <- Read10X("/path/to/iPAXseq_R2_10x")

IPAXR2.PA <- as.sparse(read.delim(file = "/path/to/IPAXR2_CITE_counts", 
                                  header = T, 
                                  row.names = 1))

IPAXR2.PLA <- as.sparse(read.delim(file = "/path/to/IPAXR2_IPAX_coounts", 
                                   header = T, 
                                   row.names = 1))

common_cells <- base::intersect(colnames(IPAXR2.PA), colnames(IPAXR2.data)) %>% base::intersect(colnames(IPAXR2.PLA))


IPAXR2.PA <- IPAXR2.PA[, common_cells]
IPAXR2.PLA <- IPAXR2.PLA[, common_cells]

IPAXR2.PA.scaled <- t(apply(IPAXR2.PA, 1, scale_by_percentile))

PA_R2_assay <- CreateAssay5Object(counts = IPAXR2.PA, data = IPAXR2.PA.scaled)
PLA_R2_assay <- CreateAssay5Object(counts = IPAXR2.PLA)

IPAXR2 <- CreateSeuratObject(IPAXR2.data)
IPAXR2 <- subset(IPAXR2, cells = common_cells)

IPAXR2[["PA"]] <- PA_R2_assay
IPAXR2[["PLA"]] <- PLA_R2_assay

IPAXR2[["percent.mt"]] <- PercentageFeatureSet(IPAXR2, pattern = "^MT-")

## create common object
# Label replicates
IPAXR1$replicate <- "Rep1"
IPAXR2$replicate <- "Rep2"

# Create merged object
IPAX_merged <- merge(IPAXR1, IPAXR2)


## Remove low quality cells (figure 7a) ## ----
IPAX_QC_data <- FetchData(IPAX_merged, vars = c("replicate", "nFeature_RNA", "nFeature_PA", "nFeature_PLA", "percent.mt")) # gather metadata

rep_colors <- c("#3FC08C", "#008A5A")


Feature_RNA_plot <- ggplot(IPAX_QC_data, aes(x = replicate, y = nFeature_RNA, fill = replicate))+
  geom_violin()+
  geom_hline(yintercept = 200, linetype = "dashed", color = "black")+
  scale_fill_manual(values = rep_colors)+
  labs(x = "", y = "RNA Features")+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.5))

percent.mt_plot <- ggplot(IPAX_QC_data, aes(x = replicate, y = percent.mt, fill = replicate))+
  geom_violin()+
  scale_fill_manual(values = rep_colors)+
  labs(x = "", y = "Mitochondrial RNA %")+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.5))

Feature_PLA_plot <- ggplot(IPAX_QC_data, aes(x = nFeature_PLA, fill = replicate))+
  geom_histogram()+
  geom_vline(xintercept = 9, linetype = "dashed")+
  scale_fill_manual(values = rep_colors)+
  labs(x = "PLA Features", y = "Count", fill = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.83, 0.85),
        panel.border = element_rect(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.5))

Feature_PA_plot <- ggplot(IPAX_QC_data, aes(x = nFeature_PA, fill = replicate))+
  geom_histogram()+
  scale_fill_manual(values = rep_colors)+
  labs(x = "PA Features", y = "Count", fill = NULL)+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.5))

Feature_RNA_plot + percent.mt_plot + Feature_PLA_plot + Feature_PA_plot


# filter out low quality cells
IPAX_merged <- subset(IPAX_merged, subset = nFeature_RNA > 200 | nFeature_PLA > 10)


## doublet identification ## ----
## Prior to doublet identification, the standard seurat workflow is applied
IPAX_merged <- NormalizeData(IPAX_merged)
IPAX_merged <- FindVariableFeatures(IPAX_merged)
IPAX_merged <- ScaleData(IPAX_merged)
IPAX_merged <- RunPCA(IPAX_merged, features = VariableFeatures(IPAX_merged))
ElbowPlot(IPAX_merged)
IPAX_merged <- FindNeighbors(IPAX_merged, dims = 1:15)
IPAX_merged <- FindClusters(IPAX_merged, resolution = 0.3)

## DoubletFinder
IPAX_merged <- JoinLayers(IPAX_merged)

## Find optimal neighborhood size
sweep.res.list <- paramSweep(IPAX_merged, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

ggplot(bcmvn, aes(x = pK, y = BCmetric)) +
  geom_point() 


# compute expected amount of dublets
prop <- modelHomotypic(IPAX_merged@meta.data$RNA_snn_res.0.3) # expected rate of homotypic doublets given the amount of cells in the identified clusters
nExp_poi <- round(0.075*ncol(IPAX_merged)) # expected amount of doublets with an assumed doublet formation rate of 0.075 multiplied by the total number of cells
nExp_poi.adj <- round(nExp_poi*(1-prop)) # adjust expected amount of doublets to account for low homotypic doublet detection

# identify dublets
IPAX_merged <- doubletFinder(IPAX_merged, PCs = 1:10, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) # find doublets using the estimated parameters

# Construct doublet UMAP (figure 7c)
IPAX_merged <- RunUMAP(IPAX_merged, reduction.name = "Doublet_UMAP", dims = 1:15) # Construct UMAP to visualize doublets

DimPlot(IPAX_merged, reduction = "Doublet_UMAP", group.by = "DF.classifications_0.25_0.29_941", cols = c(rep_colors[2], "grey"))+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(linewidth = 0.3))+
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)


# Inspect singlet vs doublet metadata (figure 7b)
metadata_doublets <- FetchData(IPAX_merged, vars = c("nFeature_RNA", "nCount_RNA", "nFeature_PLA", "nCount_PLA", "DF.classifications_0.25_0.29_941"))

ggplot(metadata_doublets, aes(x = nFeature_RNA, y = nCount_RNA, color = DF.classifications_0.25_0.29_941))+
  geom_point(size = 0.8)+
  scale_color_manual(values = c(rep_colors[2], "grey"))+
  geom_smooth(se = F, linetype = "dashed", color = "black", method = "lm", linewidth = 0.5)+
  labs(x = "RNA Feature Count", y = "RNA Count", title = "RNA")+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))+
  
  ggplot(metadata_doublets, aes(x = nFeature_PLA, y = nCount_PLA, color = DF.classifications_0.25_0.29_941))+
  geom_point(size = 0.8)+
  scale_color_manual(values = c(rep_colors[2], "grey"))+
  geom_smooth(se = F, linetype = "dashed", color = "black", method = "lm", linewidth = 0.5)+
  labs(x = "PLA Feature Count", y = "PLA Count", title = "PLA")+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

# remove identified doublets
IPAX_merged <- subset(IPAX_merged, subset = DF.classifications_0.25_0.29_941 == "Singlet")

## cell clustering and cell line identification
# Run clustering algorithm and UMAP once more on filtered cells
IPAX_merged <- FindNeighbors(IPAX_merged, dims = 1:15)
IPAX_merged <- FindClusters(IPAX_merged, resolution = 0.5, cluster.name = "Clusters")
IPAX_merged <- RunUMAP(IPAX_merged, reduction.name = "RNA_UMAP", dims = 1:15)

colors <- c("#6EA5CC", "#FF6C6C", 
            "#4E86AC","#8EC4EC","#BB3037", "#DC4F4F")

# Draw UMAP with clusters (figure 7d)
DimPlot(IPAX_merged, reduction = "RNA_UMAP", group.by = "Clusters", label = TRUE, cols = colors, label.box = FALSE, label.size = 6, label.color = "black")+
  theme(legend.position = "",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(linewidth = 0.3))+
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)

## Run processing steps on CITE (PA) assay
# Standard Workflow (Normalization performed using scale_by_percentile function above)
DefaultAssay(IPAX_merged) <- "PA"
IPAX_merged <- FindVariableFeatures(IPAX_merged)
IPAX_merged <- ScaleData(IPAX_merged)
IPAX_merged <- RunPCA(IPAX_merged)
ElbowPlot(IPAX_merged)
IPAX_merged <- FindNeighbors(IPAX_merged, dims = 1:9)
IPAX_merged <- FindClusters(IPAX_merged, resolution = 0.3)
IPAX_merged <- RunUMAP(IPAX_merged, reduction.name = "PA_UMAP", dims = 1:9)

## Run processing steps on IPAX (PLA) assay
# Standard Workflow
DefaultAssay(IPAX_merged) <- "PLA"
IPAX_merged <- NormalizeData(IPAX_merged, normalization.method = "CLR")
IPAX_merged <- FindVariableFeatures(IPAX_merged)
IPAX_merged <- ScaleData(IPAX_merged)
IPAX_merged <- RunPCA(IPAX_merged, features = VariableFeatures(IPAX_merged))
ElbowPlot(IPAX_merged)
IPAX_merged <- FindNeighbors(IPAX_merged, dims = 1:6)
IPAX_merged <- FindClusters(IPAX_merged, resolution = 0.3, cluster.name = "PLA_Clusters")
IPAX_merged <- RunUMAP(IPAX_merged, reduction.name = "PLA_UMAP", dims = 1:6)
IPAX_merged <- JoinLayers(IPAX_merged)

## cell line identification (figure 7e-g)
FeaturePlot(IPAX_merged, reduction = "PLA_UMAP", features = "ESR1")
FeaturePlot(IPAX_merged, reduction = "PLA_UMAP", features = "ER")
FeaturePlot(IPAX_merged, reduction = "PLA_UMAP", features = "ER:ER")

IPAX_merged$Cell_Line <- ifelse(IPAX_merged$Clusters %in% c("0", "2", "3"), "MCF7", "BT549") # add cell line identity to seurat object
DimPlot(IPAX_merged, reduction = "RNA_UMAP", group.by = "Cell_Line") # visualize cell line

## Protein complex DE analysis (figure h)
MCF7_markers_PLA <- FindMarkers(IPAX_merged, ident.1 = "MCF7", ident.2 = "BT549", group.by = "Cell_Line") # find differentially abundant complexes across the two cell lines
MCF7_complexes <- MCF7_markers_PLA %>% 
  rownames_to_column("Complex") %>% 
  filter(avg_log2FC > 0 & !str_detect(Complex, "FLAG") & !str_detect(Complex, "IgG")) %>%  # remove negative controls
  slice_min(p_val_adj, n = 10, with_ties = FALSE) %>% 
  pull(Complex) # pull the top10 most significant MCF7 complexes


BT549_markers_PLA <- FindMarkers(IPAX_merged, ident.1 = "BT549", ident.2 = "MCF7", group.by = "Cell_Line")
BT549_complexes <- BT549_markers_PLA %>% 
  rownames_to_column("Complex") %>% 
  filter(avg_log2FC > 0 & !str_detect(Complex, "FLAG") & !str_detect(Complex, "IgG")) %>%  # remove negative controls
  slice_min(p_val_adj, n = 10) %>% 
  pull(Complex) # pull top10 most significant BT549 complexes

# construct heatmap data
alt_heatmap <- FetchData(IPAX_merged, vars = c(MCF7_complexes, BT549_complexes, "Clusters")) %>% # pull abundance of identified complexes from seurat object along with cluster identity
  group_by(Clusters) %>% 
  summarise_all(.funs = mean) %>% # calculate mean abundance for each complex per cluster
  ungroup() %>% 
  mutate(Clusters = factor(Clusters, levels = c(0, 2, 3, 1, 4, 5))) %>% # convert cluster variable to factor allow for arrangment according to cell line
  arrange(Clusters) %>% 
  column_to_rownames("Clusters") %>% 
  as.matrix() %>% 
  t()


annotation_data <- data.frame(Cell_Line = rep(c("MCF7", "BT549"), each = 3))
rownames(annotation_data) <- c(0, 2, 3, 1, 4, 5)


annotation_colors <- list(Cell_Line = c("MCF7" = "#6EA5CC", "BT459" = "#BB3037"))

pheatmap::pheatmap(alt_heatmap, # draw heatmap
                   scale = "row", 
                   angle_col = 0, 
                   cluster_cols = F,
                   annotation_col = annotation_data, 
                   annotation_colors = annotation_colors, 
                   annotation_names_col = F, 
                   annotation_legend = F, 
                   gaps_col = 3,
                   color = colorRampPalette(c("#6EA5CC", "white", "#BB3037"))(100),
                   breaks = seq(-2, 2, length.out = 101),
                   cutree_rows = 1,
                   border_color = NA)


## Cell line gene set enrichment  (figure 7i) ## ----
# libraries
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(DOSE)

# Identify differentially expressed genes
Cell_Line_DE <- FindMarkers(IPAX_merged, ident.1 = "MCF7", ident.2 = "BT549", group.by = "Cell_Line")

# pull MCF7 and BT549 enriched genes
MCF7_genes <- Cell_Line_DE %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  rownames()

BT549_genes <- Cell_Line_DE %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  rownames()

genes_list <- list(MCF7_genes, BT459_genes)

# pull hallmark gene set from MSigDB
Hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Create empty list to store results in
Hallmark_results_cell_line <- list()

# run enrichment analysis for each gene list
for (i in seq_along(genes_list)) {
  Hallmark_results_cell_line[[i]] <- enricher(gene = genes_list[[i]], 
                                              pvalueCutoff = 0.05, 
                                              pAdjustMethod = "BH",
                                              TERM2GENE = Hallmark_sets %>% dplyr::select(gs_name, gene_symbol))
}



cell_lines <- c("MCF7\n(Cluster 0, 2, 3)", "BT549\n(1, 4, 5)")

# Combine the data frames and add cell line column
combined_Hallmark_results_cell_line <- bind_rows(lapply(seq_along(Hallmark_results_cell_line), function(i) {
  df <- as.data.frame(Hallmark_results_cell_line[[i]])
  df$ID <- cell_lines[i]
  return(df)
}))

# prepare heatmap data
combined_Hallmark_heatmap_cell_lines <- combined_Hallmark_results_cell_line %>% 
  mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>% # compute -log10(p.adj)
  dplyr::select(`-log10(p.adjust)`, Description, ID) %>% 
  pivot_wider(names_from = ID, values_from = `-log10(p.adjust)`) %>% 
  filter(Description %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING")) %>%  # select hallmarks of interest
  column_to_rownames("Description")

combined_Hallmark_heatmap_cell_lines[is.na(combined_Hallmark_heatmap_cell_lines)] <- 0 # remove NAs

# draw heatmap
pheatmap::pheatmap(combined_Hallmark_heatmap_cell_lines, color = colorRampPalette(c("white", "#2D7AB5"))(200), angle_col = 0, cluster_rows = F, cluster_cols = F, gaps_row = 4, gaps_col = 1)


# save seurat object
saveRDS(IPAX_merged, file = "files/IPAX_merged.rds")