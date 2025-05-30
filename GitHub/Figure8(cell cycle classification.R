## libraries and functions ## ----
library(Seurat)
library(tidyverse)

## load seurat object
IPAX_merged <- readRDS("files/IPAX_merged.rds")

## Categorical phase assignment (seurat) ## ----
## Cell Cycle Classification (Seurat) ----
DefaultAssay(IPAX_merged) <- "RNA"

# Cell Cycle Scoring
s.genes <- cc.genes.updated.2019$s.genes #grab S and G2M genes from Seurat
g2m.genes <- cc.genes.updated.2019$g2m.genes

IPAX_merged <- CellCycleScoring(IPAX_merged, # Calculate module scores for each gene set
                                s.features = s.genes, 
                                g2m.features = g2m.genes, 
                                set.ident = TRUE)

# visualize categorical phase assignment based on module scores (figure 8a)
CC_scores <- FetchData(IPAX_merged, vars = c("G2M.Score", "S.Score", "Phase"))

cell_cycle_colors <- c("#4062bb",  # G1 
                       "#7eb2dd",  # S 
                       "#f45b69",  # G2M
                       "#ebebeb")  # Other

ggplot(CC_scores, aes(x = G2M.Score, y = S.Score, color = factor(Phase, levels = c("G1", "S", "G2M"))))+
  scale_color_manual(values = cell_cycle_colors[1:3])+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  geom_point(size = 0.5)+
  geom_abline(slope = 1, linetype = "dashed")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = "G2M-Phase Module Score", y = "S-Phase Module Score", color = NULL)+
  theme(panel.border = element_rect(linewidth = 0.4),
        legend.text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85))

# visualize phase assignments in UMAP (figure 8b)
DimPlot(IPAX_merged, reduction = "RNA_UMAP", group.by = "Phase", cols = cell_cycle_colors[c(1,3,2)])+
  theme(legend.position = "",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(linewidth = 0.3))+
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)


# visualize distribution of categorical phases in MCF7 clusters (figure 8c)
Discrete_CC_clusters <- FetchData(IPAX_merged, vars = c("Cell_Line", "Phase", "Clusters")) %>% #pull data
  add_column(y = 1) %>% 
  filter(Cell_Line == "MCF7") %>% 
  dplyr::select(-Cell_Line) %>% 
  group_by(Clusters, Phase) %>%
  dplyr::count()

# plot data
ggplot(Discrete_CC_clusters, aes(x = Clusters, y = n, fill = factor(Phase, levels = c("G1", "S", "G2M"))))+
  geom_col(position = "fill", color = "white")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = cell_cycle_colors[-4])+
  labs(y = "Frequency", fill = "Phase")+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

# compute high confidence assignments
IPAX_merged$HC_Phase <- case_when(
  IPAX_merged$Clusters == "0" & IPAX_merged$Phase == "S" ~ "HC_S",
  IPAX_merged$Clusters == "2" & IPAX_merged$Phase == "G1" ~ "HC_G1",
  IPAX_merged$Clusters == "3" & IPAX_merged$Phase == "G2M" ~ "HC_G2M",
  TRUE ~ "Low Confidence"
)

# visualize high confidence assignments in UMAP (figure 8d)
DimPlot(IPAX_merged, reduction = "RNA_UMAP", group.by = "HC_Phase", cols = cell_cycle_colors[c(1,3,2,4)])+
  theme(legend.position = "",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(linewidth = 0.3))+
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)

## visualize marker gene expression (figure 8i, top)
VlnPlot(IPAX_merged, group.by = "HC_Phase", features = c("CDKN1A", "MCM4", "CDK1"), pt.size = 0)


## Identify differentially expressed genes across high confidence phase assignments (figure 8j)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(DOSE)

# Load Hallmarks
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")


# Prepare gene lists (identify DEGs)
HC_G1_markers <- FindMarkers(IPAX_merged, ident.1 = "HC_G1", ident.2 = c("HC_S", "HC_G2M"), group.by = "HC_Phase")
HC_S_markers <- FindMarkers(IPAX_merged, ident.1 = "HC_S", ident.2 = c("HC_G1", "HC_G2M"), group.by = "HC_Phase")
HC_G2M_markers <- FindMarkers(IPAX_merged, ident.1 = "HC_G2M", ident.2 = c("HC_G1", "HC_S"), group.by = "HC_Phase")

# Extract significant genes and arrange acording to log2FC
G1_genes <- HC_G1_markers %>% 
  rownames_to_column("GeneSymbol") %>% 
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(GeneSymbol, avg_log2FC) %>% 
  deframe()

S_genes <- HC_S_markers %>% 
  rownames_to_column("GeneSymbol") %>% 
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(GeneSymbol, avg_log2FC) %>% 
  deframe()

G2M_genes <- HC_G2M_markers %>% 
  rownames_to_column("GeneSymbol") %>% 
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(GeneSymbol, avg_log2FC) %>% 
  deframe()

# Perform gene set enrichment analysis using MSigDB's hallmark set
G1_gsea_result <- GSEA(G1_genes, 
                       TERM2GENE = hallmark_sets %>% dplyr::select(gs_name, gene_symbol),
                       pvalueCutoff = 0.05,
                       nPermSimple = 10000,
                       eps = 0)

S_gsea_result <- GSEA(S_genes, 
                      TERM2GENE = hallmark_sets %>% dplyr::select(gs_name, gene_symbol),
                      pvalueCutoff = 0.05,
                      nPermSimple = 100000,
                      eps = 0)

G2M_gsea_result <- GSEA(G2M_genes, 
                        TERM2GENE = hallmark_sets %>% dplyr::select(gs_name, gene_symbol),
                        pvalueCutoff = 0.05,
                        nPermSimple = 100000,
                        eps = 0)

# Merge results into single data frame
all_gsea_results <- list(G1 = G1_gsea_result@result, S = S_gsea_result@result, G2M = G2M_gsea_result@result)
all_pathways <- unique(unlist(lapply(all_gsea_results, function(res) res$ID)))

merged_gsea <- Reduce(function(df1, df2) full_join(df1, df2, by = "ID"), 
                      lapply(all_gsea_results, function(res) 
                        res %>% dplyr::select(ID, NES) %>% rename_with(~paste0(.x, "_NES"), -ID))) %>%
  mutate(across(contains("_NES"), ~replace_na(.x, 0)))
colnames(merged_gsea) <- c("ID", "G1", "S", "G2M")

# filter selected hallmarks
merged_gsea_long <- merged_gsea %>%
  filter(G1 + S + G2M > -1) %>% 
  column_to_rownames("ID") %>% 
  .[c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),]

# visualize results
pheatmap::pheatmap(merged_gsea_long,
                   color = colorRampPalette(c("#2D7AB5", "white", "#BB3037"))(200), 
                   angle_col = 0,legend_breaks = c(-2, 0, 2), cutree_rows = 3
            )

## continuous predictions: peco ## ----
library(peco)
library(SingleCellExperiment)

# convert seurat object to sce
sce_IPAX_merged <- as.SingleCellExperiment(IPAX_merged, assay = "RNA")

## Convert gene symbols to ENSEMBLIDs
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_symbols <- rownames(sce_IPAX_merged)  # Extract gene names

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # extract ensembl IDs

conversion_table <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = mart
)

# Create a named vector for mapping
symbol_to_ensembl <- setNames(conversion_table$ensembl_gene_id, conversion_table$hgnc_symbol)

# Replace row names with Ensembl IDs (keeping unmatched genes as is)
rownames(sce_IPAX_merged) <- ifelse(rownames(sce_IPAX_merged) %in% names(symbol_to_ensembl),
                                    symbol_to_ensembl[rownames(sce_IPAX_merged)],
                                    rownames(sce_IPAX_merged))  # Keeps original name if not found


## Load training data and cyclic genes ##
data("sce_top101genes")
data("training_human")

# Identify cell cycle genes from IPAX data
cc_pred_genes <- base::intersect(rownames(sce_top101genes), rownames(sce_IPAX_merged))

# Use only cell cycle genes to predict theta
sce_IPAX_merged <- sce_IPAX_merged[cc_pred_genes,]

# quantile normalization
sce_IPAX_merged <- data_transform_quantile(sce_IPAX_merged)

# Theta prediction
pred_IPAX_merged <- cycle_npreg_outsample(
  Y_test = assay(sce_IPAX_merged, "cpm_quantNormed"),
  sigma_est = training_human$sigma[rownames(sce_IPAX_merged),],
  funs_est = training_human$cellcycle_function[rownames(sce_IPAX_merged)],
  method.trend = "trendfilter",
  ncores = 1,
  get_trend_estimates = FALSE
)

# Sanity check (cyclic CDK1 expression)
plot(y=pred_IPAX_merged$Y["ENSG00000170312",],
     x=pred_IPAX_merged$cell_times_est, main = "CDK1",
     ylab = "quantile normalized expression")
points(y=training_human$cellcycle_function[["ENSG00000170312"]](seq(0,2*pi, length.out=100)),
       x=seq(0,2*pi, length.out=100), col = "blue", pch =16)



# Add cell cycle prediction to seurat object
IPAX_merged$CC_Prediction <- pred_IPAX_merged$cell_times_est[Cells(IPAX_merged)]

# visualize theta score (figure 8g)
FeaturePlot(IPAX_merged, features = "CC_Prediction", reduction = "RNA_UMAP")

# sanity checks (smoothed phase-marker expression, figure 8i bottom)
Peco_sanity_data <- FetchData(IPAX_merged, vars = c("CC_Prediction", "CDKN1A","CDK1","MCM4")) %>% 
  pivot_longer(-CC_Prediction, names_to = "Gene", values_to = "Expression")

ggplot(Peco_sanity_data, aes(x = CC_Prediction, y = Expression, color = Gene))+
  geom_vline(xintercept = c(2, 3.7, 5.8), linetype = "dashed", color = "grey", linewidth = 0.3)+
  facet_wrap(.~factor(Gene, levels = c("CDKN1A", "MCM4", "CDK1")), scale = "free")+
  geom_smooth()+
  scale_color_manual(values = cell_cycle_colors[c(3, 1, 2)])+
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", "π", "2π"))+
  labs(fill = "Phase", x = "Predicted Theta")+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

# Seurat comparison (figure 8h)
peco_seurat <- FetchData(IPAX_merged, vars = c("CC_Prediction", "HC_Phase")) %>% 
  filter(HC_Phase != "Low Confidence")

ggplot(peco_seurat, aes(x = CC_Prediction, fill = HC_Phase))+
  geom_vline(xintercept = c(2, 3.7, 5.8), linetype = "dashed", color = "grey", linewidth = 0.3)+
  facet_wrap(.~HC_Phase, ncol = 1)+
  geom_density()+
  scale_y_continuous(breaks = c(0, 0.8))+
  scale_fill_manual(values = cell_cycle_colors[c(1, 3, 2)], labels = c("G1", "G2M", "S"))+
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", "π", "2π"))+
  labs(fill = "Phase", x = "Predicted Theta")+
  theme(strip.text = element_blank(),
        panel.border = element_rect(linewidth = 0.6))

## slingshot trajectory ## ----
library(slingshot)

IPAX_seurat = subset(IPAX_merged, subset = Cell_Line == "MCF7") %>% JoinLayers() # extract MCF7 cells

#Convert dataset into a SingleCellExperiment Object
IPAX_seurat_sce = as.SingleCellExperiment(IPAX_seurat, assay = "RNA")
reducedDims(IPAX_seurat_sce)$UMAP <- IPAX_seurat@reductions$RNA_UMAP@cell.embeddings ## attach RNA_UMAP to sce data

# Run Slingshot analysis
IPAX_seurat_sce = slingshot(IPAX_seurat_sce, clusterLabels = "Phase", reducedDim = "RNA_UMAP", dist.method = 'mnn', approx_points = FALSE, start.clus = "G1", end.clus = "G2M", extend = "n")


## Visualise results of slingshotanalysis (figure 8e)
# Extract PCA/UMAP coordinates
embedding <- reducedDims(IPAX_seurat_sce)$RNA_UMAP 
UMAP_df <- as.data.frame(embedding)
UMAP_df$pseudotime <- slingPseudotime(IPAX_seurat_sce)[,1]  # Extract first lineage pseudotime

# Add cluster labels
UMAP_df$cluster <- factor(IPAX_seurat_sce$Phase)

# Extract trajectory
curve_data <- slingCurves(IPAX_seurat_sce)[[1]]

# Convert curve data to a data frame
trajectory_df <- data.frame(x = curve_data$s[curve_data$ord, 1], 
                            y = curve_data$s[curve_data$ord, 2])

ggplot(UMAP_df %>% filter(RNAUMAP_1 < 4), aes(x = RNAUMAP_1, y = RNAUMAP_2, color = cluster))+
  geom_point(size = 0.8, alpha = 0.7, show.legend = F)+
  scale_color_manual(values = cell_cycle_colors[c(1, 3, 2, 4)])+
  scale_x_continuous(limits = c(-11, 0))+
  geom_path(data = trajectory_df, aes(x = x, y = y), color = "#2A3D45", size = 1.2)+
  labs(color = "Phase", x = "RNA UMAP 1", y = "RNA UMAP 2")+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_blank(),
        axis.text = element_blank())


## add pseudotime to seurat object
pseudotime <- slingPseudotime(IPAX_seurat_sce)  # Pseudotime
IPAX_seurat$pseudotime <- pseudotime[, 1] 


pseudotime_sanity_check <- FetchData(IPAX_seurat, vars = c("HC_Phase", "pseudotime", "CDK1", "MCM4", "CDKN1A")) %>% 
  pivot_longer(CDK1:CDKN1A, names_to = "Gene", values_to = "Expression") %>% 
  mutate(HC_Phase = gsub("HC_", "", HC_Phase))

## marker-gene expression over pseudotime (figure 8i, middle)
ggplot(pseudotime_sanity_check, aes(x = pseudotime, y = Expression, color = Gene))+
  facet_wrap(.~factor(Gene, levels = c("CDKN1A", "MCM4", "CDK1")), scales = "free")+
  geom_smooth()+
  scale_color_manual(values = cell_cycle_colors[c(3, 1, 2)])+
  scale_x_continuous(breaks = seq(0, 12, 3))+
  labs(x = "Pseudo Time")+
  theme(legend.position = "",
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))

## distribution of categorical phases over pseudotime (figure 8g)
ggplot(pseudotime_sanity_check %>% filter(HC_Phase != "Low Confidence"), aes(x = pseudotime, fill = HC_Phase))+
  facet_wrap(.~factor(HC_Phase, levels = c("G1", "S", "G2M")), ncol = 1, strip.position = "right")+
  geom_density()+
  scale_fill_manual(values = cell_cycle_colors[c(1, 3, 2)])+
  labs(x = "Pseudo Time", y = "Density")+
  scale_y_continuous(breaks = c(0, 0.5))+
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12))+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))


## save final seurat object
saveRDS(IPAX_seurat, file = "files/IPAX_merged_CC.rds")