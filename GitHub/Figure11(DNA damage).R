## libraries ## ----
library(Seurat)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
  
## load data ## ----
IPAX_seurat <- readRDS("files/IPAX_merged_CC.rds")
DefaultAssay(IPAX_merged) <- "RNA"

## ER expression across pseudotime (figure 11a) ## ----
# extract data 
ER_modality_data <- FetchData(IPAX_seurat, vars = c("pseudotime", "rna_ESR1", "pa_ER", "pla_ER:ER")) %>% 
  pivot_longer(rna_ESR1:`pla_ER:ER`, names_to = "feature", values_to = "count") 

ggplot(ER_modality_data, aes(x = pseudotime, y = count))+
  facet_wrap(.~feature, scales = "free")+
  geom_smooth()

## ER activity phase-specific hallmarks (figure 11b) ## ----
## S
PIPE_MCF7_S <- subset(IPAX_seurat, subset = Phase == "S") # select S-phase cells
PIPE_MCF7_S$P3001_class <- ifelse(FetchData(PIPE_MCF7_S, vars = "ER:ER") > 2, "high", "low") # bin ER:ER intraction

# Identify DEGs across the bins
Diff <- FindMarkers(
  PIPE_MCF7_S,
  ident.1 = "high",
  ident.2 = "low",
  assay = "RNA",
  test.use = "wilcox"
)

# import hallmark set
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)

# pull significant genes
significant_genes <- Diff[Diff[["p_val"]] < 0.05 & Diff$avg_log2FC > 0,]

# perform gene set enrichment analysis
clust1 <- enricher(rownames(significant_genes), TERM2GENE=m_t2g, pvalueCutoff =1, maxGSSize = 10000)
clust1 <- as.data.frame(clust1@result)

# rearrange data frame
clust1 <- clust1[,c("ID","pvalue")]
colnames(clust1) <- c("pathway","cluster1_S")

# grab hallmarks of interest
new_data_S <- clust1[grepl("GOBP_DNA_REPAIR|GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE|GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_DOUBLE_STRAND_BREAK_REPAIR|GOBP_BASE_EXCISION_REPAIR|GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING", clust1$pathway), ]

## repeat procedure for G1
PIPE_MCF7_G1 <- subset(IPAX_seurat, subset = Phase == "G1")
PIPE_MCF7_G1$P3001_class <- ifelse(FetchData(PIPE_MCF7_G1, vars = "ER:ER") > 2, "high", "low")

Diff <- FindMarkers(
  PIPE_MCF7_G1,
  ident.1 = "high",
  ident.2 = "low",
  assay = "RNA", 
  test.use = "wilcox"
)

significant_genes <- Diff[Diff[["p_val"]] < 0.05 & Diff$avg_log2FC > 0,]
clust1 <- enricher(rownames(significant_genes), TERM2GENE=m_t2g, pvalueCutoff =1, maxGSSize = 10000)
clust1 <- as.data.frame(clust1@result)

clust1 <- clust1[,c("ID","pvalue")]
colnames(clust1) <- c("pathway","cluster1_G1")


new_data_G1 <- clust1[grepl("GOBP_DNA_REPAIR|GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE|GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_DOUBLE_STRAND_BREAK_REPAIR|GOBP_BASE_EXCISION_REPAIR|GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING", clust1$pathway), ]

# repeat once more for G2M
PIPE_MCF7_G2M <- subset(IPAX_seurat, subset = Phase == "G2M")
PIPE_MCF7_G2M$P3001_class <- ifelse(FetchData(PIPE_MCF7_G2M, vars = "ER:ER") > 2, "high", "low")


Diff <- FindMarkers(
  PIPE_MCF7_G2M,
  ident.1 = "high",
  ident.2 = "low",
  assay = "RNA",
  test.use = "wilcox"
)

significant_genes <- Diff[Diff[["p_val"]] < 0.05 & Diff$avg_log2FC > 0,]
clust1 <- enricher(rownames(significant_genes), TERM2GENE=m_t2g, pvalueCutoff =1, maxGSSize = 10000)
clust1 <- as.data.frame(clust1@result)

clust1 <- clust1[,c("ID","pvalue")]
colnames(clust1) <- c("pathway","cluster1_G2M")


new_data_G2M <- clust1[grepl("GOBP_DNA_REPAIR|GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE|GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS|GOBP_DOUBLE_STRAND_BREAK_REPAIR|GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING", clust1$pathway), ]


# merge results from all 3 phases
new_data <- merge(new_data_S, new_data_G1, by = "pathway") %>% merge(new_data_G2M, by = "pathway")
rownames(new_data) <- gsub(x = new_data$pathway, pattern = "GOBP_", replacement = "")
new_data <- new_data[,c("cluster1_G1","cluster1_S", "cluster1_G2M")]
new_data <- as.matrix(new_data[,c("cluster1_G1", "cluster1_S", "cluster1_G2M")])

# visualize results
library(pheatmap)
pheatmap::pheatmap(-log10(new_data), scale = "none", color = colorRampPalette(c("white","#2D7AB5"))(200), name = "pvalue", cluster_cols = F)

## DRIP-seq tag distribution (figure 11c) ## ----
# load data
DRIP_peaks <- non_overlap_drip <- read.delim("files/DRIP/E2_peaks_annotated.txt") %>%
  dplyr::select("Vehicle" = Vehicle_TagDir.Tag.Count.in.given.bp..11021921.0.Total..normalization.factor...0.91..effective.total...10000000.,
                "E2" = E2_TagDir.Tag.Count.in.given.bp..12317594.0.Total..normalization.factor...0.81..effective.total...10000000.,
                everything()) %>% 
  mutate(Annotation = gsub(" \\(.*", "", Annotation), # remove details from annotation
         Annotation = ifelse(is.na(Annotation), "Other", Annotation), # remove resulting NAs
         Vehicle = ifelse(is.na(Vehicle), 0, Vehicle),
         E2 = ifelse(is.na(E2), 0, E2))

# visualize distribution
ggplot(DRIP_peaks %>% pivot_longer(Vehicle:E2, names_to = "Condition", values_to = "Count"), aes(x = factor(Condition, levels = c("Vehicle", "E2")), y = Count, fill = factor(Annotation, levels = c("promoter-TSS", "5' UTR", "exon", "intron", "3' UTR", "TTS", "non-coding", "Intergenic", "Other"))))+
  geom_col(position = "fill")+
  labs(x = NULL, y = "Tag Count", fill = NULL)

## DRIP-seq average plots (figure 11d) ## ----
# load data
overlap_drip <- read.delim("files/DRIP/overlap_DRIP.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.E2_ER_overlaps.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Vehicle" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  pivot_longer(2:3, names_to = "Condition", values_to = "Count") %>% 
  add_column(Site = "ER (n = 107)")

non_overlap_drip <- read.delim("files/DRIP/non_overlap_DRIP.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.E2_exclusive.bed.hg38..d.Vehicle_TagDir.E2_TagDir..hist.50..size.5000.,
                "Vehicle" = Vehicle_TagDir.Coverage,
                "E2" = E2_TagDir.Coverage) %>% 
  pivot_longer(2:3, names_to = "Condition", values_to = "Count") %>% 
  add_column(Site = "Non-ER (n = 15597)")

# merge data sets
DRIP_hist_merged <- rbind(overlap_drip, non_overlap_drip)

# visualize
ggplot(DRIP_hist_merged, aes(x = DFC, y = Count, color = Condition))+
  facet_wrap(.~Site)+
  geom_line(show.legend = F)+
  scale_color_manual(values = c(rep_colors[2], "grey"))+
  labs(x = NULL, y = "Tag Count")+
  scale_x_continuous(breaks = c(-2000, 0, 2000), labels = c("-2kb", "Center", "2kb"))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))