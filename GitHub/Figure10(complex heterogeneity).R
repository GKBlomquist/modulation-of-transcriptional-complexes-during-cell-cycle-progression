## libraries ## ----
library(tidyverse)
library(patchwork)
greys <- gray.colors(5, start = 0.9, end = 0.6)

## load GAM splines ## ----
complex_splines <- read.csv("files/complex_splines.csv") %>% 
  column_to_rownames(var = "X")

# create additional data fram with scaled splines
complex_splines_scaled <- read.csv("files/complex_splines.csv") %>% 
  column_to_rownames(var = "X")
  t(scale(t(complex_splines))) %>% 
  data.frame() %>% 
  rownames_to_column(var = "Complex")
colnames(complex_splines_scaled) <- c("Complex", seq(1, 12.5, length.out = 100))

complex_splines_scaled <- complex_splines_scaled %>% 
  pivot_longer(2:101, names_to = "pseudotime", values_to = "Smooth") %>% 
  mutate(pseudotime = as.numeric(pseudotime))

# reformat complex_splines
complex_splines <- complex_splines %>% 
  pivot_longer(2:101, names_to = "pseudotime", values_to = "Smooth") %>% 
  mutate(pseudotime = gsub("lineage_", "", pseudotime),
         pseudotime = as.numeric(pseudotime))

## ER complex spline plots (figure 10a) ----
# FOXA1:ER
ggplot(complex_splines %>% filter(Complex %in% c("FOXA1:ER", "FLAG:FOXA1", "FLAG:ER")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.7, 1, 1.4))+
  scale_y_break(breaks = c(0.12, 0.7), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-FOXA1")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# ER:ER
ggplot(complex_splines %>% filter(Complex %in% c("ER:ER", "FLAG:ER")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("darkgrey", "#008A5A")))+
  scale_y_continuous(breaks = c(0, 0.12, 3, 6, 9))+
  scale_y_break(breaks = c(0.12, 3), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-ER")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# ER:MED1
ggplot(complex_splines %>% filter(Complex %in% c("ER:MED1-p", "FLAG:ER", "FLAG:MED1-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("darkgrey", "lightgrey","#008A5A")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.2, 0.3, 0.4))+
  scale_y_break(breaks = c(0.12, 0.2), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-MED1")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# ER:BRG1
ggplot(complex_splines %>% filter(Complex %in% c("ER:BRG1", "FLAG:ER", "FLAG:BRG1")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("darkgrey", "lightgrey","#008A5A")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.2, 0.3, 0.4))+
  scale_y_break(breaks = c(0.12, 0.2), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-BRG1")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# ER:MYC
ggplot(complex_splines %>% filter(Complex %in% c("MYC-p:ER", "FLAG:ER", "FLAG:MYC-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.25, 0.37, 0.5))+
  scale_y_break(breaks = c(0.12, 0.25), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-MYC")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# ER:P300
ggplot(complex_splines %>% filter(Complex %in% c("P300:ER", "FLAG:ER", "FLAG:P300")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.5, 1, 1.5))+
  scale_y_break(breaks = c(0.12, 0.5), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-P300")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

## Scaled ER complexes plots (figure 10b)
ggplot(complex_splines_scaled %>% filter(Complex %in% c("ER:ER", "ER:MED1-p", "P300:ER", "ER:BRG1","MYC-p:ER")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line(linewidth = 0.3)+
  scale_color_manual(values = greys)+
  geom_line(data = complex_splines_scaled %>% filter(Complex %in% c("FOXA1:ER")), color = "#DC4F4F")+
  labs(x = "Pseudotime", y = "Z score", title = "ER-FOXA1")+
  theme(legend.position = "",
        axis.line = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0))+
  
  ggplot(complex_splines_scaled %>% filter(Complex %in% c("ER:ER", "ER:MED1-p", "P300:ER", "FOXA1:ER","MYC-p:ER")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line(linewidth = 0.3)+
  scale_color_manual(values = greys)+
  geom_line(data = complex_splines_scaled %>% filter(Complex %in% c("ER:BRG1")), color = "#DC4F4F")+
  labs(x = "Pseudotime", y = "Z score", title = "ER-BRG1")+
  theme(legend.position = "",
        axis.line = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0))+
  
  ggplot(complex_splines_scaled %>% filter(Complex %in% c("ER:ER", "ER:MED1-p", "FOXA1:ER", "ER:BRG1","MYC-p:ER")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line(linewidth = 0.3)+
  scale_color_manual(values = greys)+
  geom_line(data = complex_splines_scaled %>% filter(Complex %in% c("P300:ER")), color = "#DC4F4F")+
  labs(x = "Pseudotime", y = "Z score", title = "ER-P300")+
  theme(legend.position = "",
        axis.line = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0))

### MYC complexes (figure 10c) ----
# ER:MYC
ggplot(complex_splines %>% filter(Complex %in% c("MYC-p:ER", "FLAG:ER", "FLAG:MYC-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.12, 0.25, 0.37, 0.5))+
  scale_y_break(breaks = c(0.12, 0.25), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "ER-MYC")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.35),
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# MYC:BRD4
ggplot(complex_splines %>% filter(Complex %in% c("MYC-p:BRD4-p", "FLAG:BRD4-p", "FLAG:MYC-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))+
  scale_y_break(breaks = c(0.025, 0.05), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "MYC-BRD4")+
  theme(legend.position = "inside",
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# MYC:MED1
ggplot(complex_splines %>% filter(Complex %in% c("MYC-p:MED1-p", "FLAG:MED1-p", "FLAG:MYC-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = rev(c("#008A5A", "darkgrey", "lightgrey")))+
  scale_y_continuous(breaks = c(0, 0.03, 0.1, 0.12, 0.18, 0.24))+
  scale_y_break(breaks = c(0.03, 0.12), scales = "fixed", space = 0.3)+
  labs(x = "Pseudotime", y = "Normalized Count", color = NULL, title = "MYC-MED1")+
  theme(legend.position = "inside",
        panel.border = element_rect(linewidth = 0),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))

# summary
ggplot(complex_splines_scaled %>% filter(Complex %in% c("MYC-p:ER", "MYC-p:BRD4-p")), aes(x = pseudotime, y = Smooth, color = Complex))+
  geom_line()+
  scale_color_manual(values = c("#DC4F4F", "#6EA5CC"))+
  labs(x = "Pseudotime", y = "Z score")+
  theme(legend.position = "",
        axis.line = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0))

# protein expression (figure 10d)
DefaultAssay(IPAX_seurat) <- "PA"
protein_expression_data <- FetchData(IPAX_seurat, vars = c("pseudotime", "FOXA1", "ER", "MED1-p", "MYC-p", "BRG1", "P300", "BRD4-p", "FLAG", "IgG")) %>% 
  pivot_longer(FOXA1:IgG, names_to = "Protein", values_to = "Scaled Count")

Proteins_order <- c("FOXA1", "ER", "MED1-p", "MYC-p", "BRG1", "P300", "BRD4-p", "FLAG", "IgG")

ggplot(protein_expression_data, aes(x = pseudotime, y = `Scaled Count`))+
  facet_wrap(.~factor(Protein, levels = Proteins_order))+
  scale_x_continuous(breaks = seq(0, 12, 3))+
  geom_smooth(color = rep_colors[2])+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.2))

## MYC ChIP-seq average plot (figure 10e) ----
# load histogram data computed by HOMER (see files directory)
promoter_hist <- read.delim("files/ChIP/figure10/Fulvestrant_hist_MYC_promoters.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.MYC_promoters.bed.hg38..d.MYC_FULV_TagDir.MYC_VEH_TagDir..hist.50..size.5000.,
                "FULV" = MYC_FULV_TagDir.Coverage,
                "VEH" = MYC_VEH_TagDir.Coverage) %>% 
  pivot_longer(-DFC, names_to = "Condition", values_to = "TagCount")

enhancer_hist <- read.delim("files/ChIP/figure10/Fulvestrant_hist_MYC_ER_enhancers.txt") %>% 
  dplyr::select("DFC" = Distance.from.Center..cmd.annotatePeaks.pl.ER_enhancers.bed.hg38..d.MYC_FULV_TagDir.MYC_VEH_TagDir..hist.50..size.5000.,
                "FULV" = MYC_FULV_TagDir.Coverage,
                "VEH" = MYC_VEH_TagDir.Coverage) %>% 
  pivot_longer(-DFC, names_to = "Condition", values_to = "TagCount")

# plot histograms
ggplot(promoter_hist, aes(x = DFC, y = TagCount, color = Condition))+
  geom_line(linewidth = 1, show.legend = F)+
  scale_color_manual(values = c("#008A5A", "grey"))+
  labs(x = NULL, y = "MYC Tag Count", title = "MYC Promoter")+
  scale_x_continuous(labels = c("-2kb", "Center", "2kb"), breaks = c(-2000, 0, 2000))+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))

ggplot(enhancer_hist, aes(x = DFC, y = TagCount, color = Condition))+
  geom_line(linewidth = 1, show.legend = F)+
  scale_color_manual(values = c("#008A5A", "grey"))+
  labs(x = NULL, y = "MYC Tag Count", title = "ER Enhancer")+
  scale_x_continuous(labels = c("-2kb", "Center", "2kb"), breaks = c(-2000, 0, 2000))+
  theme(panel.border = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth = 0.4))
