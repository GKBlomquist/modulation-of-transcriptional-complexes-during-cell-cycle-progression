library(survival)
library(survminer)

# The Metabrix gene expression matrix can be downloaded from cBioportal webpage. 
mb_exp2 <- read.delim("path/to/expression_matrix", sep = " ")

#Import metastasis data from Rueda et al 2019. This is clinical data for the metabric cohort downloaded as suppl table S6 from Rueda et al, Nature 2019. The death column is the same as from cbioportal as expected. Info on each column can be found in 41586_2019_1007_MOESM10_ESM_Metabric_variableDescription.xls.
relapse <- read.delim("path/to/metastasis_data",h=T) 

# Cluster gene vectors (i.e., Cluster1_genes, Cluster2_genes, etc.) are established in "Figure12(ER targets).R" please run that script prior to this one
# Alternatively you can import cluster genes from my analysis using the following line of code for each cluster: Cluster1_genes <- readLines("files/Clust1_genes.txt")

## Cluster 1 
# Subset genes from the metabric dataset and add clinical information
gene_set_C1 <- mb_exp2 %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% Cluster1_genes) %>% 
  pull(Gene)

km_C1 <- as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% gene_set_C1,])) %>% #Get gene expression for selected genes
  rownames_to_column("METABRIC.ID") %>% 
  left_join(relapse) %>% # Join to relapse data
  mutate(Mean = rowMeans(across(all_of(gene_set_C1)), na.rm = TRUE)) %>%  # calculate mean expression of selected gene set 
  filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-") %>% # Subset ER+ patients
  mutate(type = cut(Mean, breaks = 3, labels = c("Low", "Medium", "High"))) %>%  # Divide patients into "Low", "Medium", and "High" groups in terms of mean gene expression from provided gene list
  filter(type != "Medium") # Remove "Medium" patients


#Prepare data for KM plot. Analyse distant metastasis-free survival.
km_fit <- survfit(Surv(TDR, DR) ~ type, data=km_C1) 
surv_pvalue(km_fit)
ggsurvplot(km_fit, pval = T, pval.method = F, conf.int = F,
           pval.coord = c(0, 0.03),
           title = "Cluster 4", size = 1.5,
           xlab = "Years",
           ylab = "Metastasis-free Survival",
           legend = "none", 
           censor = F,
           ggtheme = theme(
             axis.line = element_line(linewidth = 0.4),
             panel.border = element_rect(linewidth = 0.4), legend.position = ""), 
           palette = c("grey", rep_colors[2]), size = 2)$plot +
  scale_x_continuous(breaks = c(0, 365.25*5, 365.25*10 ,365.25*15, 20*365.25, 25*365.25), labels = seq(0, 25, 5))+
  scale_y_continuous(labels = scales::percent)

## Cluster 2
# Subset genes from the metabric dataset and add clinical information
gene_set_C2 <- mb_exp2 %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% Cluster2_genes) %>% 
  pull(Gene)

km_C2 <- as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% gene_set_C2,])) %>% #Get gene expression for selected genes
  rownames_to_column("METABRIC.ID") %>% 
  left_join(relapse) %>% # Join to relapse data
  mutate(Mean = rowMeans(across(all_of(gene_set_C2)), na.rm = TRUE)) %>%  # calculate mean expression of selected gene set 
  filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-") %>% # Subset ER+ patients
  mutate(type = cut(Mean, breaks = 3, labels = c("Low", "Medium", "High"))) %>%  # Divide patients into "Low", "Medium", and "High" groups in terms of mean gene expression from provided gene list
  filter(type != "Medium") # Remove "Medium" patients


#Prepare data for KM plot. Analyse distant metastasis-free survival.
km_fit <- survfit(Surv(TDR, DR) ~ type, data=km_C2) 
surv_pvalue(km_fit)
ggsurvplot(km_fit, pval = T, pval.method = F, conf.int = F,
           pval.coord = c(0, 0.03),
           title = "Cluster 4", size = 1.5,
           xlab = "Years",
           ylab = "Metastasis-free Survival",
           legend = "none", 
           censor = F,
           ggtheme = theme(
             axis.line = element_line(linewidth = 0.4),
             panel.border = element_rect(linewidth = 0.4), legend.position = ""), 
           palette = c("grey", rep_colors[2]), size = 2)$plot +
  scale_x_continuous(breaks = c(0, 365.25*5, 365.25*10 ,365.25*15, 20*365.25, 25*365.25), labels = seq(0, 25, 5))+
  scale_y_continuous(labels = scales::percent)

## Cluster 3
# Subset genes from the metabric dataset and add clinical information
gene_set_C3 <- mb_exp2 %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% Cluster3_genes) %>% 
  pull(Gene)

km_C3 <- as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% gene_set_C3,])) %>% #Get gene expression for selected genes
  rownames_to_column("METABRIC.ID") %>% 
  left_join(relapse) %>% # Join to relapse data
  mutate(Mean = rowMeans(across(all_of(gene_set_C3)), na.rm = TRUE)) %>%  # calculate mean expression of selected gene set 
  filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-") %>% # Subset ER+ patients
  mutate(type = cut(Mean, breaks = 3, labels = c("Low", "Medium", "High"))) %>%  # Divide patients into "Low", "Medium", and "High" groups in terms of mean gene expression from provided gene list
  filter(type != "Medium") # Remove "Medium" patients


#Prepare data for KM plot. Analyse distant metastasis-free survival.
km_fit <- survfit(Surv(TDR, DR) ~ type, data=km_C3) 
surv_pvalue(km_fit)
ggsurvplot(km_fit, pval = T, pval.method = F, conf.int = F,
           pval.coord = c(0, 0.03),
           title = "Cluster 4", size = 1.5,
           xlab = "Years",
           ylab = "Metastasis-free Survival",
           legend = "none", 
           censor = F,
           ggtheme = theme(
             axis.line = element_line(linewidth = 0.4),
             panel.border = element_rect(linewidth = 0.4), legend.position = ""), 
           palette = c("grey", rep_colors[2]), size = 2)$plot +
  scale_x_continuous(breaks = c(0, 365.25*5, 365.25*10 ,365.25*15, 20*365.25, 25*365.25), labels = seq(0, 25, 5))+
  scale_y_continuous(labels = scales::percent)

## Cluster 4
# Subset genes from the metabric dataset and add clinical information
gene_set_C4 <- mb_exp2 %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% Cluster4_genes) %>% 
  pull(Gene)

km_C2 <- as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% gene_set_C4,])) %>% #Get gene expression for selected genes
  rownames_to_column("METABRIC.ID") %>% 
  left_join(relapse) %>% # Join to relapse data
  mutate(Mean = rowMeans(across(all_of(gene_set_C2)), na.rm = TRUE)) %>%  # calculate mean expression of selected gene set 
  filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-") %>% # Subset ER+ patients
  mutate(type = cut(Mean, breaks = 3, labels = c("Low", "Medium", "High"))) %>%  # Divide patients into "Low", "Medium", and "High" groups in terms of mean gene expression from provided gene list
  filter(type != "Medium") # Remove "Medium" patients


#Prepare data for KM plot. Analyse distant metastasis-free survival.
km_fit <- survfit(Surv(TDR, DR) ~ type, data=km_C4) 
surv_pvalue(km_fit)
ggsurvplot(km_fit, pval = T, pval.method = F, conf.int = F,
           pval.coord = c(0, 0.03),
           title = "Cluster 4", size = 1.5,
           xlab = "Years",
           ylab = "Metastasis-free Survival",
           legend = "none", 
           censor = F,
           ggtheme = theme(
             axis.line = element_line(linewidth = 0.4),
             panel.border = element_rect(linewidth = 0.4), legend.position = ""), 
           palette = c("grey", rep_colors[2]), size = 2)$plot +
  scale_x_continuous(breaks = c(0, 365.25*5, 365.25*10 ,365.25*15, 20*365.25, 25*365.25), labels = seq(0, 25, 5))+
  scale_y_continuous(labels = scales::percent)