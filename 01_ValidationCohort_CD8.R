## Author: Nikolaj Pagh Kristensen
## Date: 23-03-2026
## Title: Reanalyzing Kavaka et al 2024 "Twin study identifies early immunological and 
##   metabolic dysregulation of CD8 T cells in multiple sclerosis". 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(readxl)
library(RColorBrewer)
library(scCustomize)
library(cowplot)

#global theme
theme_set(theme_pubr(base_size = 7))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ---------------------------------------------------------------¨
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The preprocessed seurat data object was downloaded from GSM8492283
obj <- readRDS("04_data_objects/03_KavakaSciImmunol_TwinStudy/GSM8492283_Validation_cohort_pbmc.rds/GSM8492283_Validation_cohort_pbmc.rds")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Annotate the metadata with a CSF-biased score  ------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Make sure we are using the RNA assay
DefaultAssay(obj) <- "SymphonyQuery"

# get all object metadata
obj_meta <- obj[[]]

# get list of metadata columns
colnames(obj_meta)
obj_meta %>% select(diagnosis, diagnosis_simp, sample, samplenumb) %>% unique


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%  Calculate a score based on the genes you found previously  %% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEG_biased = read_csv("06_csv_outputs/DEG/CD8_CSF_biased_filtered_DEGs.csv") %>% select(names) %>% unique %>% pull
obj <- AddModuleScore(
  object = obj,
  features = list(DEG_biased),
  slot = "data", 
  ctrl = 50,                 #Number of control genes set to 50 to mirror Scanpy's implementation of the same function
  name = 'csf_features'
)

PB_geneset = read_csv("06_csv_outputs/DEG/CD8_PB_biased_filtered_DEGs.csv") %>% select(names) %>% unique %>% pull
obj <- AddModuleScore(
  object = obj,
  features = list(PB_geneset),
  slot = "data", 
  ctrl = 50,
  name = 'pb_features'
)

#Update the metadata shorthand
obj_meta <- obj[[]]
obj_meta

# Show the biased gene set calculations on a umap
DimPlot(obj, reduction="umap", pt.size = 1)
# FeaturePlot(obj, features = "csf_features1", cols = c("grey95", "grey50", "#700320"), order=TRUE)

## FeaturePlot improved by scCustimize
# Set color palette
base_cols <- brewer.pal(11, "RdBu")
pal <- colorRampPalette(rev(base_cols))(10)

CSF_score_overlay = FeaturePlot_scCustom(seurat_object = obj, features = "csf_features1", colors_use = pal)
PB_score_overlay = FeaturePlot_scCustom(seurat_object = obj, features = "pb_features1", colors_use = pal)


# Boxplot - mean intraclonal % of cells in CSF ----------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%  Calculate the total expansion and compartment-wise freq  %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CSF TCR frequency is already in there, see column "CSF_tcr_frequency"
# See prepro steps here: https://github.com/beltranLab/peripheral-cd8-mstwinstudy/blob/main/Validation_Cohort/Mapping_symphony_PBMC_Github.ipynb
obj_meta = obj_meta %>% as_tibble() %>% 
  mutate(csf_tcr_frequency = ifelse(csf_tcr_frequency == "FALSE", "0", csf_tcr_frequency)) %>% 
  mutate(csf_tcr_frequency = as.numeric(csf_tcr_frequency)) %>% 
  mutate(total_tcr_freq = csf_tcr_frequency + TCR_frequency_corrected) %>%  #I use TCR_frequency_corrected, which the authors also use in their code. 
  mutate(CSF = 100*csf_tcr_frequency/total_tcr_freq) %>% 
  mutate(PB = 100*TCR_frequency_corrected/total_tcr_freq)


## Define 5 data frames
### df0    dataframe without any expansion threshold. What is the average clonal %CSF frequency
df0 = obj_meta %>% 
  filter(total_tcr_freq > 0) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  group_by(samplenumb, diagnosis) %>% 
  summarise(mean_csf_freq = mean(CSF)) %>% 
  ungroup %>% 
  mutate(expan_thres = "> 0") 

### df5    dataframe with 5 > expansion threshold
df5 = obj_meta %>% 
  filter(total_tcr_freq > 5) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  group_by(samplenumb, diagnosis) %>% 
  summarise(mean_csf_freq = mean(CSF)) %>% 
  ungroup %>% 
  mutate(expan_thres = "> 5") 

### df10    dataframe with 10 > expansion threshold
df10 = obj_meta %>% 
  filter(total_tcr_freq > 10) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  group_by(samplenumb, diagnosis) %>% 
  summarise(mean_csf_freq = mean(CSF)) %>% 
  ungroup %>% 
  mutate(expan_thres = "> 10") 

### df15    dataframe with 15 > expansion threshold
df15 = obj_meta %>% 
  filter(total_tcr_freq > 15) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  group_by(samplenumb, diagnosis) %>% 
  summarise(mean_csf_freq = mean(CSF)) %>% 
  ungroup %>% 
  mutate(expan_thres = "> 15") 

### df15    dataframe with 15 > expansion threshold
df20 = obj_meta %>% 
  filter(total_tcr_freq > 20)%>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  group_by(samplenumb, diagnosis) %>% 
  summarise(mean_csf_freq = mean(CSF)) %>% 
  ungroup %>% 
  mutate(expan_thres = "> 20") 

## Bind the dataframes together
df_total = rbind(df0, df5, df10, df15, df20)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%   Plot the mean % CSF per clone over expan. thresholds    %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##   notes: I first make the expan_thres into factors to easier order the x-axis
df_total$expan_thres = factor(df_total$expan_thres, levels = c("> 0", "> 5", "> 10", "> 15", "> 20"))

expan_csf_freq_plot = ggplot(df_total, aes(x = expan_thres, y = mean_csf_freq+0.005, fill = diagnosis))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(dodge.width = 0.8, width = 0.03, size = 1.5)+
  geom_hline(yintercept = 0.005, linetype="dashed")+
  scale_y_log10(guide = "axis_logticks")+
  ylab("Mean % of cells in CSF per clone")+
  theme_pubr()
expan_csf_freq_plot

## Statistics. y as as function of x (y ~ x)
df_total %>% group_by(expan_thres) %>% 
  wilcox_test(mean_csf_freq ~ diagnosis) %>% 
  arrange(p)

## Simplified comparison - control vs ms
stats_df = df_total %>% 
  mutate(catdiag = ifelse(diagnosis == "IIH", "Ctrl", "MS"))
stats_df %>% group_by(expan_thres) %>% 
  wilcox_test(mean_csf_freq ~ catdiag) %>% 
  arrange(p)


# %%%%%%%%% 
# Export
# %%%%%%%%%
# Arrange with cowplot
plt_fin = plot_grid(expan_csf_freq_plot, ncol = 2, nrow = 4)

ggsave(filename = "05_plots/06_ValidationCohort/2026_03 - ExpansioDegree per clone in CSF.pdf",
       plot = plt_fin, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")


# Boxplot - % cells in CSF ----------------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%  Compartment-wise % of all cells  %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CSF TCR frequency is already in there, see column "CSF_tcr_frequency"
# See prepro steps here: https://github.com/beltranLab/peripheral-cd8-mstwinstudy/blob/main/Validation_Cohort/Mapping_symphony_PBMC_Github.ipynb
obj_meta = obj_meta %>% as_tibble() %>% 
  mutate(csf_tcr_frequency = ifelse(csf_tcr_frequency == "FALSE", "0", csf_tcr_frequency)) %>% 
  mutate(csf_tcr_frequency = as.numeric(csf_tcr_frequency)) %>% 
  mutate(total_tcr_freq = csf_tcr_frequency + TCR_frequency_corrected) %>%  #I use TCR_frequency_corrected, which the authors also use in their code. 
  mutate(CSF = 100*csf_tcr_frequency/total_tcr_freq) %>% 
  mutate(PB = 100*TCR_frequency_corrected/total_tcr_freq)

## Define 5 data frames
### df0    dataframe without any expansion threshold. What is the average clonal %CSF frequency
df0 = obj_meta %>% 
  filter(total_tcr_freq > 0) %>% 
  select(diagnosis, samplenumb, Sample_Clono, csf_tcr_frequency, TCR_frequency_corrected, total_tcr_freq) %>% 
  unique %>% 
  group_by(diagnosis, samplenumb) %>% 
  summarise(nTotal = sum(total_tcr_freq), 
            nCSF = sum(csf_tcr_frequency), 
            nPB = sum(TCR_frequency_corrected)) %>% 
  ungroup() %>% 
  mutate(CSF = 100*nCSF/nTotal,
         PB = 100*nPB/nTotal) %>% 
  mutate(expan_thres = "> 0") 

### df5    dataframe with 5 > expansion threshold
df5 = obj_meta %>% 
  filter(total_tcr_freq > 5) %>% 
  select(diagnosis, samplenumb, Sample_Clono, csf_tcr_frequency, TCR_frequency_corrected, total_tcr_freq) %>% 
  unique %>% 
  group_by(diagnosis, samplenumb) %>% 
  summarise(nTotal = sum(total_tcr_freq), 
            nCSF = sum(csf_tcr_frequency), 
            nPB = sum(TCR_frequency_corrected)) %>% 
  ungroup() %>% 
  mutate(CSF = 100*nCSF/nTotal,
         PB = 100*nPB/nTotal) %>% 
  mutate(expan_thres = "> 5") 

### df10    dataframe with 10 > expansion threshold
df10 = obj_meta %>% 
  filter(total_tcr_freq > 10) %>% 
  select(diagnosis, samplenumb, Sample_Clono, csf_tcr_frequency, TCR_frequency_corrected, total_tcr_freq) %>% 
  unique %>% 
  group_by(diagnosis, samplenumb) %>% 
  summarise(nTotal = sum(total_tcr_freq), 
            nCSF = sum(csf_tcr_frequency), 
            nPB = sum(TCR_frequency_corrected)) %>% 
  ungroup() %>% 
  mutate(CSF = 100*nCSF/nTotal,
         PB = 100*nPB/nTotal) %>% 
  mutate(expan_thres = "> 10") 

### df15    dataframe with 15 > expansion threshold
df15 = obj_meta %>% 
  filter(total_tcr_freq > 15) %>% 
  select(diagnosis, samplenumb, Sample_Clono, csf_tcr_frequency, TCR_frequency_corrected, total_tcr_freq) %>% 
  unique %>% 
  group_by(diagnosis, samplenumb) %>% 
  summarise(nTotal = sum(total_tcr_freq), 
            nCSF = sum(csf_tcr_frequency), 
            nPB = sum(TCR_frequency_corrected)) %>% 
  ungroup() %>% 
  mutate(CSF = 100*nCSF/nTotal,
         PB = 100*nPB/nTotal) %>% 
  mutate(expan_thres = "> 15") 

### df15    dataframe with 15 > expansion threshold
df20 = obj_meta %>% 
  filter(total_tcr_freq > 20) %>% 
  select(diagnosis, samplenumb, Sample_Clono, csf_tcr_frequency, TCR_frequency_corrected, total_tcr_freq) %>% 
  unique %>% 
  group_by(diagnosis, samplenumb) %>% 
  summarise(nTotal = sum(total_tcr_freq), 
            nCSF = sum(csf_tcr_frequency), 
            nPB = sum(TCR_frequency_corrected)) %>% 
  ungroup() %>% 
  mutate(CSF = 100*nCSF/nTotal,
         PB = 100*nPB/nTotal) %>% 
  mutate(expan_thres = "> 20") 

## Bind the dataframes together
df_total = rbind(df0, df5, df10, df15, df20)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%   Plot the mean % CSF per clone over expan. thresholds    %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##   notes: I first make the expan_thres into factors to easier order the x-axis
df_total$expan_thres = factor(df_total$expan_thres, levels = c("> 0", "> 5", "> 10", "> 15", "> 20"))

expan_csf_freq_plot = ggplot(df_total, aes(x = expan_thres, y = CSF+0.005, fill = diagnosis))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(dodge.width = 0.8, width = 0.03, size = 1.5)+
  scale_y_log10(guide = "axis_logticks")+
  ylab("% of cells in CSF")+
  theme_pubr()
expan_csf_freq_plot

## Statistics. y as as function of x (y ~ x)
### Note: Normally i would use a test like Dunn's test here, but because
##         we are using a mannwhitney u test in the main figures calc in python
##         I opted to use the Wilcox test here instead to use the same statistics in both figures

df_total %>% group_by(expan_thres) %>% 
  wilcox_test(CSF ~ diagnosis) %>% 
  arrange(p)

## Simplified comparison - control vs ms
stats_df = df_total %>% 
  mutate(catdiag = ifelse(diagnosis == "IIH", "Ctrl", "MS"))
stats_df %>% group_by(expan_thres) %>% 
  wilcox_test(CSF ~ catdiag) %>% 
  arrange(p)


# %%%%%%%%% 
# Export
# %%%%%%%%%
# Arrange with cowplot
plt_fin = plot_grid(expan_csf_freq_plot, ncol = 2, nrow = 4)

ggsave(filename = "05_plots/06_ValidationCohort/2026_03 - ExpansioDegree and cells in CSF.pdf",
       plot = plt_fin, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")




# Numbers for writing and annotation --------------------------------------
obj_meta %>% 
  filter(total_tcr_freq > 12) %>% 
  select(diagnosis, Sample_Clono) %>% 
  unique %>% 
  group_by(diagnosis) %>% 
  tally



# Column plot - CSF vs PB -------------------------------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% Export .csv file with cell counts in CSF vs PB for stats  %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## I need to set a threshold for this dataset on what constitutes
##  a CSF-biased and a PB-biased clonotype. 
##  A separate threshold is needed as their methods assay specificity and 
##  sensitivity is likely different from ours

summed_counts = obj_meta %>% 
  ## First calculate the pb_tcr_frequency in counts
  mutate(pb_tcr_frequency = total_tcr_freq - csf_tcr_frequency) %>% 
  select(diagnosis, samplenumb, csf_tcr_frequency, pb_tcr_frequency, Sample_Clono) %>% 
  unique() %>% 
  ## Second summarize the total counts in PB and CSF and derive a % split
  summarise(sum_csf = sum(csf_tcr_frequency), sum_pb = sum(pb_tcr_frequency))

## Frequency split in PB and CSf
CSF_counts = summed_counts %>% select(sum_csf) %>% pull
PB_counts = summed_counts %>% select(sum_pb) %>% pull
total_counts = CSF_counts + PB_counts
CSF_split = CSF_counts/total_counts
CSF_split
PB_split = PB_counts/total_counts
PB_split


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% At n > 12 cells 0.230769 fraction of cells in CSF is sig. %%% 
# %%%                        CSF-biased                         %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Calculated in python
## Note that PB-biased are not identifiable until n = 80 cells


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%         The column plots for hyperexpanded cells          %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#For hyperexpanded clones - show compartment-wise expansion
obj_meta_compartents = obj_meta %>% 
  filter(total_tcr_freq > 12) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique 

# Dimension of the reported condition (i.e. we have 1 IIH and 1 relapse without hyperexpanded clones
obj_meta %>% select(diagnosis, samplenumb) %>% unique %>% group_by(diagnosis) %>% tally()
obj_meta_compartents %>% select(diagnosis, samplenumb) %>% unique %>% group_by(diagnosis) %>% tally()

# Order of clones for the eventual column plots
order_clonotypes = obj_meta_compartents %>% arrange(PB) %>% select(Sample_Clono) %>% pull()

plot_df_compartments = obj_meta_compartents %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct")

# To ensure the order of clonotypes are uphold - define the Sample_Clono column as factors with a specific order
plot_df_compartments$Sample_Clono = factor(plot_df_compartments$Sample_Clono, levels = order_clonotypes)
plot_df_compartments$Compartment = factor(plot_df_compartments$Compartment, levels = c("PB", "CSF"))

P1 = ggplot(plot_df_compartments, aes(x = Sample_Clono, y = CSF_pct, fill = Compartment, width=1))+
  geom_col()+
  geom_hline(yintercept = 23.0769)+ #Cacluated in python
  scale_fill_manual(values = list("CSF" = "#4B8AC1", "PB" = "#C06A6A"))+
  theme_pubclean(base_size = 6)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_grid(~diagnosis, scale = "free", space = "free")

P1

# %%%%%%%%% 
# Export
# %%%%%%%%%
ggsave(filename = "05_plots/06_ValidationCohort/2025_12 - column_plot_CSFvPB.pdf",
       plot = P1, 
       device = cairo_pdf, 
       width = 50, 
       height = 21, 
       units = "mm")


# Boxplot - % CSF-biased clones ----------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%   Is the overal frequency of CSF-biased clones larger?    %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#For hyperexpanded clones - show compartment-wise expansion
plot_df_compartments = obj_meta %>% 
  filter(total_tcr_freq > 12) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct")

plot_df = plot_df_compartments %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral))

ggplot(plot_df, aes(x = diagnosis, y = pct_csfbiased))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(dodge.width = 0.8, width = 0.03, size = 1.5)

## Wilcox test
plot_df %>% 
  wilcox_test(pct_csfbiased ~ diagnosis) %>% 
  arrange(p)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%       % CSF-biased clones across expansion thresholds      %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df0 = obj_meta %>% 
  filter(total_tcr_freq > 0) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct") %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral)) %>% 
  mutate(Source = "0")

df5 = obj_meta %>% 
  filter(total_tcr_freq > 5) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct") %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral)) %>% 
  mutate(Source = "5")

df10 = obj_meta %>% 
  filter(total_tcr_freq > 10) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct") %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral)) %>% 
  mutate(Source = "10")

df15 = obj_meta %>% 
  filter(total_tcr_freq > 15) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct") %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral)) %>% 
  mutate(Source = "15")

df20 = obj_meta %>% 
  filter(total_tcr_freq > 20) %>% 
  select(diagnosis, samplenumb, CSF, PB, Sample_Clono) %>% 
  unique %>% 
  pivot_longer(., cols = c("CSF", "PB"), names_to = "Compartment", values_to = "CSF_pct") %>% 
  filter(Compartment == "CSF") %>% 
  mutate(biased_ann = ifelse(CSF_pct > 23.0769, "CSF_biased", "Neutral")) %>% 
  select(diagnosis, samplenumb, Sample_Clono, biased_ann) %>% 
  unique() %>% 
  group_by(diagnosis, biased_ann, samplenumb) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(diagnosis, samplenumb), names_from = biased_ann, values_from = n, values_fill = 0) %>% 
  mutate(pct_csfbiased = 100*CSF_biased / (CSF_biased + Neutral)) %>% 
  mutate(Source = "20")


## Bind the dataframes together
df_total = rbind(df0, df5, df10, df15, df20)

df_total

### %%%%
## Boxplot
### %%%% 

df_total$Source = factor(df_total$Source, levels = c("0", "5", "10", "15", "20"))

csf_biased_pct_plot = ggplot(df_total, aes(x = Source, y = pct_csfbiased+0.05, fill = diagnosis))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(dodge.width = 0.8, width = 0.03, size = 1.5)+
  scale_y_log10(guide = "axis_logticks", limits=c(0.05, 100))+
  ylab("% CSF-biased clones")+
  theme_pubr()
csf_biased_pct_plot

## Statistics. y as as function of x (y ~ x)
### Note: Normally i would use a test like Dunn's test here, but because
##         we are using a mannwhitney u test in the main figures calc in python
##         I opted to use the Wilcox test here instead to use the same statistics in both figures

df_total %>% group_by(Source) %>% 
  wilcox_test(pct_csfbiased ~ diagnosis, p.adjust.method = "BH") %>% 
  arrange(p)

# %%%%%%%%% 
# Export
# %%%%%%%%%
# Arrange with cowplot
plt_fin = plot_grid(csf_biased_pct_plot, ncol = 2, nrow = 4)

ggsave(filename = "05_plots/06_ValidationCohort/2026_08 - pct_CSFbiased_clones.pdf",
       plot = plt_fin, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")





# Correlation with our CSF-biased score -----------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%   Correlation of CSF-biased score and frac.cells in CSF    %%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## For hyperexpanded clones (n > 20)
obj_meta_small = obj_meta %>% 
  filter(total_tcr_freq > 20) %>% 
  group_by(Sample_Clono) %>% 
  mutate(avg_csf_exp_score = mean(csf_features1)) %>% 
  ungroup %>% 
  group_by(samplenumb) %>% 
  mutate(avg_csf_exp_score_donor = mean(csf_features1)) %>% 
  select(diagnosis, samplenumb, avg_csf_exp_score, avg_csf_exp_score_donor, Sample_Clono, csf_tcr_frequency, total_tcr_freq) %>% 
  unique %>% 
  ungroup %>% 
  arrange(desc(csf_tcr_frequency)) 

## Correlation for hyper-expanded clones in PB and/or CSF
# Create scatter plot with regression line and Spearman correlation
P1 = ggscatter(obj_meta_small, x = "csf_tcr_frequency", y = "avg_csf_exp_score",
          add = "reg.line",           # add linear regression line
          size = 1,
          conf.int = TRUE,            # confidence interval
          cor.coef = TRUE,            # show correlation coefficient
          cor.method = "spearman",    # use Spearman correlation
          cor.coef.coord = c(1, 0.3),   # position of the correlation label
          xlab = "csf_tcr_frequency",
          ylab = "avg_csf_exp_score") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_pubr(base_size = 7)

# %%%%%%%%% 
# Export
# %%%%%%%%%
# Arrange with cowplot
plt_fin = plot_grid(P1, ncol = 4, nrow = 6)

ggsave(filename = "05_plots/06_ValidationCohort/2026_03 - SpearmanCorrelation_CSF_fraction.pdf",
       plot = plt_fin, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")

#Double check stats placed on teh figure
x_values = obj_meta_small$csf_tcr_frequency
y_values = obj_meta_small$avg_csf_exp_score

cor.test(x = x_values, y = y_values, method = "spearman")

