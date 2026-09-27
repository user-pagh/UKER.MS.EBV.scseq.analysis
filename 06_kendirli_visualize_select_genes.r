## Visualize select genes from Kendirli et l
## auhtor: Nikolaj pagh kristensen
## Date: 10-8-2026


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(gprofiler2)

# Load data ---------------------------------------------------------------
## -- Blood vs Meninges -- 
BM_raw = read_xlsx(path = "03_annotation/06_CRISPR_sets/Kendirli_Meniges_vs_blood_validation.xlsx")

## Load the CSF-biased score
csf_permissive = read_csv(file = "06_csv_outputs/DEG/CSF_permisiveness_score.csv") %>% 
  select(names) %>% pull()

## Load data ---------------------------------------------------------------

## -- Blood vs Meninges --
BM_raw = read_xlsx(
  path = "03_annotation/06_CRISPR_sets/Kendirli_Meniges_vs_blood_validation.xlsx"
)

## Load the CSF-biased score
csf_permissive = read_csv(
  file = "06_csv_outputs/DEG/CD8_CSF_biased_filtered_DEGs.csv"
) %>%
  select(names) %>%
  pull()


## Convert human genes to rat orthologs ------------------------------------
rat_mapping <- gorth( query = csf_permissive, source_organism = "hsapiens", target_organism = "rnorvegicus" )

## Extract rat gene symbols
csf_permissive_rat <- rat_mapping %>% filter(!is.na(ortholog_name)) %>% pull(ortholog_name) %>% unique()


## Load data ---------------------------------------------------------------

## -- Blood vs Meninges --
BM_raw = read_xlsx(
  path = "03_annotation/06_CRISPR_sets/Kendirli_Meniges_vs_blood_validation.xlsx"
)

## Select genes that overlapped with the CSF biased score
BM_raw = BM_raw %>% 
  filter(id %in% csf_permissive_rat)

BM_raw


# Column plot -------------------------------------------------------------
Kendirli_plot = ggplot(BM_raw, aes(x = id %>% reorder(., `neg|p-value`), y = -log(`neg|p-value`)))+
  geom_col()+
  geom_hline(yintercept = -log(0.05), colour = "red")+
  theme_pubr()
Kendirli_plot

ggsave(filename = "05_plots/06_ValidationCohort/2026_08_Kendirli_Mengiges_vs_blood.pdf",
       plot = Kendirli_plot, 
       device = cairo_pdf, 
       width = 105, 
       height = 50, 
       units = "mm")
  