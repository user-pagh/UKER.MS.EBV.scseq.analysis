
#TCRseq of T-LCL cultures
#2025.02.17


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)
library(readxl)
library(RColorBrewer)

# Working path ------------------------------------------------------------

setwd("W:/mikrobiologie/archiv/AG_Schober/Nikolaj/R.files/MS.Schober")


# TRB Data --------------------------------------------------------------------
Dat1 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M1B1_Patient1.ATCAC_ATCAC.clones_TRB.tsv") %>% mutate(Donor = "1-IIH")
Dat2 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M2B2_Patient4.CGATG_CGATG.clones_TRB.tsv") %>% mutate(Donor = "4-MS")
Dat3 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M3B3_Patient5.TTAGG_TTAGG.clones_TRB.tsv") %>% mutate(Donor = "5-CU")
Dat4 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M4B4_Patient8.TGACC_TGACC.clones_TRB.tsv") %>% mutate(Donor = "8-IIH")
Dat5 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M5B5_Patient11.ACAGT_ACAGT.clones_TRB.tsv") %>% mutate(Donor = "11-MS")
Dat6 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M6B6_Patient12.GCCAA_GCCAA.clones_TRB.tsv") %>% mutate(Donor = "12-MS")
Dat7 = read_tsv(file = "02_raw_data/bulk_TCRseq_mixcrouts/e67/M7B7_Jurkat.CAGAT_CAGAT.clones_TRB.tsv") %>% mutate(Donor = "Jurkat")

raw = rbind(Dat1, Dat2, Dat3, Dat4, Dat5, Dat6, Dat7)  
raw %>% select(Donor) %>% unique
#Relabel the "5-MS /HC" sample to "5-CU"
raw = raw %>% 
  mutate(Donor = ifelse(Donor == "5-MS / HC",  "5-CU", Donor))

#How many differen CDR3s?
raw %>% 
  select(Donor, aaSeqCDR3) %>% 
  unique %>% 
  group_by(Donor) %>% 
  tally() %>% 
  ggplot(., aes(x = Donor %>% reorder(., n), y = n))+
  geom_col()+
  geom_text( aes(x = Donor %>% reorder(., n), y = n+750, label = n))+
  theme_pubr(base_size = 12)+
  labs(y = "Number of unique TRB CDR3 amino acid sequences", x = "Samples")

raw_wide = raw %>% 
  pivot_wider(names_from = "Donor", values_from = "uniqueMoleculeFraction") %>% 
  select(aaSeqCDR3, `1-IIH`, `4-MS`, `5-CU`, `8-IIH`, `11-MS`, `12-MS`, Jurkat) %>% 
  replace(is.na(.), 0)

P1 = ggplot(raw_wide, aes(x = `1-IIH`, y = `4-MS`))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: 4-MS, which mainly grew T cells)")+
  theme_pubr()

P2= ggplot(raw_wide, aes(x = `1-IIH`, y = `5-CU`))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: 5-MS/HC, which mainly grew B cells)")+
  theme_pubr()

P3= ggplot(raw_wide, aes(x = `1-IIH`, y = `8-IIH`))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: 8-IIH, which mainly grew B cells)")+
  theme_pubr()


P4= ggplot(raw_wide, aes(x = `1-IIH`, y = `11-MS`))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: 11-MS, which mainly grew T cells)")+
  theme_pubr()

P5= ggplot(raw_wide, aes(x = `1-IIH`, y = `12-MS`))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: 12-MS, which mainly grew T cells)")+
  theme_pubr()

P6= ggplot(raw_wide, aes(x = `1-IIH`, y = Jurkat))+
  geom_point(size = 2.5, alpha = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Fraction of UMIs\n(sample: 1-IIH, which mainly grew B cells)",
       y = "Fraction of UMIs\n(sample: Jurkat control)")+
  theme_pubr()

plot_grid(P1, P2, P3, P4, P5, P6, ncol = 3, nrow = 2)

# Visualize distribution of UMIs ------------------------------------------
raw$Donor = factor(raw$Donor, levels = c("1-IIH", "5-CU", "8-IIH", "4-MS", "11-MS", "12-MS"))


#Most patients have dominating clones
raw %>% 
  filter(Donor=="4-MS") %>%
  mutate(Donor="RRMS-4") %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeCount) %>%  
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_count = sum(uniqueMoleculeCount)) %>% 
  ggplot(., aes(x = aaSeqCDR3 %>% reorder(., umi_count) , y = umi_count, color = Donor))+
  geom_point(size = 1)+
  geom_text(data = . %>% filter(umi_count > 5000), aes(label = aaSeqCDR3), hjust = 1)+
  scale_x_discrete(expand=c(0.03,0.03))+
  labs(y = "TRB, UMI count", x = "6028 TRB sequences", caption = "Day 37 after EBV transformation")+
  theme_pubr()+
  theme(axis.text.x = element_blank())+
  facet_grid(~Donor, scale = "free_x", space="free_x")

raw %>% 
  filter(Donor=="11-MS") %>%
  mutate(Donor="RRMS-11") %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeCount) %>%  
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_count = sum(uniqueMoleculeCount)) %>% 
  ggplot(., aes(x = aaSeqCDR3 %>% reorder(., umi_count) , y = umi_count, color = Donor))+
  geom_point(size = 1)+
  geom_text(data = . %>% filter(umi_count > 250), aes(label = aaSeqCDR3), hjust = 1)+
  scale_x_discrete(expand=c(0.03,0.03))+
  # scale_y_continuous(limits = c(0,30000))+
  scale_y_continuous(limits = c(0,5000))+
  labs(y = "TRB, UMI count", x = "6028 TRB sequences", caption = "Day 37 after EBV transformation")+
  theme_pubr()+
  theme(axis.text.x = element_blank())+
  facet_grid(~Donor, scale = "free_x", space="free_x")


#Controls follow a different pattern
raw %>% 
  filter(Donor=="5-CU") %>%
  mutate(Donor="5-CU") %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeCount) %>%  
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_count = sum(uniqueMoleculeCount)) %>% 
  ggplot(., aes(x = aaSeqCDR3 %>% reorder(., umi_count) , y = umi_count, color = Donor))+
  geom_point(size = 1)+
  geom_text(data = . %>% filter(umi_count > 5000), aes(label = aaSeqCDR3), hjust = 1)+
  scale_x_discrete(expand=c(0.03,0.03))+
  # scale_y_continuous(limits = c(0,30000))+
  labs(y = "TRB, UMI count", x = "999 TRB sequences", caption = "Day 37 after EBV transformation")+
  theme_pubr()+
  theme(axis.text.x = element_blank())+
  facet_grid(~Donor, scale = "free_x", space="free_x")

#Controls follow a different pattern
raw %>% 
  filter(Donor=="8-IIH") %>%
  mutate(Donor="8-NID") %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeCount) %>%  
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_count = sum(uniqueMoleculeCount)) %>% 
  ggplot(., aes(x = aaSeqCDR3 %>% reorder(., umi_count) , y = umi_count, color = Donor))+
  geom_point(size = 1)+
  geom_text(data = . %>% filter(umi_count > 5000), aes(label = aaSeqCDR3), hjust = 1)+
  scale_x_discrete(expand=c(0.03,0.03))+
  # scale_y_continuous(limits = c(0,30000))+
  labs(y = "TRB, UMI count", x = "725 TRB sequences", caption = "Day 37 after EBV transformation")+
  theme_pubr()+
  theme(axis.text.x = element_blank())+
  facet_grid(~Donor, scale = "free_x", space="free_x")


raw %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeCount) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_count = sum(uniqueMoleculeCount)) %>% 
  group_by(umi_count) %>% 
  tally() %>% 
  filter(n < 100) %>% 
  ggplot(., aes(x = umi_count %>% reorder(., -n) , y = n))+
  geom_col()+
  labs(y = "Members of each UMI bin", x = "Binds from 1 to 100 UMIs")+
  scale_x_discrete(expand=c(0.03,0.03))+
  #scale_y_continuous(limits = c(0,0.3))+
  theme_pubr()+
  theme(axis.text.x = element_blank())

raw %>% 
  filter(!Donor=="Jurkat", uniqueMoleculeCount > 2) %>%
  select(Donor, aaSeqCDR3, uniqueMoleculeFraction) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_frac = sum(uniqueMoleculeFraction)) %>% 
  ggplot(., aes(x = aaSeqCDR3 %>% reorder(., umi_frac) , y = umi_frac, color = Donor))+
  geom_point()+
  scale_x_discrete(expand=c(0.03,0.03))+
  labs(y = "Fraction of UMIs for each TRB chain", x = "TRB chains")+
  #scale_y_continuous(limits = c(0,0.3))+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_blank(),
        legend.position = "None")+
  facet_wrap(~Donor, scale = "free_x")


# Top clones --------------------------------------------------------------
raw %>% 
  filter(!Donor=="Jurkat", uniqueMoleculeCount > 2) %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeFraction) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_frac = sum(uniqueMoleculeFraction)) %>% 
  group_by(Donor) %>% 
  top_n(n = 5, wt = umi_frac) %>% 
  ungroup() %>% 
  ggplot(., aes(x = Donor, y = umi_frac))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  theme_pubr(base_size = 12)+
  labs(title = "Top 5 clones from each donor", y = "Fraction of UMIs")


raw %>% 
  filter(!Donor=="Jurkat", uniqueMoleculeCount > 2) %>% 
  select(Donor, aaSeqCDR3, uniqueMoleculeFraction) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  summarise(umi_frac = sum(uniqueMoleculeFraction)) %>% 
  group_by(Donor) %>% 
  top_n(n = 3, wt = umi_frac) %>% 
  ungroup() %>% 
  writexl::write_xlsx(x = ., path = "06_csv_outputs/e67_TCRseq/e67_TCRseq_top3clones.xlsx")


# Load in the information from Scanpy -------------------------------------
## Get a new dataframe with clone_ids and cell counts
## 1) Remake airr_export but for both CD4 and CD8 T cells
## 2) Include a label like "CD4" or "CD8"
# airr_raw = read_csv(file = "06_csv_outputs/TCR_selection/2025_02_AllPaired_with_marker_information.csv") 
airr_raw = read_csv(file = "06_csv_outputs/airr_exports/202507_airrexport_allcells.csv")
airr_raw

#scseq data
airr_select = airr_raw %>% 
  select(Donor = donor, clone_handle, Compartment, clone_count = type_x_clone_count, 
         aaSeqCDR3 = IR_VDJ_1_junction_aa) %>% 
  unique %>% 
  pivot_wider(id_cols = c("Donor", "clone_handle", "aaSeqCDR3"), names_from = "Compartment", values_from = "clone_count") %>% 
  replace_na(., replace = list(PB = 0, CSF = 0)) %>% 
  unique %>% 
  mutate(total = CSF+PB) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  top_n(n = 1, wt = total) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  ungroup()

#bulk data
raw1 = raw %>% 
  select(Donor, uniqueMoleculeFraction, aaSeqCDR3) %>% 
  unique %>% 
  mutate(Donor = ifelse(Donor == "1-IIH", "Pt1", 
                        ifelse(Donor == "4-MS", "Pt4", 
                               ifelse(Donor == "5-CU", "Pt5", 
                                      ifelse(Donor == "8-IIH", "Pt8", 
                                             ifelse(Donor == "11-MS", "Pt11", 
                                                   ifelse(Donor == "12-MS", "Pt12", 
                                                          ifelse(is.na(Donor), "Jurkat", "Error")))))))) %>% 
  replace_na(., replace = list(Donor = "Jurkat")) %>% 
  group_by(Donor, aaSeqCDR3) %>% 
  top_n(n = 1, wt = uniqueMoleculeFraction) %>% 
  unique %>% 
  ungroup()

airr_merged = airr_select %>% left_join(x = ., y = raw1) %>% 
  replace_na(., replace = list(uniqueMoleculeFraction = 0)) %>% 
  mutate(detection_limit = ifelse(uniqueMoleculeFraction == 0, "Undetected", "d37 of LCL-culturing")) %>% 
  ungroup() %>% 
  filter(Donor %in% c("Pt1", "Pt4", "Pt5", "Pt8", "Pt11", "Pt12"))

ggplot(airr_merged, aes(x = CSF, y = PB, color = detection_limit, size=uniqueMoleculeFraction))+
  geom_point()+
  scale_color_manual(values = c(
    "Undetected" = "grey50",
    "d37 of LCL-culturing" = "red"
  ))+
  labs(y = "Cell count in sc-seq (PB)", x = "Cell count in sc-seq (CSF)")+
  scale_size(range = c(0.1, 6), labels = c(0.001, 0.01, 0.1, 0.3), breaks = c(0.001, 0.01, 0.1, 0.3))+
  theme_bw(base_size = 12)+
  facet_wrap(~Donor, scales = "free")

## Export top 3 clones from each donor Pt1, Pt11, Pt12, Pt4
LCL_exported = airr_merged %>% filter(Donor %in% c("Pt1", "Pt11", "Pt12", "Pt4")) %>% 
  filter(detection_limit != "Undetected") %>% 
  group_by(Donor) %>% 
  top_n(n = 408, wt = uniqueMoleculeFraction) %>% 
  ungroup() %>% 
  unique %>% 
  group_by(Donor) %>% 
  arrange(uniqueMoleculeFraction %>% desc()) %>% 
  mutate(rank = dense_rank(desc(uniqueMoleculeFraction))) %>% 
  mutate(cOverlap = ifelse(CSF >0 & PB > 0, "shared", 
                           ifelse(CSF > 0, "CSF", "PB"))) %>% 
  mutate(csf_pct = CSF / (PB+CSF)*100) %>% 
  mutate(pb_pct = PB / (PB+CSF)*100)

LCL_exported %>% 
  write_csv(., file = "06_csv_outputs/airr_exports/202507_LCL_clones_export.csv")

#Association between cell count and uniqueMoleculeFraction

