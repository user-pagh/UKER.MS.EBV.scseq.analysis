## Immune monitoring using TRB sequencing 
## Author: Nikolaj P Kristensen
## Date: 22-04-2026
## Version 2


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(scales)
library(ggbeeswarm)
library(cowplot)
library(ggrastr)
library(readxl)
library(ggalluvial)

# Load data ---------------------------------------------------------------

#Generate a list of file paths
files_to_read = list.files(path = "02_raw_data/bulk_TCRseq_mixcrouts/e78-80/", pattern="\\.tsv$", full.names = T)

#"list apply" the read.table function with tab delimiter to make a list of data frames
all_files = lapply(files_to_read, function(x) {
  read.table(file = x, 
             sep = "\t", 
             header = TRUE)
})

#bind_rows to merge the data frames. Generate a numerical ID. 
raw = bind_rows(all_files, .id="ID")

#Take the list of paths, and add them to the data data according to the numerical ID 
raw %>% select(ID) %>% unique
file_ids = as_tibble(files_to_read) %>% 
  mutate(ID = row.names(.))

raw = raw %>% left_join(x = ., y = file_ids)

#Convert to tibble and continue processing
raw1 = raw %>% 
  as_tibble() %>% 
  select(-ID) %>% 
  rename("ID" = value) %>% 
  relocate(ID, .before="cloneId") %>% 
  mutate(Library = ID %>% str_remove(., pattern = "...*")) %>% 
  mutate(ID = ID %>% str_remove(string = ., pattern = "02_raw_data/bulk_TCRseq_mixcrouts/e78-80/...-......")) %>% 
  mutate(ID = ID %>% str_remove(string = ., pattern = "-Pt.*")) %>%  
    #What the reges means:
    ### -Pt   matches "Pt"
    ### .    matches any character
    ### *    means "zero or more". In this case, means zero or more of any character after -Pt 
  mutate(ID = ID %>% str_remove(string = ., pattern = "^-")) %>%
    ### ^-   matches - at the start of the strand %>% 
  # mutate(ID = ifelse(str_detect(., "Jurkat"), "Jurkat", ID)) %>% 
  mutate(ID = ID %>% str_remove(string = ., pattern = "-Jurkat_ctrl-M7A7B7.CAGAT_CAGAT.clones_TRB.tsv")) %>%
  mutate(ID = ID %>% str_remove(string = ., pattern = "rkat-")) %>% 
  rename(TRBV = allVHitsWithScore, TRBJ = allJHitsWithScore) %>% 
  mutate(TRBV = TRBV %>% str_remove(string = ., pattern = "\\(.*$")) %>%
  mutate(TRBJ = TRBJ %>% str_remove(string = ., pattern = "\\(.*$")) %>%
  select(ID, TRBV, TRBJ, aaSeqCDR3, uniqueMoleculeFraction, uniqueMoleculeCount)

#Add filtering arguments to remove nonesense chains
## %%%%%%%%%%%%%%
## Important filtering decisions. 
##  I am removing reads with less than 2 UMIs, because the frequency estimation is very difficult with 
##  rare TRB chains. This is mentioned in the TCRseq benchmarking paper, where the threshold is 
##  mentioned as 2-4 "reads" per UMI. The mixcr webpages doesn't say anything about a UMI 
##  requiring multiple reads, so I am convinced to some degree that what they mean is that 
##  each chain needs more than two unique UMIs to be identified to adequately count it. 
## %%%%%%%%%%%%
raw2 = raw1 %>% 
  filter(!str_detect(string = aaSeqCDR3, pattern = "\\*")) %>% 
  filter(!str_detect(string = aaSeqCDR3, pattern = "_")) %>% 
  ## What this filtering step does:
  ## It removes aaSeqCDR3 entries containing \\* or _  
  filter(uniqueMoleculeCount > 2) %>% 
  ## This removes UMI counts below 2 as the frequency here cannot be adequately estimated. 
  mutate(Donor = str_extract(string = ID, pattern = "..")) %>% 
  mutate(Donor = ifelse(Donor == "Ju", "Jurkat", Donor)) %>% 
  mutate(Timepoint = str_extract(string = ID, pattern = "-.*")) %>%
  mutate(Timepoint = str_remove(Timepoint, "-")) %>% 
  replace_na(., replace = list(Timepoint = "Control")) %>% 
  mutate(Clone_handle = paste0("Pt",Donor, "_", aaSeqCDR3))


# Correct for non-synonymous mutations by summarizing the frequency pre aaSeqCDR3 --------
raw2 = raw2 %>% 
  group_by(ID, aaSeqCDR3, Donor, Timepoint, Clone_handle) %>% 
  summarise(uniqueMoleculeFraction = sum(uniqueMoleculeFraction)) %>% 
  ungroup

raw2 %>% select(Donor) %>% unique

raw2$Donor = factor(raw2$Donor, levels = c("14", "17", "19", "18", "20", "21", "Jurkat"))
raw2$Timepoint = factor(raw2$Timepoint, levels = c("D0", "D19", "D37", "Control"))

#Color annotate outliers
## Define a list of CDR3s that you would like to annotate
clone_list = raw2 %>%
  filter(Timepoint == "D37") %>% 
  group_by(Donor) %>% 
  top_n(., n = 5, wt = uniqueMoleculeFraction) %>%  #Outlier = top 5
  ungroup() %>% 
  select(Clone_handle) %>% 
  pull

raw2 = raw2 %>% 
  mutate(COLOR = ifelse(Clone_handle %in% clone_list, "red", "black"))


## For those that you highlight, make sure that it is visible that they are undetected at D0
## First, generate the highlight dataframe and complete the missing combinations of timepoint and clone handle 
## for each donor. 
raw3 = raw2 %>% 
  filter(Timepoint != "Control")

raw3$Timepoint = factor(raw3$Timepoint, levels = c("D0", "D19", "D37"))


# Define d19 expanded clones -------------------
####
# Increasing frequency is here defined by increasing fold change
# comparing D0 vs D19, D0 vs D37, and D19 vs D37
###
raw3 = raw2 %>% 
  filter(Timepoint != "Control")

# Kilian: Please exclude D37
raw3v1 = raw3 %>% filter(Timepoint != "D37")

#Backbone_df contains all combinations of clone_handle and timepoint. 
backbone_df = raw3v1  %>% 
  group_by(Donor) %>% 
  expand(Timepoint, Clone_handle) %>%  #This expands Timepoint and Clone_handle combinations
  ungroup %>% 
  mutate(aaSeqCDR3 = str_extract(string = Clone_handle, pattern = "(?<=_).*"))

#Regex explained:
#(?<=-)   Checks behind to see whether "-" is present
#.*       Means any number(*) of characters(.)

###
# Complete the backbone_df with UMI values
### 
value_df = raw3v1 %>% 
  select(Clone_handle, aaSeqCDR3, Donor, Timepoint, uniqueMoleculeFraction)

# Left join
complete_df = backbone_df %>% left_join(x = ., y = value_df) %>% 
  replace_na(replace = list(uniqueMoleculeFraction = 1e-6))

# Make a wide dataframe for each clone handle and time point
clone_list2 = 
  complete_df %>% pivot_wider(., id_cols = c("Donor", "Clone_handle"), names_from = Timepoint, values_from = uniqueMoleculeFraction)  %>% 
  mutate(FC0_19 = D19 / D0) %>% 
  group_by(Donor) %>% 
  top_n(n = 5, wt = FC0_19) %>% 
  ungroup %>% 
  select(Clone_handle) %>% 
  unique %>% 
  pull


# Export d19 expanded clones -------------------
clone_export = 
  complete_df %>% pivot_wider(., id_cols = c("Donor", "Clone_handle", "aaSeqCDR3"), names_from = Timepoint, values_from = uniqueMoleculeFraction)  %>% 
  mutate(FC0_19 = D19 / D0) %>% 
  filter(FC0_19 > 1)

clone_export %>% write_csv(., "06_csv_outputs/Clones/2026_4_BLCL_reactive_d19.csv")


# Numbers for writing -----------------------------------------------------
raw3v1 %>% 
  group_by(Timepoint, Donor) %>% 
  tally %>% 
  ungroup %>% 
  pivot_wider(names_from = Timepoint, values_from = n)



# Plot d19 expanded clones -------------------
## For those that you highlight, make sure that it is visible that they are undetected at D0
## First, generate the highlight dataframe and complete the missing combinations of timepoint and clone handle 
## for each donor. 
raw3v1$Timepoint = factor(raw3v1 $Timepoint, levels = c("D0", "D19"))

backbone_df = raw3v1 %>% 
  filter(Clone_handle %in% clone_list2) %>% 
  group_by(Donor) %>% 
  expand(Timepoint, Clone_handle) %>%  #This expands Timepoint and Clone_handle combinations
  ungroup %>% 
  mutate(aaSeqCDR3 = str_extract(string = Clone_handle, pattern = "(?<=_).*"))
#Regex explained:
#(?<=-)   Checks behind to see whether "-" is present
#.*       Means any number(*) of characters(.)

# Values data frame containing uniqueMoleculeFraction and a series of feature columns for left joining
value_df = raw3v1 %>% 
  select(Clone_handle, aaSeqCDR3, Donor, Timepoint, uniqueMoleculeFraction) %>% 
  mutate(COLOR = "Expanded")

# left join
highlight_df = backbone_df %>% left_join(x = ., y = value_df) %>% 
  replace_na(replace = list(uniqueMoleculeFraction = 1e-6))

## Second, remove those highlight entries from the original dataframe
plot_df1 = raw3v1 %>% 
  mutate(Timepoint = Timepoint %>% as.character()) %>% 
  mutate(Donor = Donor %>% as.character()) %>% 
  filter(!Clone_handle %in% clone_list2) %>% 
  filter(Timepoint != "Control") %>% 
  mutate(COLOR = "Background") %>% 
  select(Clone_handle, aaSeqCDR3, Donor, Timepoint, uniqueMoleculeFraction, COLOR) 

## Third, bind the highlight_df 
plot_df2 = rbind(plot_df1, highlight_df) %>% 
  replace_na(replace = list(COLOR = "Expanded")) ## Corresponds to the missing levels for highlighted dots

## Further annotation
patient_list = c("4", "11", "12", "14", "17", "20")
control_list = c("1", "5", "8", "18", "19", "21")
plot_df2 = plot_df2 %>% 
  mutate(catpat = ifelse(Donor %in% patient_list, "RRMS",
                         ifelse(Donor %in% control_list, "Control", "error")))

#Remake the factor levels
plot_df2$Donor = factor(plot_df2$Donor, levels = c("14", "17",  "18", "19", "20", "21"))
plot_df2$Timepoint = factor(plot_df2$Timepoint, levels = c("D0", "D19"))

P3 = plot_df2 %>% 
  arrange(COLOR) %>%   #serves to arrange dot order
  ggplot(., aes(x = Timepoint, y = 100*uniqueMoleculeFraction, 
                color = COLOR, group = aaSeqCDR3, 
                alpha = COLOR, linewidth = COLOR, size = COLOR))+
  ggrastr::rasterise(geom_line(), dpi=400)+
  ggrastr::rasterise(geom_point(), dpi=400)+
  geom_hline(yintercept = 5e-4)+
  ylab("% of UMIs")+
  scale_y_log10(guide="axis_logticks")+
  scale_color_manual(values = c(
    "Expanded" = "red", 
    "Background" = "black"
  ))+
  scale_alpha_manual(values = c(
    "Expanded" = 0.5, 
    "Background" = 0.1
  ))+
  scale_linewidth_manual(values = c(
    "Expanded" = 0.25, 
    "Background" = 0.1
  ))+
  scale_size_manual(values = c(
    "Expanded" = 1, 
    "Background" = 0.2
  ))+
  theme_pubr(base_size = 8)+
  theme(     panel.background = element_rect(fill = NA, colour = NA),
             plot.background  = element_rect(fill = NA, colour = NA))+
  facet_wrap(catpat~Donor, axis.labels = "all_x", nrow = 1)

P3

# Export plot -------------------------------------------------------------
to_export_plot = plot_grid(P3, 
                      ncol = 2, rel_widths = c(1,1), 
                      nrow = 3, rel_heights = c(1,2,3))

ggsave(filename = "05_plots/04_TCRseq_visualized/2026_04 - Longitudinal per donor - expanding clones.pdf",
       plot = to_export_plot, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")





# Multimer+ clones -----------------------------------
ag_ann = readxl::read_xlsx(path = "03_annotation/11_multimer_experiment_ann/clone_x_multimer_assignments.xlsx") %>% 
  filter(!is.na(clone_handle), !is.na(Virus)) %>% 
  select(aaSeqCDR3 = IR_VDJ_1_junction_aa, Virus, Antigen, HLA, clone_handle) %>% 
  mutate(donor = clone_handle %>% str_extract(., "^Pt\\d+(?=-)")) %>% 
  ## Regex explanation
  #  ^Pt means start of the string with "Pt"
  #  \\d+ → matches one or more digits
  #  (?=-) → ensures the match is followed by "-"
  unique %>% 
  mutate(Clone_handle = paste0(donor, "_", aaSeqCDR3)) %>% 
  select(Clone_handle, Virus, Antigen)


#Left join annotations
raw3 = raw2 %>% 
  filter(Timepoint != "Control") %>% 
  filter(Timepoint != "D37") 
raw4 = raw3 %>% left_join(x = ., y = ag_ann)

#All matched clones with epitope specificity
clone_list3 = raw4 %>% filter(!is.na(Virus)) %>% select(Clone_handle) %>% unique %>% pull


# Alluvial plot  -----------------------------------------------------------------
## See the last example on vaccination surveys here: https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
plot_df = raw4 %>%
  mutate(Timepoint = Timepoint %>% as.character()) %>% 
  filter(Timepoint != "D37", Clone_handle %in% clone_list3) %>% 
  select(-ID) %>% 
  pivot_wider(., names_from = Timepoint, values_from = uniqueMoleculeFraction)  %>% 
  replace_na(replace = list(D0 = 1e-6, D19 = 1e-6, D37 = 1e-6)) %>% 
  mutate(FC0_19 = D19 / D0) %>% 
  mutate(catexpand = ifelse(FC0_19 > 1, "Expanding", "Background")) 


allu_df1 <- plot_df %>%
  distinct(Virus, Clone_handle) %>%
  count(Virus, name = "count") %>%
  mutate(source = "All_clones")

allu_df1 %>% summarise(sum = sum(count))

allu_df2 <- plot_df %>%
  filter(catexpand == "Background") %>%
  distinct(Virus, Clone_handle) %>%
  count(Virus, name = "count") %>%
  mutate(source = "Background")
allu_df2 %>% summarise(sum = sum(count))

allu_df3 <- plot_df %>%
  filter(catexpand == "Expanding") %>%
  distinct(Virus, Clone_handle) %>%
  count(Virus, name = "count") %>%
  mutate(source = "Expanding")
allu_df3 %>% summarise(sum = sum(count))

  #x = input_distribution
  #y = frac
  #stratum = virus
  #alluvum = subject

allu_df4 = rbind(allu_df1, allu_df2, allu_df3 )

#Normalize
allu_df4 <- allu_df4 %>%
  group_by(source) %>%
  mutate(frac = count / sum(count)) %>%
  ungroup()

allu_df4 <- transform(allu_df4,
                          response = factor(Virus, rev(levels(Virus))))

P5 = ggplot(allu_df4,
       aes(x = source, stratum = Virus, alluvium = Virus,
           y = frac,
           fill = Virus, label = Virus)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 1.5) +
  ylab("Fraction of clones")+
  theme(legend.position = "none")+
  theme_pubr(base_size = 7)

P5

# Export  -----------------------------------------------------------------
to_export_plot = plot_grid(P5, ncol = 4) %>% 
  plot_grid(., nrow = 5)

ggsave(filename = "05_plots/04_TCRseq_visualized/2026_04 - Longitudinal per donor - alluvial.pdf",
       plot = to_export_plot, 
       device = cairo_pdf, 
       width = 210,   #A4
       height = 297, 
       units = "mm")


# Across Ag --------------------------------------------------------
####
# Increasing frequency is here defined by increasing fold change
# comparing D0 vs D19 only
###
raw3 = raw2 %>% 
  filter(Timepoint != "Control")

#Left join annotations
raw4 = raw3 %>% left_join(x = ., y = ag_ann)

# Make a wide dataframe for each clone handle and time point
plot_df = raw4 %>%
  select(-ID) %>% 
  filter(!is.na(Virus)) %>% 
  pivot_wider(., names_from = Timepoint, values_from = uniqueMoleculeFraction)  %>% 
  replace_na(replace = list(D0 = 1e-6, D19 = 1e-6, D37 = 1e-6)) %>% 
  mutate(FC0_19 = D19 / D0, 
         FC0_37 = D37 / D0, 
         FC19_37 = D37 / D19) %>% 
  mutate(catexpand = ifelse(FC0_19 > 1, "Expanding", "Background")) 

plot_df$catexpand = factor(plot_df$catexpand, levels = c("Expanding", "Background"))

plot_df$Antigen = factor(plot_df$Antigen, levels = c("IE1","pp65", "pp50", #CMV
                                                     "BaRF0", "BaRF1", "BFRF3", "BORF2", "BMLF1", "BRLF1", "BZLF1", 
                                                     "EBNA1", "EBNA3", "EBNA6", "LMP1", "LMP2",   #EBV
                                                     "MP",  #FLU-A
                                                     "VP1", #NWV
                                                     "S", "rep",    #Sars-cov-2
                                                     "ICP0",        #HHV-1
                                                     "RIR2"         #VZV
                                                     ))
 


# Make a quick stacked bargraph visualizing counts
P1 = plot_df %>% 
  ggplot(., aes(x = Antigen, fill = catexpand))+
  geom_bar(color = "black")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Number of clones")+
  theme_pubr(base_size = 7)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill = guide_legend(title = " "))
P1

# normalized
P2 = plot_df %>% 
  ggplot(., aes(x = Antigen, fill = catexpand))+
  geom_bar(position="fill", color = "black")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Number of clones")+
  theme_pubr(base_size = 7)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill = guide_legend(title = " "))

# Per donor - count
plot_df3 = plot_df %>% 
  group_by(Antigen, catexpand, Donor) %>% 
  tally() %>% 
  ungroup %>% 
  pivot_wider(., names_from = catexpand, values_from = n) %>% 
  replace_na(., replace = list(Expanding=0, Background=0)) %>% 
  # mutate(total = Expanding+Background,
  #        Expanding = 100*Expanding/total, 
  #        Background = 100*Background/total) %>% 
  select(Donor, Antigen, Expanding, Background) %>% 
  pivot_longer(cols = c("Expanding", "Background"), names_to = "catexpand", values_to = "counts")  

plot_df3$catexpand = factor(plot_df3$catexpand, levels = c("Expanding", "Background"))

P3 = plot_df3 %>% 
  ggplot(., aes(x = Antigen, y = counts, fill=catexpand))+
  geom_boxplot(outlier.shape=NA,  aes(fill = catexpand))+
  geom_point(position=position_jitterdodge())+
  ylab("Number of clones")+
  theme_pubr(base_size = 4)+
  guides(fill = guide_legend(title = " "))

# Per catpat
patient_list = c(14, 17, 20)
control_list = c(18, 19, 21)

plot_df3 = plot_df3 %>% 
  mutate(catpat = ifelse(Donor %in% patient_list, "MS", "NIC"))

P4 = plot_df3 %>% 
  filter(catexpand == "Expanding") %>% 
  ggplot(., aes(x = Antigen, y = counts, fill=catpat))+
  geom_boxplot(outlier.shape=NA,  aes(fill = catpat))+
  geom_point(position=position_jitterdodge())+
  scale_fill_brewer(palette = "Set2")+
  labs(caption = "Of expanding clones only") +
  ylab("Number of clones")+
  theme_pubr(base_size = 4)+
  guides(fill = guide_legend(title = " "))


# Writing -----------------------------------------------------------------
plot_df %>% 
  group_by(Antigen) %>% 
  tally() %>% 
  arrange(desc(n))



# Export  -----------------------------------------------------------------
to_export_plot = plot_grid(P1, P2, P3, P4, nrow = 2) %>% 
  plot_grid(., nrow = 3)

ggsave(filename = "05_plots/04_TCRseq_visualized/2026_04 - Longitudinal per donor - Across antigen.pdf",
       plot = to_export_plot, 
       device = cairo_pdf, 
       width = 210,   #A4
       height = 297, 
       units = "mm")








