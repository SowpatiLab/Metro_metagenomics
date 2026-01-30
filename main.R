# R packages required to run this code. 
# Refer SessionInfo file to check package version information
library(igraph)
library(tibble)
library(tidyverse)
library(dplyr)
library(janitor)
library(reshape2)
library(RColorBrewer)
library(ALDEx2)
library(abdiv)
library(ggalt)
library(paletteer)
library(ggpubr)
library(tidyplots)
library(ggsankey)
library(cowplot)
library(vegan)
library(ecodist)
library(UpSetR)
library(ggplotify)

# Function to generate top 10 id per city
# Input: A dataframe with id, normalized_count and city_state as column
# Output: 1. Dotplot of top 10 ids across cities
#          2. List of city specific stacked bar plot showing monthly abundance fraction
bar_plot_top10 <- function(dataframe) 
{
  top_all = dataframe %>%
    group_by(city_state, id) %>%
    summarize(sum_normalized_count = sum(normalized_count)) %>% ungroup() %>%
    group_by(id) %>%
    mutate(sum_sum_normalized_count = sum(sum_normalized_count)) %>% ungroup() %>%
    arrange(desc(sum_sum_normalized_count)) %>%
    mutate(id = factor(id, levels = unique(id))) %>% 
    group_by(city_state) %>%
    top_n(10,wt = sum_normalized_count) %>% ungroup() %>% droplevels()
  
  top_all_20 = dataframe %>%
    group_by(city_state, id) %>%
    summarize(sum_normalized_count = sum(normalized_count)) %>% ungroup() %>%
    group_by(city_state) %>%
    top_n(20,wt = sum_normalized_count) %>% ungroup() %>%
    dplyr::select(city_state,id) %>%
    mutate(present = "yes") %>%
    filter(id %in% top_all$id)
  
  top_all_color = setNames(paletteer_d("ggsci::default_igv")[1:length(levels(top_all$id))],levels(top_all$id))
  
  table_plot = top_all %>%
    ggplot(aes(y = (id), x = city_state)) +
    geom_point(size = 5, shape = 21, aes(fill = id), stroke = 0.75) +
    geom_point(data = top_all_20, size = 5, shape = 21, aes(color = present), stroke = 0.75) +
    theme_bw(25) +
    scale_fill_manual(values = top_all_color, name = "") +
    scale_color_manual(values = c("black", "white"), name = "") +
    scale_x_discrete(position = "top", labels = c("Delhi" = "D","Kolkata" = "K","Chennai" = "C","Mumbai" = "M")) +
    theme(axis.text.x=element_text(angle=0, hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_y_discrete(limits = rev, position = "right", label = function(x) str_trunc(x, width = 25)) +
    xlab("") +
    ylab("") 
  
  city_top = list()
  for (city in levels(as.factor(dataframe$city_state)))
  {
    top = top_all %>%
      filter(city_state == city) %>%
      arrange(desc(sum_normalized_count)) %>%
      mutate(id = factor(id, levels = id))
    
    city_top[[city]] = dataframe %>% 
      filter(city_state == city) %>%
      group_by(month_label) %>%
      mutate(fraction_normalized_count = normalized_count/ sum(normalized_count)) %>%
      ungroup() %>%
      filter(id %in% top$id) %>%
      mutate(id = factor(id, levels = levels(top$id))) %>%
      ggplot(aes(x = month_label, y = fraction_normalized_count, fill = fct_inorder(id), group = id)) +
      geom_bar(position="stack", stat="identity") +
      scale_fill_manual(values = top_all_color, name = "") +
      scale_x_discrete(drop = F) +
      theme_bw(25) +
      labs(x = "", y = "", title = city) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.text.x=element_text(angle=90, vjust=0.3, size = 20),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      guides(fill=guide_legend(ncol =1))
  }
  return(list(table_plot, city_top))
}

# Function to generate Alpha diversity (richness, Simpson, Shannon) plot
# Input: A dataframe with month_label, normalized_count and city_state as column
# Output: ggplot object of points faceted by diversity metric
alpha_diversity_plot <- function(dataframe)
{
  dataframe %>%
    group_by(month_label, city_state) %>%
    summarise_at(vars(normalized_count), alpha_diversities) %>%
    dplyr::select(month_label, city_state, richness, simpson, shannon) %>%
    melt(id.vars = c('month_label', 'city_state')) %>%
    ggplot(aes(x = as.factor(month_label), y = value, fill = city_state, group = city_state)) +
    geom_point(size = 5, shape = 21, show.legend = T, alpha = 1) +
    geom_xspline(size = 1, aes(color = city_state), show.legend = T, alpha = 0.7) +
    scale_color_brewer(palette = 'Dark2', name = "City") +
    scale_fill_brewer(palette = 'Dark2', name = "City") +
    theme_minimal(25)  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    labs(x="Collection month", y="", title="") +
    facet_wrap(~variable, scales = "free_y", nrow = 3)
}

# Function to generate PCoA plot based on Bray-Curtis distance
# Input: A dataframe with id, fileName (sample identifier), normalized_count and city_state as column
# Output: ggplot object PCoA plot using first two PcoA axes
bray_curtis_plot <- function(dataframe)
{
  #https://microbiome.github.io/course_2021_radboud/beta-diversity.html
  dataframe_wider = dataframe %>%
    pivot_wider(id_cols = id, names_from = fileName, values_from = normalized_count, values_fill = 0) %>% 
    column_to_rownames(var = 'id')
  bray_curtis_dist <- vegdist(t(dataframe_wider), method = "bray")
  bray_curtis_pcoa <- pco(bray_curtis_dist)
  eigenvalues = bray_curtis_pcoa$values
  percentVar = (eigenvalues/ sum(eigenvalues)) * 100
  bray_curtis_tibble = tibble(pcoa1 = bray_curtis_pcoa$vectors[,1],
                              pcoa2 = bray_curtis_pcoa$vectors[,2]) %>%
    mutate(fileName = rownames(bray_curtis_pcoa$vectors)) %>%
    left_join(dataframe)
  bray_curtis_tibble_centroid = bray_curtis_tibble %>%
    group_by(city_state) %>%
    summarize(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2)) 
  bray_curtis_tibble %>%
    left_join(dataframe %>%
                group_by(month_label, city_state) %>%
                summarise_at(vars(normalized_count), alpha_diversities) %>%
                dplyr::select(month_label, city_state, shannon)) %>%
    ggplot(aes(x = pcoa1, y= pcoa2)) +
    geom_point(shape = 21, aes(fill = as.factor(city_state), size = shannon)) +
    geom_point(data = bray_curtis_tibble_centroid, aes(x = pcoa1, y = pcoa2, fill = city_state), 
               size = 15, color = "black", shape = 22, stroke = 2, alpha = 0.75) +
    scale_size_continuous(range = c(2,10), transform = "reverse") +
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2') +
    xlab(paste0("PCo1: (",round(percentVar[1], 2),"%)")) +
    ylab(paste0("PCo2: (",round(percentVar[2], 2),"%)")) +
    theme_bw(30) +
    guides(fill = guide_legend(override.aes = list(shape=21, size = 7))) +
    theme(axis.text.x=element_text(angle=0, vjust=0.3),
          axis.text.y=element_text(),
          plot.title=element_text(hjust = 0.5)) +
    labs(fill = "City", size = "Shannon") +
    stat_ellipse(aes(color = as.factor(city_state)), alpha = 1,show.legend = F)
}


# Function to generate PCoA plot based on Aitchison distance
# Input: A dataframe with id, fileName (sample identifier), normalized_count and city_state as column
# Output: ggplot object PCoA plot using first two PcoA axes
# Pseudocount is added which is 0.65% of detection limit
aitchison_plot = function(dataframe)
{
  dataframe_wider = dataframe %>%
    pivot_wider(id_cols = id, names_from = fileName, values_from = normalized_count, values_fill = 0) %>% 
    column_to_rownames(var = 'id')
  pseudocount = min(dataframe_wider[dataframe_wider > 0]) * 0.0065
  aitchison_dist <- vegdist(t(dataframe_wider), method = "aitchison", pseudocount = pseudocount)
  aitchison_pcoa <- pco(aitchison_dist)
  eigenvalues = aitchison_pcoa$values
  percentVar = (eigenvalues/ sum(eigenvalues)) * 100
  aitchison_tibble = tibble(pcoa1 = aitchison_pcoa$vectors[,1],
                            pcoa2 = aitchison_pcoa$vectors[,2]) %>%
    mutate(fileName = rownames(aitchison_pcoa$vectors)) %>%
    left_join(dataframe) 
  aitchison_tibble_centroid = aitchison_tibble %>%
    group_by(city_state) %>%
    summarize(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2)) 
  aitchison_tibble %>%
    left_join(dataframe %>%
                group_by(month_label, city_state) %>%
                summarise_at(vars(normalized_count), alpha_diversities) %>%
                dplyr::select(month_label, city_state, shannon)) %>%
    ggplot(aes(x = pcoa1, y= pcoa2)) +
    geom_point(shape = 21, aes(fill = as.factor(city_state), size = shannon)) +
    geom_point(data = aitchison_tibble_centroid, aes(x = pcoa1, y = pcoa2, fill = city_state), 
               size = 15, color = "black", shape = 22, stroke = 2, alpha = 0.75) +
    scale_size_continuous(range = c(2,10), transform = "reverse") +
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2') +
    xlab(paste0("PCo1: (",round(percentVar[1], 2),"%)")) +
    ylab(paste0("PCo2: (",round(percentVar[2], 2),"%)")) +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    theme(axis.text.x=element_text(angle=0, vjust=0.3),
          axis.text.y=element_text(),
          plot.title=element_text(hjust = 0.5)) +
    labs(fill = "City", size = "Shannon diversity", color = "City") +
    stat_ellipse(aes(color = as.factor(city_state))) +
    theme_bw(30)
}




# Reading data from the saved dataframes
# These dataframes are combined table of output from various tools

# inspect file of k2_nt database downloaded on 2nd May 2023
k2_nt = read_delim(file = "data/k2_nt_20230502.inspect.species.krona.phylotree.txt", delim = "\t")

# Meta file for the data describing sampling city and collection month
meta = read_delim(file = "data/meta.tsv", delim = "\t", col_names = T) %>% 
  clean_names() %>%
  arrange(year, month) %>%
  mutate(month_label = factor(month_label, levels = unique(month_label))) %>%
  mutate(city_state = factor(city_state, levels = c("Kolkata", "Delhi", "Chennai", "Mumbai")))

# Bracken estimated read at the phyla level
# https://github.com/jenniferlu717/Bracken
# https://benlangmead.github.io/aws-indexes/k2
bracken_phylum = read_delim(file = "data/bracken_phylum.tsv", delim = "\t", col_names = T)
# Bracken estimated read at the genera level
bracken_genus = read_delim(file = "data/bracken_genus.tsv", delim = "\t", col_names = T)
# Bracken estimated read at the species level
bracken_species = read_delim(file = "data/bracken_species.tsv", delim = "\t", col_names = T)
# File with total read count (kraken) and domain level read count (bracken)
# https://ccb.jhu.edu/software/kraken/
bacterial_read_df = read_delim(file = "data/bacterial_read_df.tsv", delim = "\t", col_names = T)
# RGI output file
rgi_bowtie = read_delim(file = "data/rgi_bowtie.tsv", delim = "\t", col_names = T)

# CheckM2 quality report for all the MAGs
checkm2 = read_delim("data/checkm2_quality_report.csv", delim = "\t")

# geNomad output containing plasmid score for high quality MAGs
# https://github.com/apcamargo/genomad
genomad_plasmid = read_delim(file = paste0("data/genomad_plasmid_HQ_MAG.tsv"), delim = "\t", col_names = T)
# BAT output for high quality MAGs classification
# https://github.com/MGXlab/CAT_pack
bat_drep = read_delim(file = paste0("data/bat_HQ_MAG.tsv"), delim = "\t", col_names = T)
# BacAnt output for hiqh quality MAGs identifying AMR signatures
# https://github.com/xthua/bacant
amr_drep = read_delim(file = paste0("data/amr_HQ_MAG.tsv"), delim = "\t", col_names = T)
# BacAnt output for hiqh quality MAGs identifying transposon signatures
transposon_drep = read_delim(file = paste0("data/transposon_HQ_MAG.tsv"), delim = "\t", col_names = T)

# coverage of HQ MAG for each samples
coverage = read_delim(file = "data/sample_to_dereplicated_bin_coverage.tsv", delim = "\t", col_names = T)


# Read level analysis

# Merging multiple location data from the same city for the read level analysis

meta_select = meta %>% 
  dplyr::select(city_state, month, year, month_label) %>%
  distinct() %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bacterial_read_df_select = bacterial_read_df %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label) %>%
  summarise(total_reads = sum(total_reads), bacterial_reads = sum(bacterial_reads), total_mapped_reads = sum(total_mapped_reads), eukaryota_reads = sum(eukaryota_reads)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bracken_phylum_select = bracken_phylum  %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, name) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bracken_genus_select = bracken_genus %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, name) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bracken_species_select = bracken_species %>% 
  #mutate(across(c(new_est_reads), ~ifelse(. == 0, (.65*species_detection_limit), .))) %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, name) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
rgi_bowtie_selected = rgi_bowtie %>% 
  filter(`Reference Allele(s) Identity to CARD Reference Protein (%)` >= 90) %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, `ARO Term`, `Drug Class`) %>%
  summarise(`All Mapped Reads` = sum(`All Mapped Reads`)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))


bacterial_read_df %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  dplyr::select(city_state, total_reads, total_mapped_reads, bacterial_reads) %>%
  dplyr::rename(c(`Total reads` = total_reads, `Bacterial reads` = bacterial_reads, `Mapped reads` = total_mapped_reads)) %>%
  melt() %>%
  ggplot(aes(x = variable, y = value, fill = city_state)) +
  geom_point(aes(fill = city_state), size = 2.5, position = position_jitterdodge(), shape = 21, alpha = 0.7, stroke = 1) +
  geom_boxplot(alpha = 0.5, outliers = FALSE, linewidth = 1, median.linewidth = 1, show.legend = F) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  scale_size_manual(values = c(2.5, 10)) +
  labs(title="") +
  labs(x = "", y = "Read counts", fill = "") +
  theme_classic(40) +
  theme(axis.text.x=element_text(angle=0, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5),
        legend.position = "top") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 21, size = 7)))


# Bacterial phylums

bracken_phylum_select_bacterial = bracken_phylum_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(phylum) %>% distinct(), by = c("name" = "phylum")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(normalized_count = (new_est_reads/ bacterial_reads) * 1e06) %>%
  dplyr::rename(id = name)

top_bacterial_phylum = bar_plot_top10(bracken_phylum_select_bacterial)

plot_grid(plotlist = top_bacterial_phylum[[2]]) %>%
  annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                          bottom = text_grob("Collection month", color = "black",
                                             hjust = 0.5, x = 0.5, , size = 30, y = 1.5),
                          left = text_grob("Fraction of bacterial reads", color = "black", rot = 90, size = 30, x = 1.5)) %>%
  plot_grid(top_bacterial_phylum[[1]],  rel_widths = c(5, 1.5))


# Bacterial Genus

bracken_genus_select_bacterial = bracken_genus_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(genus) %>% distinct(), by = c("name" = "genus")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(normalized_count = (new_est_reads/ bacterial_reads) * 1e06) %>%
  dplyr::rename(id = name)

alpha_diversity_plot(bracken_genus_select_bacterial)

bray_curtis_plot(bracken_genus_select_bacterial)

aitchison_plot(bracken_genus_select_bacterial)

top_bacterial_genus = bar_plot_top10(bracken_genus_select_bacterial)

plot_grid(plotlist = top_bacterial_genus[[2]]) %>%
  annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                  bottom = text_grob("Collection month", color = "black",
                                     hjust = 0.5, x = 0.5, , size = 30, y = 1.5),
                  left = text_grob("Fraction of bacterial reads", color = "black", rot = 90, size = 30, x = 1.5)) %>%
  plot_grid(top_bacterial_genus[[1]],  rel_widths = c(5, 1.5))


# Bacterial species

bracken_species_select_bacterial = bracken_species_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(species) %>% distinct(), by = c("name" = "species")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(normalized_count = (new_est_reads/ bacterial_reads) * 1e06) %>%
  dplyr::rename(id = name)

alpha_diversity_plot(bracken_species_select_bacterial)

bray_curtis_plot(bracken_species_select_bacterial)

aitchison_plot(bracken_species_select_bacterial)

top_bacterial_species = bar_plot_top10(bracken_species_select_bacterial)

plot_grid(plotlist = top_bacterial_species[[2]]) %>%
  annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                  bottom = text_grob("Collection month", color = "black",
                                     hjust = 0.5, x = 0.5, , size = 30, y = 1.5),
                  left = text_grob("Fraction of bacterial reads", color = "black", rot = 90, size = 30, x = 1.5)) %>%
  plot_grid(top_bacterial_species[[1]],  rel_widths = c(5, 2))


# Bacterial species differential

bracken_species_bacterial_diff = bracken_species %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(species) %>% distinct(), by = c("name" = "species")) %>%
  left_join(meta, by = c("fileName" = "sample_code"))

bracken_species_bacterial_meta = bracken_species_bacterial_diff %>%
  dplyr::select(fileName, city_state, month_label) %>%
  distinct()

bracken_species_bacterial_df = bracken_species_bacterial_diff %>%
  tidyr::pivot_wider(id_cols = name, names_from = fileName, values_from = new_est_reads, values_fill = 0) %>%
  column_to_rownames('name') 

differential_up = differential_dn = list()
differential_df = tibble(month_label = character(), numerator = character(), denominator = character(), enriched = character(), depleted = character())
log2Foldchange_cutoff = 1

city_vector = levels(as.factor(bracken_species_bacterial_meta$city_state))

for (month in levels(bracken_species_bacterial_meta$month_label))
{
  for (i in c(1:(length(city_vector)-1)))
  {
    denominator = city_vector[i]
    for (j in c((i+1):length(city_vector)))
    {
      numerator = city_vector[j]
      if (!(numerator == denominator))
      {
        select_meta = bracken_species_bacterial_meta %>%
          filter(month_label == month) %>%
          filter(city_state == numerator) %>%
          bind_rows(bracken_species_bacterial_meta %>%
                      filter(month_label == month) %>%
                      filter(city_state == denominator))  %>%
          droplevels()
        if((select_meta %>% pull(city_state) %>% unique() %>% length() > 1) &
           (min(table(select_meta$city_state)) > 1))
        {
          selected_df = bracken_species_bacterial_df[, select_meta$fileName] %>% filter(rowSums(.) > 0)
          aldex_x <- aldex.clr(selected_df, as.character(select_meta$city_state), denom="all", useMC = T)
          ttest.test <- aldex.ttest(aldex_x, paired.test = F) %>%
            bind_cols(aldex.effect(aldex_x, paired.test = F))
          enriched = rownames(ttest.test %>%  filter(we.eBH < 0.05) %>% filter(diff.btw > log2Foldchange_cutoff))
          depleted = rownames(ttest.test %>%  filter(we.eBH < 0.05) %>% filter(diff.btw < -log2Foldchange_cutoff))
          differential_up[[paste0(numerator,"_",denominator, "_", month)]] = enriched
          differential_dn[[paste0(numerator,"_",denominator, "_", month)]] = depleted
          differential_df = differential_df %>%
            bind_rows(c(month_label = month, numerator = numerator, denominator = denominator,
                        enriched = length(enriched), depleted = length(depleted)))
        }
      }
    }
  }
}

differential_df_other = differential_df 
colnames(differential_df_other) = c("month_label", "denominator", "numerator", "depleted", "enriched")
differential_df_other = differential_df_other[,c("month_label", "numerator", "denominator", "enriched", "depleted")]

differential_df %>% 
  bind_rows(differential_df_other) %>%
  mutate(month_label = factor(month_label, levels = levels(meta$month_label))) %>%
  dplyr::select(month_label, numerator, denominator, enriched) %>%
  melt(id.vars = c("month_label", "numerator", "denominator")) %>%
  filter(value > 0) %>%
  mutate(denominator = factor(denominator, levels = levels(meta$city_state))) %>%
  mutate(numerator = factor(numerator, levels = levels(meta$city_state))) %>%
  {ggplot(.,aes(x = month_label, y = numerator)) +
      geom_point(aes(size = as.numeric(value)), shape = 21, alpha = 1, fill = "red4") +
      scale_size_continuous(range = c(1,15)) +
      ylab(label = "") +
      xlab(label = "Collection month") +
      labs(size = "Count") +
      facet_wrap(~denominator, scales = "free_y", ncol = 1) +
      theme_bw(25) +
      scale_y_discrete(limits = rev) +
      theme(axis.text.x=element_text(angle=90, vjust=0.3, size = 30),
            axis.text.y=element_text(size = 35),
            axis.title.x = element_text(size = 35),
            plot.title=element_text(hjust = 0.5),
            strip.background = element_rect(fill="white"),
            strip.text = element_text(size = 30),
            legend.position = "right",
            legend.title = element_text(size = 35),
            legend.text = element_text(size = 35))}

# Antimicrobial resistant genes (ARG)

# trimmed rgi, keeping only those aro for which we have Average MAPQ >10 and completely mapped reads more than 10
rgi_bowtie_trimmed = rgi_bowtie %>%
  filter(`Reference Allele(s) Identity to CARD Reference Protein (%)` >= 90) %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, `ARO Term`, `Drug Class`) %>%
  summarise(`All Mapped Reads` = sum(`All Mapped Reads`), .groups = 'drop') %>%
  clean_names() %>%
  mutate(fileName = paste0(city_state, "_", month_label)) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(normalized_count = (all_mapped_reads/total_reads) * 1e06) %>%
  dplyr::rename(id = aro_term) %>%
  filter(normalized_count > 10)

alpha_diversity_plot(rgi_bowtie_trimmed)

bray_curtis_plot(rgi_bowtie_trimmed)

aitchison_plot(rgi_bowtie_trimmed)

top_arg = bar_plot_top10(rgi_bowtie_trimmed)

plot_grid(plotlist = top_arg[[2]]) %>%
  annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                  bottom = text_grob("Collection month", color = "black",
                                     hjust = 0.5, x = 0.5, , size = 30, y = 1.5),
                  left = text_grob("Fraction of total reads", color = "black", rot = 90, size = 30, x = 1.5)) %>%
  plot_grid(top_arg[[1]],  rel_widths = c(5, 2))

# Defining betalactam and MLS drug classes
# https://www.nature.com/articles/s41467-022-29283-8
betalactam = c("penam", "cephalosporin", "carbapenem", "cephamycin", "penem", "monobactam")
MLS = c("macrolide", "lincosamide", "streptogramin", "streptogramin A", "streptogramin B")


# Defining 'Multi-drug' drug classes
rgi_bowtie_trimmed_drug_class = rgi_bowtie_trimmed %>%
  mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
  mutate(drug_class_group = drug_class,
         drug_class_multidrug = drug_class) 

for (i in c(1:nrow(rgi_bowtie_trimmed_drug_class)))
{
  drug_class_list = strsplit(rgi_bowtie_trimmed_drug_class$drug_class[i], split = "; ")[[1]]
  drug_class_list[drug_class_list %in% betalactam] <- "Beta-lactam"
  drug_class_list[drug_class_list %in% MLS] <- "MLS"
  drug_class_list = unique(drug_class_list)
  if(length(drug_class_list) >= 2)
  {rgi_bowtie_trimmed_drug_class$drug_class_multidrug[i] = "multidrug"}
  else
  {rgi_bowtie_trimmed_drug_class$drug_class_multidrug[i] = drug_class_list}
  rgi_bowtie_trimmed_drug_class$drug_class_group[i] = paste(drug_class_list, collapse = ";")
}

top_drug_class = bar_plot_top10(
  rgi_bowtie_trimmed_drug_class %>%
    group_by(drug_class_multidrug, city_state) %>%
    mutate(normalized_count = sum(normalized_count)) %>% ungroup() %>%
    dplyr::rename(arg = id, id = drug_class_multidrug))

plot_grid(plotlist = top_drug_class[[2]]) %>%
  annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                  bottom = text_grob("Collection month", color = "black",
                                     hjust = 0.5, x = 0.5, , size = 30, y = 1.5),
                  left = text_grob("Fraction of total reads", color = "black", rot = 90, size = 30, x = 1.5)) %>%
  plot_grid(top_drug_class[[1]],  rel_widths = c(5, 1.7))



multidrug_upset_list = rgi_bowtie_trimmed_drug_class %>% filter(drug_class_multidrug %in% "multidrug") %>%
  dplyr::select(city_state, drug_class_group) %>%
  mutate(drug_class_group = str_replace_all(drug_class_group, pattern = ";", replacement = "&")) %>%
  group_by(city_state) %>%
  group_map(
    .f = ~ .x %>% 
      group_by(drug_class_group) %>%
      summarise(count = n()) %>%
      filter(count >= 10) %>%
      deframe() %>%
      fromExpression() %>%
      upset(matrix.dot.alpha = 0.1, 
                    nintersects = NA, 
                    nsets = 200, 
                    order.by = "freq", 
                    decreasing = T, 
                    mb.ratio = c(0.4, 0.4),
                    number.angles = 0, 
                    text.scale = 2, 
                    point.size = 2, 
                    line.size = 1
      ) %>%
      as.grob()
  )
plot_grid(plotlist = (multidrug_upset_list), labels = levels(rgi_bowtie_trimmed_drug_class$city_state))

# ARG Differential analysis

rgi_diff = rgi_bowtie %>% 
  filter(`Reference Allele(s) Identity to CARD Reference Protein (%)` >= 90) %>%
  left_join(meta, by = c("fileName" = "sample_code"))

rgi_meta = rgi_diff %>%
  dplyr::select(fileName, city_state, month_label) %>%
  distinct()

rgi_diff_df = rgi_diff %>%
  tidyr::pivot_wider(id_cols = `ARO Term`, names_from = fileName, values_from = `All Mapped Reads`, values_fill = 0) %>%
  column_to_rownames('ARO Term') 


differential_rgi_up = differential_rgi_dn = list()
differential_rgi_df = tibble(month_label = character(), numerator = character(), denominator = character(), enriched = character(), depleted = character())

city_vector = levels(as.factor(rgi_meta$city_state))

for (month in levels(rgi_meta$month_label))
{
  for (i in c(1:(length(city_vector)-1)))
  {
    denominator = city_vector[i]
    #j = 4
    for (j in c((i+1):length(city_vector)))
    {
      numerator = city_vector[j]
      if (!(numerator == denominator))
      {
        select_meta = rgi_meta %>%
          filter(month_label == month) %>%
          filter(city_state == numerator) %>%
          bind_rows(rgi_meta %>%
                      filter(month_label == month) %>%
                      filter(city_state == denominator)) %>%
          droplevels()
        if((select_meta %>% pull(city_state) %>% unique() %>% length() > 1) &
           (min(table(select_meta$city_state)) > 1))
        {
          selected_df = rgi_diff_df[, select_meta$fileName] %>% filter(rowSums(.) > 0)
          aldex_x <- aldex.clr(selected_df, as.character(select_meta$city_state), denom="all", useMC = T)
          ttest.test <- aldex.ttest(aldex_x, paired.test = F) %>%
            bind_cols(aldex.effect(aldex_x, paired.test = F))
          enriched = rownames(ttest.test %>%  filter(we.eBH < 0.05) %>% filter(diff.btw > log2Foldchange_cutoff))
          depleted = rownames(ttest.test %>%  filter(we.eBH < 0.05) %>% filter(diff.btw < -log2Foldchange_cutoff))
          differential_rgi_up[[paste0(numerator,"_",denominator, "_", month)]] = enriched
          differential_rgi_dn[[paste0(numerator,"_",denominator, "_", month)]] = depleted
          differential_rgi_df = differential_rgi_df %>%
            bind_rows(c(month_label = month, numerator = numerator, denominator = denominator,
                        enriched = length(enriched), depleted = length(depleted)))
        }
      }
    }
  }
}

differential_rgi_df_other = differential_rgi_df 
colnames(differential_rgi_df_other) = c("month_label", "denominator", "numerator", "depleted", "enriched")
differential_rgi_df_other = differential_rgi_df_other[,c("month_label", "numerator", "denominator", "enriched", "depleted")]

differential_rgi_df %>% 
  bind_rows(differential_rgi_df_other) %>%
  dplyr::select(month_label, numerator, denominator, enriched) %>%
  melt(id.vars = c("month_label", "numerator", "denominator")) %>%
  tibble() %>%
  mutate(month_label = factor(month_label, levels = levels(meta$month_label))) %>%
  mutate(denominator = factor(denominator, levels = levels(meta$city_state))) %>%
  mutate(numerator = factor(numerator, levels = levels(meta$city_state))) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  {ggplot(.,aes(x = month_label, y = numerator)) +
      geom_point(aes(size = as.numeric(value)), shape = 21, alpha = 1, fill = "red4") +
      scale_size_continuous(range = c(2,15)) +
      ylab(label = "") +
      xlab(label = "Collection month") +
      scale_x_discrete(drop=FALSE) +
      scale_y_discrete(limits = rev) +
      labs(size = "Number of differenitally enriched ARGs") +
      facet_wrap(~denominator, scales = "free_y", nrow = 2) +
      theme_bw(30) +
      theme(axis.text.x=element_text(angle=90, vjust=0.3, size = 20),
            axis.text.y=element_text(size = 35),
            plot.title=element_text(hjust = 0.5),
            strip.background = element_rect(fill="white"),
            strip.text = element_text(size = 50),
            legend.position = "top",
            legend.title = element_text(size = 35),
            legend.text = element_text(size = 35)) +
      guides(size = guide_legend(title.position = "top", title.hjust=0.5))}


# MAG level analysis

# Reading checkm2 quality score of MAGs and retaining High quality MAGs
checkm2 = checkm2 %>%
  separate(Name, into = c("sample_code", "extra"), sep = "_", remove = F) %>%
  left_join(meta) %>% 
  #filter(Completeness > 50, Contamination < 10) # Medium Quality
  filter(Completeness > 90, Contamination < 5) # High Quality

# MAGs taxonomic classification
bat_drep_resolved = bat_drep %>%
  mutate(phylum = case_when(phylum ==  "no support" ~ "Unresolved MAG",
                            TRUE ~ "Resolved MAG"),
         class = case_when(class ==  "no support" ~ "Unresolved MAG",
                           TRUE ~ "Resolved MAG"),
         order = case_when(order ==  "no support" ~ "Unresolved MAG",
                           TRUE ~ "Resolved MAG"),
         family = case_when(family ==  "no support" ~ "Unresolved MAG",
                            TRUE ~ "Resolved MAG"),
         genus = case_when(genus ==  "no support" ~ "Unresolved MAG",
                           TRUE ~ "Resolved MAG"),
         species = case_when(species ==  "no support" ~ "Unresolved MAG",
                             TRUE ~ "Resolved MAG")) %>%
  left_join(checkm2) %>%
  dplyr::select(month_label, city_state, phylum, class, order, family, genus, species) %>%
  melt(id.vars = c("month_label", "city_state"))

# Plotting percentage of MAGs resolved at each taxonomical level
bat_drep_resolved %>%  
  group_by(city_state, variable, value) %>% summarise(count = n(), .groups = "drop") %>%
  left_join(bat_drep_resolved %>%  
              group_by(city_state, variable) %>% summarise(total = n(), .groups = "drop")) %>%
  mutate(percentage = round((count*100)/total, 2)) %>%
  mutate(value = fct_relevel(value, c("Unresolved MAG", "Resolved MAG"))) %>%
  filter(value == "Resolved MAG") %>%
  ggplot(aes(y = variable, x = percentage)) +
  geom_bar(aes(fill = city_state), color = "black", position = "dodge", stat = "identity", alpha = 1, linewidth = 1) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic(45) +
  scale_y_discrete(limits=rev) +
  guides(fill = guide_legend(override.aes = list(shape=21)),
         color = guide_legend(override.aes = list(shape=21))) +
  theme(axis.text.x=element_text(angle=90, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5)) +
  labs(color = "City", fill = "City") +
  ylab(paste0("")) +
  xlab(paste0("Percentage of MAGs resolved by BAT")) +
  ggtitle("")

# Comparing stats of MAG that were resolved with unresolved MAGs at the species level
bat_drep %>%
  left_join(checkm2) %>%
  mutate(Log_contig_N50 = log(Contig_N50)) %>%
  mutate(Genome_Size = Genome_Size/1000000) %>%
  dplyr::select(species, Completeness, Contamination, Log_contig_N50, Genome_Size, Total_Contigs) %>%
  dplyr::rename(c(`log(Contigs N50)` = Log_contig_N50), `Genome Size (Mb)` = Genome_Size, `Total contigs` = Total_Contigs) %>%
  mutate(species = case_when(species ==  "no support" ~ "Novel",
                             TRUE ~ "Known")) %>%
  melt(id.vars = c("species")) %>%
  tidyplot(x = species, y = value, color = species, alpha = 1) %>%
  add_boxplot(dodge_width = 0.2, box_width = 0.9, outlier.size = 0.2, outlier.alpha = 0.2, alpha = 0.5) %>%
  add_test_pvalue(ref.group = 1, p.adjust.method = "BH", hide_info = F, label.size = 5,bracket.nudge.y = .2) %>%
  theme_tidyplot(fontsize = 15) %>%
  remove_x_axis_title() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_title() %>%
  adjust_legend_title("MAGs at the\nspecies taxa") %>%
  split_plot(by = variable)

'
bat_drep %>%
left_join(checkm2) %>%
mutate(Log_contig_N50 = log(Contig_N50)) %>%
mutate(Genome_Size = Genome_Size/1000000) %>%
dplyr::select(species, Completeness, Contamination, Log_contig_N50, Genome_Size, Total_Contigs) %>%
dplyr::rename(c(`log(Contigs N50)` = Log_contig_N50), `Genome Size (Mb)` = Genome_Size, `Total contigs` = Total_Contigs) %>%
mutate(species = case_when(species ==  "no support" ~ "Novel",
TRUE ~ "Known")) %>%
melt(id.vars = c("species")) %>%
group_by(variable) %>%
summarise(p.value = formatC(t.test(value ~ species)$p.value, format = "e", digits = 2)) %>%
ungroup()
'
  
# Sankey plot for bin classification
bat_drep_sankey = bat_drep %>% left_join(checkm2)

bat_drep_df = bat_drep_sankey %>%
  separate(superkingdom, into = c("superkingdom", "superkingdom_score"), sep = ":") %>%
  separate(phylum, into = c("phylum", "phylum_score"), sep = ":") %>%
  separate(class, into = c("class", "class_score"), sep = ":") %>%
  separate(order, into = c("order", "order_score"), sep = ":") %>%
  separate(family, into = c("family", "family_score"), sep = ":") %>%
  separate(genus, into = c("genus", "genus_score"), sep = ":") %>%
  separate(species, into = c("species", "species_score"), sep = ":") %>%
  dplyr::select(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "city_state", "Name"))
bat_drep_df_contig_count = bat_drep_df %>% group_by(city_state) %>% summarise(contig_count = n())
prop_cutoff = 0.02

df = bat_drep_df %>%
  group_by(city_state) %>%
  group_modify(~
                 .x %>%
                 ggsankey::make_long(superkingdom, phylum,class,order, family, genus,species) %>% 
                 filter(!((node == "no support") )) %>%
                 filter(!((is.na(node)) | (is.na(next_node)))) %>%
                 group_by(node, next_node) %>% mutate(value = n()) %>% 
                 ungroup() %>% distinct() %>% 
                 mutate(prop = value/(bat_drep_df_contig_count %>% filter(city_state == .y$city_state) %>% pull(contig_count))) %>%
                 filter(prop >= prop_cutoff) %>%
                 mutate(city = .y$city_state)
  ) %>%
  ungroup()

df %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node, value = prop)) +
  geom_sankey(flow.alpha = 1,
              node.color = "gray30") +
  geom_sankey_label(size = 4.1, color = "gray30", fill = "white") +
  theme_sankey(base_size = 25) +
  labs(x = NULL) +
  theme(legend.position = "NULL",
        plot.title = NULL) +
  facet_wrap(~city)

# Mobile genetic elements analysis (Plasmid + Transoposon)
genomad_plasmid = genomad_plasmid %>% filter(plasmid_score >= 0.9) %>% 
  left_join(checkm2) %>% 
  mutate(bin_contig = paste0(Name, "_", seq_name))  

amr_drep_checkm2 = amr_drep %>% mutate(name_contig = paste0(Name, "_", X2)) %>% left_join(checkm2, by = c("Name")) %>% filter(X11 >= 90)
transposon_drep_checkm2 = transposon_drep %>% mutate(name_contig = paste0(Name, "_", X1)) %>% left_join(checkm2, by = c("Name")) %>% filter(X3 >= 90)

checkm2 %>%
  dplyr::select(c("city_state", "Name", "Total_Contigs")) %>%
  left_join(amr_drep_checkm2 %>% dplyr::select(X2, Name) %>% distinct() %>% 
              group_by(Name) %>% summarise(AMR_contigs = n())) %>% #Number of unique contigs per bin with AMR
  replace(is.na(.), 0) %>%
  mutate('Proportion of ARG contigs in MAGs' = (AMR_contigs/Total_Contigs)) %>%
  tidyplot(x = city_state, y = `Proportion of ARG contigs in MAGs`) %>%
  add_boxplot(dodge_width = 0.2, box_width = 0.9, show_outliers = F , alpha = 0.7, fill = "white", linewidth = 1) %>%
  add_data_points_beeswarm(shape = 21, alpha = 0.5, jitter_width = 0.6, size = 3.5, ) %>%
  add_test_pvalue(p.adjust.method = "BH", hide_info = T, label.size = 15, label = "p.adj",hide.ns = T, bracket.nudge.y = -0.725) %>%
  theme_tidyplot(fontsize = 40) %>%
  adjust_y_axis(limits = c(0, 0.1)) %>%
  remove_x_axis_title() %>%
  adjust_legend_title("") %>%
  adjust_size(350,250)

# Making df with unique contigs for each drug class, city and bin
amr_drep_checkm2_topPR_df = amr_drep_checkm2 %>% 
  mutate(X14 = replace_na(X14, "UNKNOWN")) %>%
  dplyr::select(X14, name_contig, city_state, Name) %>% 
  distinct()

# For each city top 10 drug classes with most number of contigs
top_pr = amr_drep_checkm2_topPR_df %>% 
  filter(!(is.na(X14))) %>%
  group_by(city_state, X14) %>%
  summarize(contig_count = n(), .groups = "drop") %>%
  group_by(city_state) %>% 
  top_n(10,wt = contig_count) %>% ungroup()
top_pr_order = top_pr %>%
  group_by(X14) %>% summarise(contig_count = sum(contig_count), .groups = "drop") %>% 
  arrange(desc(contig_count)) %>%
  mutate(X14 = factor(X14, levels = unique(X14))) %>% 
  mutate(X14 = forcats::fct_relevel(X14, "UNKNOWN", after = Inf))
top_pr_color = setNames(paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(as.factor(top_pr_order$X14)))],levels(as.factor(top_pr_order$X14)))

# Total AMR carrying contigs per city
city_contig_count = amr_drep_checkm2 %>% dplyr::select(name_contig, city_state, Name) %>%
  distinct() %>%
  group_by(city_state) %>%
  summarise(city_contig_count = n(), .groups = "drop")

# Top 10 Drug class for each city
amr_drep_checkm2_topPR_df %>% 
  filter(X14 %in% top_pr$X14) %>%
  group_by(X14, city_state) %>%
  summarise(contig_count_pr = n(), .groups = "drop") %>%
  left_join(city_contig_count) %>%
  mutate(contig_count_pr_prop = contig_count_pr/city_contig_count) %>%
  mutate(X14 = factor(X14, levels = levels(top_pr_order$X14))) %>%
  ggplot(aes(y = X14, x = contig_count_pr_prop)) +
  geom_bar( position = "identity", stat = "identity", aes(fill = X14), color = "black") +
  scale_fill_manual(values = top_pr_color) +
  scale_y_discrete(limits=rev) +
  labs(fill = "Drug class", 
       x = "Proportion of contigs carrying Drug class", y = "Top 10 Drug class for each city") +
  facet_wrap( ~ city_state, nrow = 1) +
  theme_classic(25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.8, "cm"),
        legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Remove duplicate AMR on same contigs
amr_drep_checkm2_topamr_df = amr_drep_checkm2 %>%
  mutate(X14 = replace_na(X14, "UNKNOWN")) %>%
  dplyr::select(X6, name_contig, city_state, Name, X14) %>% 
  distinct()

# For each city top 10 AMR with most number of contigs
top_amr = amr_drep_checkm2_topamr_df %>% 
  group_by(city_state, X6, X14) %>%
  summarize(contig_count = n(), .groups = "drop") %>%
  group_by(city_state) %>%
  top_n(10,wt = contig_count) %>%
  ungroup()
# AMR ordering first with Drug class with most contigs, followed by AMR with most contigs
top_amr_order = top_amr %>%
  group_by(X6) %>% mutate(contig_count_per_amr = sum(contig_count)) %>% ungroup() %>%
  group_by(X14) %>% mutate(contig_count_per_dc = n()) %>% ungroup() %>%
  arrange(desc(contig_count_per_dc), X14, desc(contig_count_per_amr)) %>% ungroup() %>%
  mutate(X14 = replace_na(X14, "UNKNOWN")) %>%
  mutate(X14 = factor(X14, levels = unique(X14)), 
         X6 = factor(X6, levels = unique(X6)))
top_amr_color = setNames(paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(top_amr_order$X14))],levels(top_amr_order$X14))

# Top 10 AMR per city colored by drug class
amr_drep_checkm2_topamr_df %>% 
  filter(X6 %in% top_amr$X6) %>%
  group_by(X6, city_state) %>%
  mutate(contig_count_amr = n()) %>% ungroup() %>%
  left_join(city_contig_count) %>%
  mutate(contig_count_amr_prop = contig_count_amr/city_contig_count) %>%
  mutate(X14 = replace_na(X14, "UNKNOWN")) %>%
  mutate(X6 = factor(X6, levels = levels(top_amr_order$X6)),
         X14 = factor(X14, levels = levels(top_pr_order$X14))) %>%
  ggplot(aes(y = X6, x = contig_count_amr_prop)) +
  geom_bar( position = "identity", stat = "identity", aes(fill = X14), color = "black") +
  scale_fill_manual(values = top_pr_color) +
  scale_y_discrete(limits=rev) +
  labs(fill = "Drug class", 
       x = "Proportion of contigs carrying ARG", y = "Top 10 ARGs for each city") +
  facet_wrap( ~ city_state, nrow = 1) +
  theme_classic(25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "cm"),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Unique Contig count per drug class per city
drug_class_contig_count = amr_drep_checkm2 %>% 
  mutate(X14 = replace_na(X14, "UNKNOWN")) %>%
  dplyr::select(X14, name_contig, city_state, Name) %>%
  distinct() %>% #One drug class instance per contig.
  group_by(X14, city_state) %>% 
  summarise(dc_contig_count = n(),.groups = "drop")

# Unique Contig count per drug class "CONTAINING MGE" per city
drug_class_mge_contig_count = amr_drep_checkm2 %>% 
  dplyr::select(X14, name_contig, city_state, Name) %>%
  distinct() %>%
  mutate(plasmid_transposon = case_when((name_contig %in% transposon_drep_checkm2$name_contig) &  
                                          (name_contig %in% genomad_plasmid$bin_contig) ~ "Yes",
                                        TRUE ~ "No")) %>%
  filter(plasmid_transposon == "Yes") %>%
  group_by(X14, city_state) %>%
  summarise(mge_dc_contig_count = n(), .groups = "drop")

# Unique AMR count per drug class "CONTAINING MGE" per city
drug_class_amr_count = amr_drep_checkm2 %>% 
  dplyr::select(X6, X14, city_state, Name, name_contig) %>%
  distinct() %>%
  group_by(X14, city_state) %>%
  summarise(dc_amr_count = n(), .groups = "drop")

# Plot with top 10 drug class for each city
# X-axis has proportion of `contigs with specific drug class` for each city
# Y-axis has proportion of `ARG carrying contigs` with MGE for each city
# Size of dot is number of unique ARG for the drug class
drug_class_contig_count %>% left_join(drug_class_mge_contig_count) %>% left_join(drug_class_amr_count) %>% left_join(city_contig_count) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  group_by(X14) %>% 
  mutate(contig_count_per_dc_overall = sum(dc_contig_count)) %>% ungroup() %>%
  mutate(`dc_contig_per_city_contig` = dc_contig_count/city_contig_count, #Proportion of contig for each drug class out of total contig per city
         mge_dc_contig_count_per_dc_contig_count  = mge_dc_contig_count/dc_contig_count) %>% #Proportion of contig with mge for each dc, for each city
  filter(X14 %in% top_pr$X14) %>%
  mutate(X14 = factor(X14, levels = levels(top_pr_order$X14))) %>%
  filter(dc_amr_count > 0) %>%
  ggplot(aes(x = dc_contig_per_city_contig, y = mge_dc_contig_count_per_dc_contig_count)) +
  geom_point(aes(fill = X14, size = dc_amr_count, shape = city_state), alpha = 1) +
  scale_size_continuous(range = c(2,15)) +
  scale_shape_manual(values = c(21,22,23,24)) +
  theme_bw(20) +
  labs(fill = "Drug class", size = "Unique ARG \nper Drug class", shape = "City", 
       x = "Proportion of ARG containing contigs", y = "Proportion of specific ARG contig with MGE") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)),
         shape = guide_legend(override.aes = list(fill = "white", size = 5)),
         size = guide_legend(override.aes = list(shape = 21, fill = "white"))) +
  scale_fill_manual(values = top_pr_color)

contig_amr_genomad_transposon = amr_drep_checkm2 %>% 
  dplyr::select(X6, X14, name_contig, city_state, Name) %>%
  mutate(plasmid_transposon = case_when((name_contig %in% transposon_drep_checkm2$name_contig) &  
                                          (name_contig %in% genomad_plasmid$bin_contig) ~ "Contigs with MGE",
                                        TRUE ~ "Contigs without MGE"))

# Plot with multiple arg for same drug class on a contig
contig_amr_genomad_transposon %>% 
  group_by(X14, name_contig, city_state, plasmid_transposon) %>% 
  summarise(pr_per_count = n(), .groups = "drop") %>%
  filter(!is.na(X14)) %>%
  filter(pr_per_count > 1) %>%
  mutate(X14 = factor(X14, levels = levels(top_pr_order$X14))) %>%
  ggplot(aes(x = city_state, y = pr_per_count)) +
  geom_jitter(aes(fill = X14, shape = plasmid_transposon) , alpha = 1, size = 10, width = 0.3) +
  theme_bw(35) +
  scale_shape_manual(values = c(23,21)) +
  labs(fill = "Drug class", 
       x = "", y = "Number of ARGs per contig", shape = "") +
  guides(fill = guide_legend(override.aes = list(size=10, shape = 21)),
         shape = guide_legend(override.aes = list(size = 5, fill = "grey70")),
         size = guide_legend(override.aes = list(shape = 21, fill = "grey70"))) +
  scale_fill_manual(values = top_pr_color)


# Network analysis for dereplicated high quality MAGs

# Calculating bin level coverage for each sample.
coverage_bin = coverage %>%
  left_join(meta, by = c("sampleName" = "sample_code")) %>%
  group_by(bin, city_state, month_label) %>%
  summarise(binLen = sum(contigLen), nucleotide_mapped = sum(nucleotide_mapped),.groups = "drop") %>%
  mutate(avg_nucleotide_mapped = nucleotide_mapped/binLen) %>%
  mutate(sampleName = paste0(city_state, "_", month_label))

coverage_bin_wide = coverage_bin %>%
  pivot_wider(names_from = sampleName, id_cols = bin, values_from = avg_nucleotide_mapped) %>%
  column_to_rownames(var = 'bin')

bin_with_mean_cov_0.5 = coverage_bin_wide %>% rowMeans() %>% stack() %>% filter(values > 0.5)
coverage_bin_wide_mean_cov_0.5 = coverage_bin_wide[bin_with_mean_cov_0.5$ind,]

bat_drep_network = bat_drep_df %>%
  filter(Name %in% rownames(coverage_bin_wide_mean_cov_0.5))
phylum_color = tibble(phylum = unique(bat_drep_network$phylum), 
                      color = paletteer::paletteer_d("ggsci::category20_d3")[1:length(unique(bat_drep_network$phylum))]) %>%
  left_join(bat_drep_network %>% group_by(phylum) %>% summarise(count = n())) %>% arrange(desc(count))
phylum_color$color[which(phylum_color$phylum == "no support")] = "grey50"
phylum_color$color[which(is.na(phylum_color$phylum))] = "white"
bat_drep_network = bat_drep_network %>%
  left_join(phylum_color) %>%
  relocate(Name, .before = superkingdom) %>%
  mutate(novel = ifelse(species == "no support", "Novel", "Known"))

phylum_network = bat_drep_network %>%
  dplyr::select(phylum, color, novel) %>%
  group_by(phylum, color, novel) %>%
  summarise(count = n()) %>%
  arrange((count)) %>%
  filter(!(is.na(phylum))) %>%
  mutate(phylum = as.factor(phylum))
phylum_network_color = setNames(phylum_network$color, phylum_network$phylum)
phylum_network %>%
  ggplot(aes(y = fct_inorder(phylum), x = novel, fill = phylum)) +
  geom_point(aes(shape = novel, size = count)) +
  theme_classic(25) +
  scale_size_continuous(range = c(4,12)) +
  scale_shape_manual(values = c(21,22)) +
  scale_fill_manual(values = phylum_network_color, name = "Phylum") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x=element_text(angle=90, size = 25), 
        axis.text.y=element_text(size = 25),
        legend.position = "bottom") +
  guides(fill="none", shape = "none") +
  coord_fixed(ratio = .5) +
  xlab("") +
  ylab("")

network_list = list()
clp_list = list()
for (city_selected in unique(meta$city_state))
{
  city_samples = meta_select %>% filter(city_state == city_selected) %>% pull(fileName)
  coverage_bin_wide_mean_cov_0.5_location = coverage_bin_wide_mean_cov_0.5 %>% 
    dplyr::select(city_samples)
  
  paired_pvalue = tibble(to = character(), from = character(), pvalue = numeric(), cor = numeric())
  for (i in rownames(coverage_bin_wide_mean_cov_0.5_location))
  {
    for (j in rownames(coverage_bin_wide_mean_cov_0.5_location))
    {
      if(!(i == j))
      {
        cor_test = cor.test(as.numeric(coverage_bin_wide_mean_cov_0.5_location[i,]), as.numeric(coverage_bin_wide_mean_cov_0.5_location[j,]))
        pval = cor_test$p.value
        est = as.numeric(cor_test$estimate)
        paired_pvalue = paired_pvalue %>%
          bind_rows(tibble(to = i, from = j, pvalue = pval, cor = est))
      }
    }
  }
  paired_pvalue = paired_pvalue %>%
    mutate(fdr = p.adjust(pvalue, method = "BH"))
  
  edge_list = paired_pvalue %>% filter(fdr < 0.01)
  nodes_list = bat_drep_network %>% filter(Name %in% unique(edge_list$to))
  #Network with given nodes and edge list
  network = graph_from_data_frame(vertices = nodes_list, d = edge_list, directed=F)
  network = simplify(network, remove.multiple = T, remove.loops = T, edge.attr.comb	="first")
  #computing optimal cluster
  clp = cluster_optimal(network)
  network_list[[city_selected]] = network
  clp_list[[city_selected]] = clp
  
  groups_loop = groups(clp) %>% stack()
  groups_table_loop = groups_loop %>% group_by(ind) %>% summarise(count = n()) %>% 
    mutate(color = paletteer_d("RColorBrewer::Set3")[1:nrow(.)]) %>%
    mutate(color = case_when(count < 10 ~ NA, TRUE ~ (color)))
  groups_loop = groups_loop %>% filter(ind %in% groups_table_loop$ind) %>%
    left_join(groups_table_loop) %>%
    droplevels()
  plot(clp, network, vertex.label=NA, vertex.size = 5, edge.width = 1.5,
       col = V(network)$color,
       mark.col = NA, 
       mark.border = "black", alpha = 0.1,
       vertex.shape = if_else(V(network)$species == "no support", "circle", "square"),
       edge.color = if_else(E(network)$cor < 0, "red4", "blue4"), 
       edge.width = ((abs(E(network)$cor))),
       edge.curved = T,
       vertex.frame.width = 2,
       layout= layout_with_fr(network, weights = abs(E(network)$cor)),
       frame = T, main = city_selected, rescale = T)
}


community_temporal_df = tibble()
groups_df = tibble()

for (city_selected in unique(meta_select$city_state))
{
  city_samples = meta_select %>% filter(city_state == city_selected) %>% pull(fileName)
  coverage_bin_wide_mean_cov_0.5_location = coverage_bin_wide_mean_cov_0.5 %>% 
    dplyr::select(city_samples)
  
  groups = igraph::groups(clp_list[[city_selected]]) %>% stack()
  
  groups_table = groups %>% group_by(ind) %>% summarise(count = n()) %>% 
    slice_max(order_by = count, n = 6, with_ties = F) #Selecting top 6 community with max members for each city
  
  groups = groups_table %>%
    left_join(groups) %>%
    droplevels() %>%
    left_join(bat_drep_network[,c('Name', 'phylum', 'novel', 'color')], by = c("values" = "Name"))
  
  groups_df = groups_df %>% rbind(groups %>% mutate(city = city_selected) %>% as.tibble())
  
  city_community_abundance = tibble()
  for(i in groups_table$ind)
  {
    members = groups %>% filter(ind == i) %>% pull(values)
    members_df = coverage_bin_wide_mean_cov_0.5_location[members,]
    city_community_abundance = city_community_abundance %>%
      bind_rows(colMeans(members_df))
  }
  
  bacterial_count = bacterial_read_df_select %>% filter(fileName %in% colnames(city_community_abundance)) %>% ungroup() %>%
    column_to_rownames("fileName")
  bacterial_count = bacterial_count[colnames(city_community_abundance),]
  
  city_community_abundance_norm = ((city_community_abundance) / bacterial_count$total_reads) *1e06
  community_temporal_df = community_temporal_df %>% bind_rows(melt(city_community_abundance_norm %>% mutate(community = groups_table$ind)))
}

groups_df_label = groups_df %>%
  dplyr::select(city, ind, count) %>% distinct() %>%
  arrange(city, desc(count)) %>%
  mutate(city_community = paste0(city, "_", ind)) %>%
  mutate(city_community_new = paste0(c(rep('C',6),rep('D',6),rep('K',6),rep('M',6)), rep(1:6,4))) %>%
  dplyr::select(city_community, city_community_new)

community_temporal_df = community_temporal_df %>% separate(variable, into = c("city", "month"), sep = "_") %>%
  mutate(city_community = paste0(city, "_", community)) %>%
  left_join(groups_df_label)

groups_df = groups_df %>% mutate(city_community = paste0(city, "_", ind)) %>%
  left_join(groups_df_label)

community_temporal_df_plot = community_temporal_df %>%  
  left_join(groups_df %>% dplyr::select(city_community_new, count) %>% distinct(), 
            by = c("city_community_new")) %>%
  mutate(month = factor(month, level = levels(meta$month_label))) %>%
  mutate(city = factor(city, level = levels(meta$city_state))) %>%
  arrange(city, desc(count)) %>%
  mutate(city_community_new = factor(city_community_new, levels = unique(city_community_new)))



community_temporal_df_plot_color = setNames(c(rep(paletteer_dynamic("ggthemes_ptol::qualitative", 12)[c(1:6)],4)),
                                            levels(community_temporal_df_plot$city_community_new))

#Ploting mean abundance of community derived from member nodes/MAGs
community_temporal_df_plot %>%
  mutate(label = ifelse(month == "03/24", as.character(city_community_new), "")) %>%
  ggplot(aes(x = month, y = value)) +
  geom_smooth(aes(group = city_community_new, linewidth = count, color = fct_inorder(city_community_new)), se = F, alpha = 1) +
  scale_linewidth(range = c(2,5)) +
  theme_classic(25) +
  scale_color_manual(values = community_temporal_df_plot_color, guide = F) +
  labs(linewidth = "Nodes count", title = "", color = "Community") +
  ylab(label = "Mean abundance per million total count") +
  xlab(label = "Collection month") +
  theme(axis.text.x=element_text(angle=90, hjust = 0, size = 20), 
        axis.line = element_line(),
        axis.text.y=element_text(size = 20),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank(), 
        panel.spacing.y = unit(1,"cm")) +
  facet_wrap(~city, scale = "free_y", ncol = 1, axes = "all_x", axis.labels = "margins")

groups_col = groups_df %>% dplyr::select(phylum, color) %>% distinct()
groups_col = setNames(groups_col$color, groups_col$phylum)
#Ploting Phylum composition of the communities
groups_df %>% 
  group_by(city_community_new, phylum) %>% mutate(bin.count = n()) %>% ungroup() %>%
  dplyr::select(-values, -novel) %>% distinct() %>%
  group_by(phylum, city) %>% mutate(total_bin = sum(bin.count)) %>% ungroup() %>%
  arrange(desc(total_bin)) %>%
  mutate(city = factor(city, level = levels(meta$city_state))) %>%
  ggplot(aes(y = city_community_new, x = bin.count, fill = fct_inorder(phylum), group = city_community_new)) +
  geom_bar(stat = "identity", position = position_fill(), color = "grey20") +
  scale_fill_manual(values = groups_col) +
  labs(fill = "Phylum") +
  coord_polar(theta = "y") +
  facet_wrap(~city, scale = "free_y", ncol = 2) +
  theme_minimal(25)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.position = "none")


# Correlation matrix between abundance vector of community
community_temporal_df_plot_cor = community_temporal_df_plot %>%
  dplyr::select(city_community_new, month, value) %>%
  pivot_wider(id_cols = city_community_new, values_from = value, names_from = month, values_fill = 0) %>%
  column_to_rownames("city_community_new") %>%
  t() %>%
  cor()
diag(community_temporal_df_plot_cor) = 0

# Phylum composition of community
phylum_community = groups_df %>%
  dplyr::select(phylum, city_community_new) %>%
  filter(!(is.na(phylum))) %>%
  group_by(phylum, city_community_new) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(city_community_new = factor(city_community_new, levels = unique(groups_df$city_community_new))) %>%
  arrange(city_community_new) %>%
  pivot_wider(id_cols = phylum, names_from = city_community_new, values_from = count, values_fill = 0) %>%
  column_to_rownames("phylum")

bray_curtis_dm = 1 - data.matrix(vegan::vegdist(t(phylum_community), method = "bray"))
diag(bray_curtis_dm) = 0
bray_curtis_hclust_obj = hclust(as.dist(1 - bray_curtis_dm))
bray_curtis_hclust_ordering = bray_curtis_hclust_obj$labels[bray_curtis_hclust_obj$order]

ggcorrplot::ggcorrplot(bray_curtis_dm, hc.order = T, show.diag = T, type = "lower") +
  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), labels=c(0, 1), mid = "white", high = "darkgreen", name = "1-bray curtis") +
  theme_bw(20) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle=90, hjust = 1))

ggcorrplot::ggcorrplot(community_temporal_df_plot_cor[bray_curtis_hclust_ordering, bray_curtis_hclust_ordering], show.diag = T, type = "upper") +
  scale_fill_gradient2(breaks=c(-1, 1), limit=c(-1, 1), mid = "white", high = "red3", low = "blue2", name = "Pearson correlation") +
  theme_bw(20) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle=90, hjust = 1))

# Community AMR proportion
total_contigs = checkm2 %>%
  dplyr::select(Name, Total_Contigs)
amr_contigs = amr_drep_checkm2 %>%
  dplyr::select(Name, name_contig) %>%
  distinct() %>%
  group_by(Name) %>%
  summarise(contig_with_amr_count = n())

mge_contigs = contig_amr_genomad_transposon %>%
  filter(plasmid_transposon == "Contigs with MGE") %>%
  dplyr::select(Name, name_contig) %>%
  distinct() %>%
  group_by(Name) %>%
  summarise(contig_with_mge_count = n())

groups_df %>%
  dplyr::select(city, values, city_community_new) %>% 
  dplyr::rename(Name = values) %>% 
  left_join(total_contigs) %>%
  left_join(amr_contigs) %>%
  left_join(mge_contigs) %>%
  replace(is.na(.), 0) %>%
  group_by(city_community_new, city) %>%
  summarise(total_contigs = sum(Total_Contigs),
            contig_with_amr_count = sum(contig_with_amr_count),
            contig_with_mge_count = sum(contig_with_mge_count),
            .groups = "drop") %>%
  mutate(amr_contig_prop = contig_with_amr_count/ total_contigs,
         mge_contig_prop = contig_with_mge_count/ total_contigs) %>%
  mutate(city_community_new = factor(city_community_new, levels = bray_curtis_hclust_ordering)) %>%
  mutate(city = factor(city, levels = levels(meta$city_state))) %>%
  ggplot(aes(y = (city_community_new), x = amr_contig_prop)) +
  geom_col(aes(fill = city),stat = "identity", position = "dodge", color = "black") +
  theme_bw(25) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = 'Dark2', name = "") +
  xlab(paste0("Proportion of Contigs with ARG")) +
  ylab(paste0("Communities"))+
  scale_y_discrete(limits=rev)



# Network stats
degree_network = tibble(bin = bat_drep_network$Name) %>%
  left_join(degree(network_list$Delhi) %>% stack() %>% dplyr::rename(c(Delhi = values, bin = ind))) %>%
  left_join(degree(network_list$Kolkata) %>% stack() %>% dplyr::rename(c(Kolkata = values, bin = ind))) %>%
  left_join(degree(network_list$Chennai) %>% stack() %>% dplyr::rename(c(Chennai = values, bin = ind))) %>%
  left_join(degree(network_list$Mumbai) %>% stack() %>% dplyr::rename(c(Mumbai = values, bin = ind))) %>%
  column_to_rownames('bin')
betweenness_network = tibble(bin = bat_drep_network$Name) %>%
  left_join(betweenness(network_list$Delhi) %>% stack() %>% dplyr::rename(c(Delhi = values, bin = ind))) %>%
  left_join(betweenness(network_list$Kolkata) %>% stack() %>% dplyr::rename(c(Kolkata = values, bin = ind))) %>%
  left_join(betweenness(network_list$Chennai) %>% stack() %>% dplyr::rename(c(Chennai = values, bin = ind))) %>%
  left_join(betweenness(network_list$Mumbai) %>% stack() %>% dplyr::rename(c(Mumbai = values, bin = ind))) %>%
  column_to_rownames('bin')
closeness_network = tibble(bin = bat_drep_network$Name) %>%
  left_join(closeness(network_list$Delhi) %>% stack() %>% dplyr::rename(c(Delhi = values, bin = ind))) %>%
  left_join(closeness(network_list$Kolkata) %>% stack() %>% dplyr::rename(c(Kolkata = values, bin = ind))) %>%
  left_join(closeness(network_list$Chennai) %>% stack() %>% dplyr::rename(c(Chennai = values, bin = ind))) %>%
  left_join(closeness(network_list$Mumbai) %>% stack() %>% dplyr::rename(c(Mumbai = values, bin = ind))) %>%
  column_to_rownames('bin')

network_stat_df = tibble(city = character(), singleton = numeric(), 
                         max_degree = numeric(), min_degree = numeric(), mean_degree = numeric(),
                         max_between = numeric(), min_between = numeric(), mean_between = numeric(),
                         max_close = numeric(), min_close = numeric(), mean_close = numeric(),
                         community_number = numeric(),
                         community_modularity = numeric(),
                         max_member = numeric(), min_member = numeric(), mean_member = numeric())

for (i in colnames(degree_network))
{
  community_selected = clp_list[[i]]
  modularity = community_selected$modularity
  membership_df = lapply(groups(community_selected), length) %>% stack()
  network_stat_df = network_stat_df %>% 
    bind_rows(tibble(city = i, 
                     singleton = degree_network %>% filter(get({{i}}) == 1) %>% nrow(),
                     max_degree = degree_network %>% pull(i) %>% max(na.rm = T),
                     min_degree = degree_network %>% pull(i) %>% min(na.rm = T),
                     mean_degree = degree_network %>% pull(i) %>% mean(na.rm = T),
                     max_between = betweenness_network %>% pull(i) %>% max(na.rm = T),
                     min_between = betweenness_network %>% pull(i) %>% min(na.rm = T),
                     mean_between = betweenness_network %>% pull(i) %>% mean(na.rm = T),
                     max_close = closeness_network %>% pull(i) %>% max(na.rm = T),
                     min_close = closeness_network %>% pull(i) %>% min(na.rm = T),
                     mean_close = closeness_network %>% pull(i) %>% mean(na.rm = T),
                     community_number = nrow(membership_df),
                     community_modularity = modularity,
                     max_member = max(membership_df$values),
                     min_member = min(membership_df$values),
                     mean_member = mean(membership_df$values)
    ))
}

network_stat_df_plot = network_stat_df %>%
  dplyr::select(city,mean_degree, mean_between, mean_close, community_number, community_modularity, max_member) %>%
  mutate_if(is.numeric, round, 2)

colnames(network_stat_df_plot) = c("City", "Missing \nNodes", "Mean \nDegree","Mean \nBetweenness", "Mean \nCloseness", 
                                   "Number of \ncommunity", "Community \nModularity", "Node in \nbiggest \ncommunity")
as.data.frame(network_stat_df_plot) %>% gt::gt()

