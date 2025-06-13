library(tidyverse)
library(janitor)
library(reshape2)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(ggwordcloud)
library(abdiv)
library(ggalt)
library(paletteer)
library(ggpubr)
library(tidyplots)
library(ggsankey)


bar_plot_top_list <- function(dataframe){
  top_all = dataframe %>%
    group_by(city_state, ID) %>%
    summarize(sum_count_data = sum(count_data)) %>%
    ungroup() %>%
    group_by(city_state) %>% 
    arrange(city_state, desc(sum_count_data)) %>%
    top_n(10,wt = sum_count_data)
  
  top_all_label = top_all %>% 
    pivot_wider(names_from = city_state, values_from = city_state, id_cols = ID, values_fill = "") %>%
    mutate(label = ID)
  
  top_all_color = setNames(col_vector[1:length(levels(as.factor(top_all_label$label)))],levels(as.factor(top_all_label$label)))
  top_all_order = top_all %>%
    left_join(top_all_label) %>%
    group_by(label) %>% summarise(sum_count_data = sum(sum_count_data)) %>% 
    arrange(desc(sum_count_data)) %>% ungroup %>%
    mutate(label = as.factor(label))
  
  table_plot = top_all_label %>%
    dplyr::select(-c(ID)) %>%
    melt(id.vars = 'label') %>%
    mutate(label = factor(label, levels = top_all_order$label)) %>%
    filter(!(value == "")) %>%
    ggplot(aes(y = fct_inorder(label), x = variable, fill = label)) +
    geom_point(size = 5, shape = 21) +
    theme_classic(20) +
    scale_fill_manual(values = top_all_color, name = "") +
    scale_y_discrete(limits = rev, label = function(x) str_trunc(x, width = 30)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle=90, hjust = 0, size = 20), 
          axis.text.y=element_text(size = 15),
          legend.position = "none") +
    coord_fixed(ratio = 1) +
    xlab("") +
    ylab("") 
  
  city_top = list()
  for (city in levels(as.factor(dataframe$city_state)))
  {
    top = dataframe %>%
      filter(city_state == city) %>%
      group_by(city_state, ID) %>%
      summarize(sum_count_data = sum(count_data)) %>%
      ungroup() %>%
      group_by(city_state) %>% 
      arrange(city_state, desc(sum_count_data)) %>%
      top_n(10,wt = sum_count_data) %>%
      left_join(top_all_label)
    
    top_order = top %>%
      group_by(label) %>% summarise(sum_count_data = sum(sum_count_data)) %>% 
      arrange(desc(sum_count_data)) %>% ungroup %>%
      mutate(label = as.factor(label))
    
    top$label = factor(top$label,
                       levels = top_order$label)
    
    city_top[[city]] = dataframe %>% 
      filter(city_state == city) %>%
      group_by(city_state, month_label) %>%
      mutate(fraction_count_data = count_data/ sum(count_data)) %>%
      ungroup() %>%
      filter(ID %in% top$ID) %>%
      left_join(top_all_label) %>%
      mutate(label = factor(label, levels = top_order$label)) %>%
      left_join(dataframe %>%
                  filter(city_state == city) %>%
                  group_by(month_label, city_state) %>%
                  summarise_at(vars(count_data), alpha_diversities) %>%
                  dplyr::select(month_label, city_state, simpson, richness, shannon)) %>%
      ggplot(aes(x = as.factor(month_label), y = fraction_count_data, fill = label, group = label)) +
      geom_bar(position="stack", stat="identity") +
      geom_line(aes(y = simpson, group = city_state),color = "red4") +
      geom_point(aes(y = simpson, group = city_state), color = "red4") +
      scale_fill_manual(values = top_all_color, name = "") +
      labs(x = "", y = "", title = city) +
      theme_bw(15) +
      scale_y_continuous(name = "",
                         sec.axis = sec_axis(~./1, name="")) +
      scale_x_discrete(drop = F) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.text.x=element_text(angle=90, vjust=0.3),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      guides(fill=guide_legend(ncol =1))
  }
  p = cowplot::plot_grid(plotlist = city_top) %>%
    ggpubr::annotate_figure(top = text_grob("", color = "black", size = 20, y = 0),
                            bottom = text_grob("Collection month", color = "black",
                                               hjust = 0.5, x = 0.5, , size = 20, y = 1.5),
                            left = text_grob("Fraction of read count", color = "black", rot = 90, size = 20, x = 1.5),
                            right = text_grob("Simpson diversity score", color = "black", rot = 90, size = 20, x = -.5))
  return(cowplot::plot_grid(p, table_plot,  rel_widths = c(3.5, 1.5)))
}

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


priority.pathogen.list = read_delim(file = "data/WHO_indian_priority_pathogen_list.csv", delim = "\t", col_names = F)


k2_nt = read_delim(file = "data/k2_nt_20230502.inspect.species.krona.phylotree.txt", delim = "\t") %>%
  left_join(priority.pathogen.list %>% dplyr::rename(species = X1)) %>% 
  left_join(priority.pathogen.list %>% dplyr::rename(genus = X1), by = "genus") %>% 
  unite(col = "species_genus", c("X2.x", "X2.y"), na.rm = T, remove = T) %>%
  unite(col = "WHO.priority", c("X3.x", "X3.y"), na.rm = T, remove = T) %>%
  unite(col = "Drug_class", c("X4.x", "X4.y"), na.rm = T, remove = T) %>%
  unite(col = "IPPL", c("X5.x", "X5.y"), na.rm = T, remove = T) %>%
  unite(col = "WHOPPL", c("X6.x", "X6.y"), na.rm = T, remove = T) %>%
  mutate_all(~na_if(., '')) %>%
  mutate_all(as.factor)

k2_nt = k2_nt %>% separate(Drug_class, into = c("Drug_class", "Extra"), extra = "merge") %>% dplyr::select(-Extra) %>%
  separate(species_genus, into = c("species_genus", "Extra"), extra = "merge") %>% dplyr::select(-Extra) %>%
  separate(WHO.priority, into = c("WHO.priority", "Extra"), extra = "merge") %>% dplyr::select(-Extra)

meta = read_delim(file = "data/meta.tsv", delim = "\t", col_names = T) %>% 
  clean_names() %>%
  arrange(year, month) %>%
  mutate(month_label = factor(month_label, levels = unique(month_label)))


bracken_phylum = read_delim(file = "data/bracken_phylum.tsv", delim = "\t", col_names = T)
bracken_domain = read_delim(file = "data/bracken_domain.tsv", delim = "\t", col_names = T)
bracken_genus = read_delim(file = "data/bracken_genus.tsv", delim = "\t", col_names = T)
bracken_species = read_delim(file = "data/bracken_species.tsv", delim = "\t", col_names = T)
kraken_report = read_delim(file = "data/kraken_report.tsv", delim = "\t", col_names = T)
bacterial_read_df = read_delim(file = "data/bacterial_read_df.tsv", delim = "\t", col_names = T)
rgi_bowtie = read_delim(file = "data/rgi_bowtie.tsv", delim = "\t", col_names = T)

#########################################################################################
meta_select = meta %>% 
  dplyr::select(city_state, month, year, month_label) %>%
  distinct() %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bacterial_read_df_select = bacterial_read_df %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label) %>%
  summarise(total_reads = sum(total_reads), bacterial_reads = sum(bacterial_reads), total_mapped_reads = sum(total_mapped_reads), eukaryota_reads = sum(eukaryota_reads)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
kraken_report_select = kraken_report %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, X6) %>%
  summarize(X2 = mean(X2)) %>%
  mutate(fileName = paste0(city_state, "_", month_label))
bracken_domain_select = bracken_domain  %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, name) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
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
#########################################################################################

#########################################################################################

bacterial_read_df %>% 
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  dplyr::select(city_state, total_reads, total_mapped_reads, bacterial_reads) %>%
  dplyr::rename(c(`Total reads` = total_reads, `Bacterial reads` = bacterial_reads, `Mapped reads` = total_mapped_reads)) %>%
  melt() %>%
  ggplot(aes(x = variable, y = value, fill = city_state)) +
  geom_boxplot(alpha = 1) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  labs(title="") +
  theme(axis.text.x=element_text(angle=0, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5)) +
  labs(x = "", y = "Read count", fill = "City") +
  theme_classic(20)

bracken_domain_select %>%
  left_join(bacterial_read_df_select) %>%
  mutate(new_est_reads_per_mil_total_mapped_reads = (new_est_reads/ total_mapped_reads)) %>%
  ggplot(aes(x = as.factor(month_label), y = new_est_reads_per_mil_total_mapped_reads, color = name)) +
  geom_point(alpha = 1, shape = 21, size =2) +
  geom_xspline(aes(group = name, color = name)) +
  scale_fill_viridis(discrete = TRUE) + 
  scale_color_hue(l = 25) +
  theme_bw(15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(vjust = 0.5, hjust=0.5)) +
  facet_wrap(~city_state) +
  labs(x = "Collection month", y = "Percentage of total mapped reads", 
       title = "Percentage of reads mapping to Domains",
       fill = "Domains", color = "Domains")

#####################################################################################################################
#     Bacterial phylums
#####################################################################################################################

bracken_phylum_select_bacterial = bracken_phylum_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(phylum) %>% distinct(), by = c("name" = "phylum")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(new_est_reads_per_mil_bacterial_reads = (new_est_reads/ bacterial_reads) * 1e06)

bar_plot_top_list(bracken_phylum_select_bacterial %>%
                    dplyr::rename(c(ID = name, count_data = new_est_reads_per_mil_bacterial_reads)))

#########################################################################################################################################################################
# Genus level analysis
#########################################################################################################################################################################

bracken_genus_select_bacterial = bracken_genus_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(genus) %>% distinct(), by = c("name" = "genus")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(new_est_reads_per_mil_bacterial_reads = (new_est_reads/ bacterial_reads) * 1e06)

bar_plot_top_list(bracken_genus_select_bacterial %>%
                    dplyr::rename(c(ID = name, count_data = new_est_reads_per_mil_bacterial_reads)))

####################################################################################################################
# Bacterial species abundance and differential analysis
####################################################################################################################

bracken_species_select_bacterial = bracken_species_select %>%
  inner_join(k2_nt %>% filter(kingdom == "Bacteria") %>% dplyr::select(species) %>% distinct(), by = c("name" = "species")) %>%
  left_join(bacterial_read_df_select) %>%
  mutate(new_est_reads_per_mil_bacterial_reads = (new_est_reads/ bacterial_reads) * 1e06)

bracken_species_select_bacterial %>%
  group_by(month_label, city_state) %>%
  summarise_at(vars(new_est_reads_per_mil_bacterial_reads), alpha_diversities) %>%
  dplyr::select(month_label, city_state, richness, simpson, shannon) %>%
  melt(id.vars = c('month_label', 'city_state')) %>%
  ggplot(aes(x = as.factor(month_label), y = value, fill = city_state, group = city_state)) +
  geom_point(size = 2, shape = 21, show.legend = T, alpha = 1) +
  geom_xspline(size = 1, aes(color = city_state), show.legend = T, alpha = 0.7) +
  scale_color_brewer(palette = 'Dark2', name = "City") +
  scale_fill_brewer(palette = 'Dark2', name = "City") +
  theme_classic(25)  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="Collection month", y="", title="") +
  facet_wrap(~variable, scales = "free_y", nrow = 3) %>%
  print()


#https://microbiome.github.io/course_2021_radboud/beta-diversity.html
bracken_species_dm = bracken_species_select_bacterial %>%
  pivot_wider(id_cols = name, names_from = fileName, values_from = new_est_reads_per_mil_bacterial_reads, values_fill = 0) %>% 
  column_to_rownames(var = 'name')
bray_curtis_dist <- vegan::vegdist(t(bracken_species_dm), method = "bray")
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
tibble(pcoa1 = bray_curtis_pcoa$vectors[,1],
       pcoa2 = bray_curtis_pcoa$vectors[,2]) %>%
  mutate(fileName = rownames(bray_curtis_pcoa$vectors)) %>%
  left_join(meta_select) %>%
  left_join(bracken_genus_select_bacterial %>%
              group_by(month_label, city_state) %>%
              summarise_at(vars(new_est_reads_per_mil_bacterial_reads), alpha_diversities) %>%
              dplyr::select(month_label, city_state, shannon)) %>%
  ggplot(aes(x = pcoa1, y= pcoa2)) +
  geom_point(shape = 21, aes(fill = as.factor(city_state), size = shannon)) +
  scale_size_continuous(range = c(2,10), transform = "reverse") +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  labs(title="PCoA using bray curtis distance between samples by species abundance") +
  theme_bw(30) +
  guides(fill = guide_legend(override.aes = list(shape=21, size = 5))) +
  theme(axis.text.x=element_text(angle=0, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5)) +
  labs(fill = "City", size = "Shannon diversity") +
  ggforce::geom_mark_ellipse(aes(color = as.factor(city_state), fill = as.factor(city_state)), alpha = 0.05,show.legend = F)

bray_curtis_dist <- vegan::vegdist(t(bracken_species_dm), method = "aitchison", pseudocount = 0.0000001)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
tibble(pcoa1 = bray_curtis_pcoa$vectors[,1],
       pcoa2 = bray_curtis_pcoa$vectors[,2]) %>%
  mutate(fileName = rownames(bray_curtis_pcoa$vectors)) %>%
  left_join(meta_select) %>%
  left_join(bracken_genus_select_bacterial %>%
              group_by(month_label, city_state) %>%
              summarise_at(vars(new_est_reads_per_mil_bacterial_reads), alpha_diversities) %>%
              dplyr::select(month_label, city_state, shannon)) %>%
  ggplot(aes(x = pcoa1, y= pcoa2)) +
  geom_point(shape = 21, aes(fill = as.factor(city_state)), size = 7) +
  scale_size_continuous(range = c(1,10), transform = "reverse") +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  labs(title="") +
  theme_classic(25) + 
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme(axis.text.x=element_text(angle=0, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5)) +
  labs(fill = "City", size = "Shannon diversity") +
  ggforce::geom_mark_ellipse(aes(color = as.factor(city_state),fill = as.factor(city_state)), alpha = 0.1,show.legend = F)

bar_plot_top_list(bracken_species_select_bacterial %>%
                    dplyr::rename(c(ID = name, count_data = new_est_reads_per_mil_bacterial_reads)))

bracken_species_select_bacterial_meta = bracken_species_select_bacterial %>%
  dplyr::select(fileName, city_state, month_label) %>%
  distinct()

bracken_species_select_bacterial_df = bracken_species_select_bacterial %>%
  filter(fileName %in% bracken_species_select_bacterial_meta$fileName) %>%
  group_by(name, fileName) %>% summarise(new_est_reads = sum(new_est_reads)) %>% ungroup() %>% #This is to bypass when there are two instance of the genus
  tidyr::pivot_wider(id_cols = name, names_from = fileName, values_from = new_est_reads, values_fill = 0) %>%
  column_to_rownames('name')
bracken_species_select_bacterial_df = bracken_species_select_bacterial_df[,bracken_species_select_bacterial_meta$fileName]

dds = DESeqDataSetFromMatrix(countData = bracken_species_select_bacterial_df,
                                  colData = bracken_species_select_bacterial_meta, design = ~0+month_label+city_state)
dds <- dds[rowSums(counts(dds))>10,]
dds <- DESeq(dds, parallel=T)
resultsNames(dds)

log2Foldchange_cutoff = 1

deg_list_Chennai_location_up = deg_list_Delhi_location_up = deg_list_Kolkata_location_up = deg_list_Mumbai_location_up = list()
deg_list_Chennai_location_dn = deg_list_Delhi_location_dn = deg_list_Kolkata_location_dn = deg_list_Mumbai_location_dn = list()

number_of_city = 4 
city_list = c("Kolkata", "Mumbai", "Chennai", "Delhi")

city = "Chennai"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Chennai_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Chennai_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}
city = "Delhi"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Delhi_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                             ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Delhi_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                             ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
    }
}
city = "Kolkata"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Kolkata_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Kolkata_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}
city = "Mumbai"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Mumbai_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                              ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Mumbai_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                              ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}

deg_list_Chennai_location_intersect = deg_list_Chennai_location_up %>% unlist %>% stack() %>% dplyr::rename('Chennai' = values) %>% group_by(`Chennai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Chennai')
deg_list_Delhi_location_intersect = deg_list_Delhi_location_up %>% unlist %>% stack() %>% dplyr::rename("Delhi" = values) %>% group_by(`Delhi`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Delhi')
deg_list_Kolkata_location_intersect = deg_list_Kolkata_location_up %>% unlist %>% stack() %>% dplyr::rename("Kolkata" = values) %>% group_by(`Kolkata`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Kolkata')
deg_list_Mumbai_location_intersect = deg_list_Mumbai_location_up %>% unlist %>% stack() %>% dplyr::rename("Mumbai" = values) %>% group_by(`Mumbai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Mumbai')


bracken_species_select_bacterial_df_norm = bracken_species_select_bacterial %>%
  filter(fileName %in% bracken_species_select_bacterial_meta$fileName) %>%
  group_by(name, fileName) %>% summarise(new_est_reads_per_mil_bacterial_reads = sum(new_est_reads_per_mil_bacterial_reads)) %>% ungroup() %>% #This is to bypass when there are two instance of the genus
  tidyr::pivot_wider(id_cols = name, names_from = fileName, values_from = new_est_reads_per_mil_bacterial_reads, values_fill = 0) %>%
  column_to_rownames('name')
bracken_species_select_bacterial_df_norm = bracken_species_select_bacterial_df_norm[,bracken_species_select_bacterial_meta$fileName]


bracken_species_select_bacterial_intersect_df = bracken_species_select_bacterial_df_norm %>% filter((rownames(.) %in% deg_list_Chennai_location_intersect)) %>%
  bind_rows(bracken_species_select_bacterial_df_norm %>% filter((rownames(.) %in% deg_list_Delhi_location_intersect))) %>%
  bind_rows(bracken_species_select_bacterial_df_norm %>% filter((rownames(.) %in% deg_list_Kolkata_location_intersect))) %>%
  bind_rows(bracken_species_select_bacterial_df_norm %>% filter((rownames(.) %in% deg_list_Mumbai_location_intersect)))

annotation_row = tibble("city" = "Chennai", "DEG" = deg_list_Chennai_location_intersect) %>%
  bind_rows(tibble("city" = "Delhi", "DEG" = deg_list_Delhi_location_intersect)) %>%
  bind_rows(tibble("city" = "Kolkata", "DEG" = deg_list_Kolkata_location_intersect)) %>%
  bind_rows(tibble("city" = "Mumbai", "DEG" = deg_list_Mumbai_location_intersect)) %>%
  column_to_rownames("DEG") %>%
  arrange(city) %>%
  mutate_all(as.factor)
annotation_column = tibble(column = colnames(bracken_species_select_bacterial_intersect_df)) %>% 
  separate(column, into = c("city", "date"), sep = "_", remove = F) %>% 
  column_to_rownames("column") %>% dplyr::select(c("city", "date")) %>%
  mutate_all(as.factor)

ann_colors = list(
  city = setNames(brewer.pal(n = length(levels(annotation_column$city)), name = 'Dark2'), levels(annotation_column$city)))

ComplexHeatmap::pheatmap(bracken_species_select_bacterial_intersect_df[rownames(annotation_row), rownames(annotation_column)] %>% as.matrix(), 
                         main = paste0(""), 
                         annotation_row = annotation_row , 
                         annotation_col =  annotation_column %>% dplyr::select(-date),
                         show_rownames = F, 
                         show_colnames = F, 
                         cluster_cols = F, cluster_rows = T,
                         scale = "row", 
                         color = colorRampPalette((c("blue", "white", "red4")))(100),
                         labels_col = annotation_column$date,
                         labels_row = annotation_row$label,
                         row_split = as.character(annotation_row$city), row_title = NULL, 
                         show_row_dend = F,
                         gaps_col = head(as.numeric(cumsum(table(annotation_column$city))), - 1),
                         annotation_colors = ann_colors, 
                         cellwidth = 6, cellheight = .2,
                         legend = T, fontsize = 15)

deg_list_Chennai_location_intersect = deg_list_Chennai_location_dn %>% unlist %>% stack() %>% dplyr::rename('Chennai' = values) %>% group_by(`Chennai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Chennai')
deg_list_Delhi_location_intersect = deg_list_Delhi_location_dn %>% unlist %>% stack() %>% dplyr::rename("Delhi" = values) %>% group_by(`Delhi`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Delhi')
deg_list_Kolkata_location_intersect = deg_list_Kolkata_location_dn %>% unlist %>% stack() %>% dplyr::rename("Kolkata" = values) %>% group_by(`Kolkata`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('Kolkata')
deg_list_Mumbai_location_intersect = deg_list_Mumbai_location_dn %>% unlist %>% stack() %>% dplyr::rename("A1 Mumbai" = values) %>% group_by(`A1 Mumbai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull('A1 Mumbai')


bracken_species_select_bacterial_intersect_df = bracken_species_select_bacterial_df %>% filter((rownames(.) %in% deg_list_Chennai_location_intersect)) %>%
  bind_rows(bracken_species_select_bacterial_df %>% filter((rownames(.) %in% deg_list_Delhi_location_intersect))) %>%
  bind_rows(bracken_species_select_bacterial_df %>% filter((rownames(.) %in% deg_list_Kolkata_location_intersect))) %>%
  bind_rows(bracken_species_select_bacterial_df %>% filter((rownames(.) %in% deg_list_Mumbai_location_intersect)))

annotation_row = tibble("city" = "Chennai", "DEG" = deg_list_Chennai_location_intersect) %>%
  bind_rows(tibble("city" = "Delhi", "DEG" = deg_list_Delhi_location_intersect)) %>%
  bind_rows(tibble("city" = "Kolkata", "DEG" = deg_list_Kolkata_location_intersect)) %>%
  bind_rows(tibble("city" = "Mumbai", "DEG" = deg_list_Mumbai_location_intersect)) %>%
  column_to_rownames("DEG") %>%
  arrange(city) %>%
  mutate_all(as.factor)
annotation_column = tibble(column = colnames(bracken_species_select_bacterial_intersect_df)) %>% 
  separate(column, into = c("city", "date"), sep = "_", remove = F) %>% 
  column_to_rownames("column") %>% dplyr::select(c("city", "date"))%>%
  mutate_all(as.factor)

ann_colors = list(
  city = setNames(brewer.pal(n = length(levels(annotation_column$city)), name = 'Dark2'), levels(annotation_column$city)))

ComplexHeatmap::pheatmap(bracken_species_select_bacterial_intersect_df[rownames(annotation_row), rownames(annotation_column)]  %>% as.matrix(), 
                         main = paste0("Differentially depleted species"), 
                         annotation_row = annotation_row , 
                         annotation_col =  annotation_column %>% dplyr::select(-date),
                         show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = T,
                         scale = "row", 
                         color = colorRampPalette((c("blue", "white", "red4")))(100),
                         labels_col = annotation_column$date,
                         row_split = as.character(annotation_row$city), row_title = NULL, 
                         show_row_dend = F,
                         labels_row = annotation_row$label,
                         gaps_col = head(as.numeric(cumsum(table(annotation_column$city))), - 1),
                         annotation_colors = ann_colors, 
                         cellwidth = 4, cellheight = .2,
                         legend = T)


species_list_city = deg_list_Chennai_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
  filter(values ==3) %>%
  mutate(city = "Chennai", regulation = "Enriched") %>%
  
  bind_rows(deg_list_Chennai_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Chennai", regulation = "Depleted")) %>%
  bind_rows(deg_list_Delhi_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Delhi", regulation = "Enriched")) %>%
  bind_rows(deg_list_Delhi_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Delhi", regulation = "Depleted")) %>%
  bind_rows(deg_list_Kolkata_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Kolkata", regulation = "Enriched")) %>%
  bind_rows(deg_list_Kolkata_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Kolkata", regulation = "Depleted")) %>%
  bind_rows(deg_list_Mumbai_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Mumbai", regulation = "Enriched")) %>%
  bind_rows(deg_list_Mumbai_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% 
              mutate(city = "Mumbai", regulation = "Depleted"))

species_list_city_abundance = bracken_species_select_bacterial %>%
  filter(name %in% species_list_city$ind) %>%
  group_by(name, city_state) %>%
  summarise(new_est_reads_per_mil_bacterial_reads = 1e06*(sum(new_est_reads)/sum(bacterial_reads))) %>%
  ungroup() %>% left_join(species_list_city %>% dplyr::rename(name = ind, city_state = city)) %>%
  filter(!(is.na(regulation)))

species_list_city_abundance %>%
  mutate(city_state = factor(city_state, levels = c("Chennai", "Delhi", "Kolkata", "Mumbai"))) %>%
  ggplot(aes(label = name, color = regulation, size = (new_est_reads_per_mil_bacterial_reads))) +
  geom_text_wordcloud(show.legend = T) +
  scale_size_area(max_size = 10) +
  theme_bw() +
  ggtitle(label = "Enriched/Depleted bacterial species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(size = "Average abundance") +
  facet_wrap(~city_state, scales = "free") +
  scale_color_manual(values = c("blue", "red4"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# trimmed rgi, keeping only those aro for which we have Average MAPQ >10 and completely mapped reads more than 10
rgi_bowtie_trimmed = rgi_bowtie %>% 
  filter(`Reference Allele(s) Identity to CARD Reference Protein (%)` >= 90) %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, `ARO Term`, `Drug Class`) %>%
  summarise(`All Mapped Reads` = sum(`All Mapped Reads`), .groups = 'drop') %>%
  mutate(fileName = paste0(city_state, "_", month_label)) %>%
  left_join(bacterial_read_df_select) %>%
  clean_names() %>%
  mutate(all_mapped_read_per_mil_total_reads = (all_mapped_reads/total_reads) * 1e06) %>%
  filter(all_mapped_read_per_mil_total_reads > 10)

rgi_bowtie_trimmed %>%
  group_by(month_label, city_state) %>%
  summarise_at(vars(all_mapped_read_per_mil_total_reads), alpha_diversities) %>%
  dplyr::select(month_label, city_state, richness, simpson, shannon) %>%
  melt(id.vars = c('month_label', 'city_state')) %>%
  ggplot(aes(x = as.factor(month_label), y = value, fill = city_state, group = city_state)) +
  geom_point(size = 2, shape = 21, show.legend = T) +
  geom_xspline(size = 1, alpha = 1, aes(color = city_state), show.legend = T) +
  scale_fill_brewer(palette = 'Dark2', name = "City") +
  scale_color_brewer(palette = 'Dark2', name = "City") +
  theme_classic(20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="Collection month", y="", title="") +
  facet_wrap(~variable, scales = "free_y", nrow = 3) %>%
  print()


#https://microbiome.github.io/course_2021_radboud/beta-diversity.html
rgi_bowtie_trimmed_dm = rgi_bowtie_trimmed %>%
  pivot_wider(id_cols = aro_term, names_from = file_name, values_from = all_mapped_read_per_mil_total_reads, values_fill = 0) %>% 
  column_to_rownames(var = 'aro_term')

bray_curtis_dist <- vegan::vegdist(t(rgi_bowtie_trimmed_dm), method = "aitchison", pseudocount = 0.000001)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
tibble(pcoa1 = bray_curtis_pcoa$vectors[,1],
       pcoa2 = bray_curtis_pcoa$vectors[,2]) %>%
  mutate(fileName = rownames(bray_curtis_pcoa$vectors)) %>%
  left_join(meta_select) %>%
  ggplot(aes(x = pcoa1, y= pcoa2)) +
  geom_point(size=7, aes(fill = as.factor(city_state)), shape = 21) +
  labs(title="") +
  ggforce::geom_mark_ellipse(aes(color = paste(as.factor(city_state)), fill =  paste(as.factor(city_state))), alpha = 0.1,show.legend = F) +
  scale_fill_brewer(palette = 'Dark2', name = "City") +
  scale_color_brewer(palette = 'Dark2', name = "City") +
  theme_bw(50) + 
  guides(fill = guide_legend(override.aes = list(shape=21)),
         color = guide_legend(override.aes = list(shape=21))) +
  theme(axis.text.x=element_text(angle=0, vjust=0.3),
        axis.text.y=element_text(),
        plot.title=element_text(hjust = 0.5)) +
  #facet_wrap(~location) +
  labs(fill = "City") %>%
  print()

bar_plot_top_list(rgi_bowtie_trimmed %>%
                    dplyr::rename(c(ID = aro_term, count_data = all_mapped_read_per_mil_total_reads)))


rgi_bowtie_trimmed_drug_class = rgi_bowtie_trimmed %>%
  mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
  mutate(drug_class_original = drug_class) %>%
  mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug",
                                TRUE ~ drug_class))
bar_plot_top_list(rgi_bowtie_trimmed_drug_class %>%
                    dplyr::rename(c(ID = drug_class, count_data = all_mapped_read_per_mil_total_reads)))


for (city in unique(rgi_bowtie_trimmed_drug_class$city_state))
{
  print(
  rgi_bowtie_trimmed_drug_class %>% filter(drug_class %in% "multidrug") %>% ungroup() %>%
    dplyr::select(city_state, drug_class_original) %>%
    filter(city_state == city) %>%
    mutate(drug_class_original = str_replace_all(drug_class_original, pattern = "; ", replacement = "&")) %>%
    group_by(drug_class_original) %>%
    summarise(count = n()) %>%
    filter(count >= 10) %>%
    deframe() %>%
    UpSetR::fromExpression() %>%
    UpSetR::upset(matrix.dot.alpha = 0.1,
                  nintersects = NA, 
                  nsets = 200, 
                  order.by = "freq", 
                  decreasing = T, 
                  mb.ratio = c(0.4, 0.4),
                  number.angles = 0, 
                  text.scale = 3, 
                  point.size = 4, 
                  line.size = 1, sets.x.label = city
    ))
}

#################################################################################################################
## ARG Differential analysis
#################################################################################################################

rgi_bowtie_trimmed_meta = rgi_bowtie_trimmed %>% 
  ungroup() %>%
  dplyr::select(file_name, city_state, month_label) %>%
  distinct()

rgi_bowtie_trimmed_df = rgi_bowtie_trimmed %>%
  filter(file_name %in% rgi_bowtie_trimmed_meta$file_name) %>%
  tidyr::pivot_wider(id_cols = aro_term, names_from = file_name, values_from = all_mapped_reads, values_fill = 0) %>%
  column_to_rownames('aro_term')
rgi_bowtie_trimmed_df = rgi_bowtie_trimmed_df[,rgi_bowtie_trimmed_meta$file_name]

dds_arg = DESeqDataSetFromMatrix(countData = rgi_bowtie_trimmed_df,
                             colData = rgi_bowtie_trimmed_meta, design = ~0+month_label+city_state)
dds_arg <- dds_arg[rowSums(counts(dds_arg))>10,]
dds_arg <- DESeq(dds_arg, parallel=T)
resultsNames(dds)

log2Foldchange_cutoff = 1


deg_list_Chennai_location_up = deg_list_Delhi_location_up = deg_list_Kolkata_location_up = deg_list_Mumbai_location_up = list()
deg_list_Chennai_location_dn = deg_list_Delhi_location_dn = deg_list_Kolkata_location_dn = deg_list_Mumbai_location_dn = list()

number_of_city = 4 
city_list = c("Kolkata", "Mumbai", "Chennai", "Delhi")

city = "Chennai"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Chennai_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Chennai_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}
city = "Delhi"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Delhi_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                             ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Delhi_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                             ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
    }
}
city = "Kolkata"
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Kolkata_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Kolkata_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                               ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}
city = "Mumbai"
enriched_species = c()
for (city_denominator in city_list)
{
  if(!(city == city_denominator))
  {
    deg_list_Mumbai_location_up[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                              ((padj < 0.05)&((log2FoldChange) > log2Foldchange_cutoff))))
    deg_list_Mumbai_location_dn[[paste0(city_denominator)]] = rownames(subset(results(dds_arg, alpha = 0.05, contrast = c("city_state", city, city_denominator)), 
                                                                              ((padj < 0.05)&((log2FoldChange) < -log2Foldchange_cutoff))))
  }
}

deg_list_Chennai_location_intersect = deg_list_Chennai_location_up %>% unlist %>% stack() %>% dplyr::rename(`Chennai` = values) %>% group_by(`Chennai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Chennai`)
deg_list_Delhi_location_intersect = deg_list_Delhi_location_up %>% unlist %>% stack() %>% dplyr::rename(`Delhi` = values) %>% group_by(`Delhi`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Delhi`)
deg_list_Kolkata_location_intersect = deg_list_Kolkata_location_up %>% unlist %>% stack() %>% dplyr::rename(`Kolkata` = values) %>% group_by(`Kolkata`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Kolkata`)
deg_list_Mumbai_location_intersect = deg_list_Mumbai_location_up %>% unlist %>% stack() %>% dplyr::rename(`Mumbai` = values) %>% group_by(`Mumbai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Mumbai`)


rgi_bowtie_trimmed_df_norm = rgi_bowtie_trimmed %>%
  filter(file_name %in% rgi_bowtie_trimmed_meta$file_name) %>%
  tidyr::pivot_wider(id_cols = aro_term, names_from = file_name, values_from = all_mapped_read_per_mil_total_reads, values_fill = 0) %>%
  column_to_rownames('aro_term')
rgi_bowtie_trimmed_df_norm = rgi_bowtie_trimmed_df_norm[,rgi_bowtie_trimmed_meta$file_name]


rgi_bowtie_trimmed_intersect_df = rgi_bowtie_trimmed_df_norm %>% filter((rownames(.) %in% deg_list_Chennai_location_intersect)) %>%
  bind_rows(rgi_bowtie_trimmed_df_norm %>% filter((rownames(.) %in% deg_list_Delhi_location_intersect))) %>%
  bind_rows(rgi_bowtie_trimmed_df_norm %>% filter((rownames(.) %in% deg_list_Kolkata_location_intersect))) %>%
  bind_rows(rgi_bowtie_trimmed_df_norm %>% filter((rownames(.) %in% deg_list_Mumbai_location_intersect)))

annotation_row = tibble("city" = "Chennai", "DEG" = deg_list_Chennai_location_intersect) %>%
  bind_rows(tibble("city" = "Delhi", "DEG" = deg_list_Delhi_location_intersect)) %>%
  bind_rows(tibble("city" = "Kolkata", "DEG" = deg_list_Kolkata_location_intersect)) %>%
  bind_rows(tibble("city" = "Mumbai", "DEG" = deg_list_Mumbai_location_intersect)) %>%
  left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("DEG" = "aro_term")) %>%
  mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
  mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug",
                                TRUE ~ drug_class)) %>%
  column_to_rownames("DEG") %>%
  arrange(city, drug_class) %>%
  mutate_all(as.factor)
annotation_column = tibble(column = colnames(rgi_bowtie_trimmed_intersect_df)) %>% 
  separate(column, into = c("city", "date"), sep = "_", remove = F) %>% 
  column_to_rownames("column") %>% dplyr::select(c("city", "date"))%>%
  mutate_all(as.factor)

ann_colors = list(
  city = setNames(brewer.pal(n = length(levels(annotation_column$city)), name = 'Dark2'), levels(annotation_column$city)),
  drug_class = setNames(object = paletteer_d("ggsci::default_igv")[1:length(levels(annotation_row$drug_class))], nm = levels(annotation_row$drug_class)))

ann_colors$drug_class["multidrug"] = "#000000"

ComplexHeatmap::pheatmap(rgi_bowtie_trimmed_intersect_df[rownames(annotation_row), rownames(annotation_column)] %>% as.matrix(), 
                         main = paste0("Differentially enriched ARGs"), 
                         annotation_row = annotation_row, annotation_col =  annotation_column %>% dplyr::select(-date),
                         show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                         scale = "row", 
                         color = colorRampPalette((c("blue", "white", "red4")))(100),
                         labels_col = annotation_column$date,
                         gaps_row = head(as.numeric(cumsum(table(annotation_row$city))), - 1),
                         gaps_col = head(as.numeric(cumsum(table(annotation_column$city))), - 1),
                         annotation_colors = ann_colors, 
                         legend = T, fontsize = 12, 
                         cellwidth = 8, cellheight = 4,border_color = "grey99")


deg_list_Chennai_location_intersect = deg_list_Chennai_location_dn %>% unlist %>% stack() %>% dplyr::rename(`Chennai` = values) %>% group_by(`Chennai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Chennai`)
deg_list_Delhi_location_intersect = deg_list_Delhi_location_dn %>% unlist %>% stack() %>% dplyr::rename(`Delhi` = values) %>% group_by(`Delhi`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Delhi`)
deg_list_Kolkata_location_intersect = deg_list_Kolkata_location_dn %>% unlist %>% stack() %>% dplyr::rename(`Kolkata` = values) %>% group_by(`Kolkata`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Kolkata`)
deg_list_Mumbai_location_intersect = deg_list_Mumbai_location_dn %>% unlist %>% stack() %>% dplyr::rename(`Mumbai` = values) %>% group_by(`Mumbai`) %>% summarise(count = n()) %>% filter(count == 3) %>% pull(`Mumbai`)


rgi_bowtie_trimmed_intersect_df = rgi_bowtie_trimmed_df %>% filter((rownames(.) %in% deg_list_Chennai_location_intersect)) %>%
  bind_rows(rgi_bowtie_trimmed_df %>% filter((rownames(.) %in% deg_list_Delhi_location_intersect))) %>%
  bind_rows(rgi_bowtie_trimmed_df %>% filter((rownames(.) %in% deg_list_Kolkata_location_intersect))) %>%
  bind_rows(rgi_bowtie_trimmed_df %>% filter((rownames(.) %in% deg_list_Mumbai_location_intersect)))

annotation_row = tibble("city" = "Chennai", "DEG" = deg_list_Chennai_location_intersect) %>%
  bind_rows(tibble("city" = "Delhi", "DEG" = deg_list_Delhi_location_intersect)) %>%
  bind_rows(tibble("city" = "Kolkata", "DEG" = deg_list_Kolkata_location_intersect)) %>%
  bind_rows(tibble("city" = "Mumbai", "DEG" = deg_list_Mumbai_location_intersect)) %>%
  left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("DEG" = "aro_term")) %>%
  mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
  mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug",
                                TRUE ~ drug_class)) %>%
  column_to_rownames("DEG") %>%
  arrange(city, drug_class) %>%
  mutate_all(as.factor)
annotation_column = tibble(column = colnames(rgi_bowtie_trimmed_intersect_df)) %>% 
  separate(column, into = c("city", "date"), sep = "_", remove = F) %>% 
  column_to_rownames("column") %>% dplyr::select(c("city", "date"))%>%
  mutate_all(as.factor)

ann_colors = list(
  city = setNames(brewer.pal(n = length(levels(annotation_column$city)), name = 'Dark2'), levels(annotation_column$city)),
  drug_class = setNames(object = paletteer_d("ggsci::default_igv")[1:length(levels(annotation_row$drug_class))], nm = levels(annotation_row$drug_class)))

ann_colors$drug_class["multidrug"] = "#000000"

ComplexHeatmap::pheatmap(rgi_bowtie_trimmed_intersect_df[rownames(annotation_row), rownames(annotation_column)] %>% as.matrix(), 
                         main = paste0(""), 
                         annotation_row = annotation_row, annotation_col =  annotation_column %>% dplyr::select(-date),
                         show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                         scale = "row", 
                         color = colorRampPalette((c("blue", "white", "red4")))(100),
                         labels_col = annotation_column$date,
                         #row_split = as.character(annotation_row$city), row_title = NULL, 
                         #show_row_dend = F,
                         gaps_row = head(as.numeric(cumsum(table(annotation_row$city))), - 1),
                         gaps_col = head(as.numeric(cumsum(table(annotation_column$city))), - 1),
                         annotation_colors = ann_colors, 
                         legend = T, fontsize = 12, 
                         cellwidth = 6, cellheight = 6,border_color = "grey99")


arg_list_city = deg_list_Chennai_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
  filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
  mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
  mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
  group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
  mutate(city = "Chennai", regulation = "Enriched") %>%
  
  bind_rows(deg_list_Chennai_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>%
              mutate(city = "Chennai", regulation = "Depleted")) %>%
  bind_rows(deg_list_Delhi_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Delhi", regulation = "Enriched")) %>%
  bind_rows(deg_list_Delhi_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Delhi", regulation = "Depleted")) %>%
  bind_rows(deg_list_Kolkata_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Kolkata", regulation = "Enriched")) %>%
  bind_rows(deg_list_Kolkata_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Kolkata", regulation = "Depleted")) %>%
  bind_rows(deg_list_Mumbai_location_up %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Mumbai", regulation = "Enriched")) %>%
  bind_rows(deg_list_Mumbai_location_dn %>% unlist() %>% table() %>% stack() %>% arrange(desc(values)) %>% tibble() %>% 
              filter(values ==3) %>% left_join(rgi_bowtie_trimmed %>% ungroup() %>% dplyr::select(aro_term, drug_class) %>% distinct(), by = c("ind" = "aro_term")) %>% 
              mutate(drug_class = str_replace_all(drug_class, " antibiotic", "")) %>%
              mutate(drug_class = case_when(stringr::str_count(drug_class, pattern = ";") > 1 ~ "multidrug", TRUE ~ drug_class)) %>%
              group_by(drug_class) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
              mutate(city = "Mumbai", regulation = "Depleted"))


arg_list_city %>%
  mutate(city = factor(city, levels = c("Chennai", "Delhi", "Kolkata", "Mumbai"))) %>%
  ggplot(aes(label = drug_class, color = regulation, size = (count))) +
  geom_text_wordcloud(show.legend = T) +
  scale_size_area(max_size = 20) +
  theme_bw(25) +
  ggtitle(label = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(size = "Number of ARGs", color = "") +
  facet_wrap(~city) +
  scale_color_manual(values = c("blue", "red4"))

########################################################################################################################
########################################################################################################################
#Reading and plotting checkm2 stats 
checkm2 = read_delim("data/Metro_merged_quality_report.tsv", delim = "\t") %>%
  separate(Name, into = c("sample_code", "extra"), sep = "_", remove = F) %>%
  left_join(meta) %>%
  filter(Completeness >= 50, Contamination <= 10)


checkm2 %>%
  group_by(month_label, city_state, sample_code) %>%
  summarise(bin_count = n()) %>%
  left_join(bacterial_read_df, by = c("sample_code" = "fileName")) %>%
  mutate(bin_count_per_total_reads = bin_count/total_reads) %>%
  melt(id.vars = c("city_state","month_label"), measure.vars = c("bin_count_per_total_reads")) %>%
  ggplot(aes(x = city_state, y = value)) +
  geom_boxplot(aes(color = city_state), position=position_dodge(0.75), outlier.shape = NA, width = 0.75) +
  geom_point(aes(fill = city_state), shape = 21, position = position_dodge2(0.75), size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw(20) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme(axis.text.x=element_text(size=12, angle=90, vjust=0.3),
        axis.text.y=element_text(size=12),
        plot.title=element_text(size=18, hjust = 0.5)) +
  labs(color = "City", fill = "City", shape = "City") +
  xlab(paste0("City")) +
  ylab(paste0("Number of bins/ total reads"))

coverM = read_tsv("data/coverM.tsv", col_names = T)

rgi_drep = read_delim(file = "data/rgi_drep.tsv", delim = "\t", col_names = T)
genomad_plasmid = read_delim(file = "data/genomad_plasmid.tsv", delim = "\t", col_names = T)
bat_drep = read_delim(file = "data/bat_drep.tsv", delim = "\t", col_names = T)


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

bat_drep_resolved %>%  
  group_by(city_state, variable, value) %>% summarise(count = n()) %>% ungroup() %>%
  left_join(bat_drep_resolved %>%  
              group_by(city_state, variable) %>% summarise(total = n()) %>% ungroup()) %>%
  mutate(percentage = round((count*100)/total, 2)) %>%
  mutate(value = fct_relevel(value, c("Unresolved MAG", "Resolved MAG"))) %>%
  filter(value == "Resolved MAG") %>%
  ggplot(aes(y = variable, x = percentage, group = city_state)) +
  geom_bar(aes(fill = city_state, color = city_state), position = "dodge", stat = "identity", alpha = 1) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic(25) +
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

bat_drep %>%
  left_join(checkm2) %>%
  left_join(coverM) %>%
  mutate(Log_contig_N50 = log(Contig_N50)) %>%
  mutate(log_trim_mean = log(trim_mean)) %>%
  dplyr::select(species, Completeness, Contamination, Log_contig_N50, Genome_Size, Total_Contigs, log_trim_mean) %>%
  dplyr::rename(c(`log(Contigs N50)` = Log_contig_N50), `Genome Size` = Genome_Size, `Total contigs` = Total_Contigs, `log(trimmed \nmean coverage)` = log_trim_mean) %>%
  mutate(species = case_when(species ==  "no support" ~ "Novel",
                             TRUE ~ "Known")) %>%
  melt(id.vars = c("species")) %>%
  tidyplot(x = species, y = value, color = species, alpha = 1) %>%
  add_violin(dodge_width = 0.2 ) %>%
  add_test_pvalue(ref.group = 1, p.adjust.method = "BH", hide_info = T, label.size = 7,bracket.nudge.y = .2) %>%
  theme_tidyplot(fontsize = 25) %>%
  remove_x_axis_title() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_title() %>%
  adjust_legend_title("MAGs at the\nspecies taxa") %>%
  split_plot(by = variable, widths = 100, heights = 100)

#Sankey plot for bin classification
bat_drep_sankey = bat_drep %>% left_join(checkm2)

bat_drep_df = bat_drep_sankey %>%
  separate(superkingdom, into = c("superkingdom", "superkingdom_score"), sep = ":") %>%
  separate(phylum, into = c("phylum", "phylum_score"), sep = ":") %>%
  separate(class, into = c("class", "class_score"), sep = ":") %>%
  separate(order, into = c("order", "order_score"), sep = ":") %>%
  separate(family, into = c("family", "family_score"), sep = ":") %>%
  separate(genus, into = c("genus", "genus_score"), sep = ":") %>%
  separate(species, into = c("species", "species_score"), sep = ":") %>%
  dplyr::select(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "city_state"))
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
  geom_sankey_label(size = 4.2, color = "gray30", fill = "white") +
  theme_sankey(base_size = 15) +
  labs(x = NULL) +
  theme(legend.position = "NULL",
        plot.title = element_text(hjust = .5)) +
  facet_wrap(~city)

############################################################################################################################################################################
############################################################################################################################################################################

genomad_plasmid = genomad_plasmid %>% filter(plasmid_score >= 0.9) %>% left_join(checkm2)

genomad_plasmid_bin_contig = genomad_plasmid %>% mutate(bin_contig = paste0(Name, "_", seq_name))  

amr_drep = read_delim(file = "data/amr_drep.tsv", delim = "\t", col_names = T)
transposon_drep = read_delim(file = "data/transposon_drep.tsv", delim = "\t", col_names = T)

amr_drep_checkm2 = amr_drep %>% mutate(name_contig = paste0(Name, "_", X2)) %>% left_join(checkm2, by = c("Name")) %>% filter(X11 >= 90)
transposon_drep_checkm2 = transposon_drep %>% mutate(name_contig = paste0(Name, "_", X1)) %>% left_join(checkm2, by = c("Name")) %>% filter(X3 >= 90)

bat_bin_resolved_species = bat_drep %>%
  mutate(species = case_when(species ==  "no support" ~ "Novel MAGs",
                             TRUE ~ "Known MAGs")) %>%
  dplyr::select(species, Name)

bat_drep %>% left_join(checkm2) %>%
  dplyr::select(c("city_state", "Name", "Total_Contigs")) %>%
  left_join(bat_bin_resolved_species %>% dplyr::rename(species_resolved_mag = species)) %>%
  left_join(amr_drep_checkm2 %>% dplyr::select(X2, Name) %>% distinct() %>% group_by(Name) %>% summarise(AMR_contigs = n())) %>% #Number of unique contigs per bin with AMR
  left_join(amr_drep_checkm2 %>% dplyr::select(X6, Name) %>% distinct() %>% group_by(Name) %>% summarise(AMR_unique_count_per_bin = n())) %>% #Number of unique AMRs per bin
  replace(is.na(.), 0) %>%
  mutate('Proportion of AMR contigs in MAGs' = AMR_contigs/Total_Contigs) %>%
  tidyplot(x = species_resolved_mag, y = `Proportion of AMR contigs in MAGs`, color = species_resolved_mag) %>%
  add_violin() %>%
  add_data_points_beeswarm(shape = 21, alpha = 0.2, jitter_width = 0.6) %>%
  add_test_pvalue(ref.group = 1, p.adjust.method = "BH", hide_info = T, label.size = 7.5) %>%
  theme_tidyplot(fontsize = 30) %>%
  add_title("Proportion of AMR carrying contigs in MAGs") %>%
  remove_x_axis_title() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_title() %>%
  adjust_legend_title("") %>%
  split_plot(by = city_state, widths = 100, heights = 100, guides = "collect", nrow = 2)

amr_drep_checkm2_topamr_df = amr_drep_checkm2 %>% 
  dplyr::select(X6, name_contig, city_state, Name) %>% 
  distinct() %>%
  left_join(bat_bin_resolved_species)

#For each city top 10 drug classes with most number of contigs
top_amr = amr_drep_checkm2_topamr_df %>% 
  filter(!(is.na(X6))) %>%
  group_by(city_state, X6) %>%
  summarize(contig_count = n()) %>%
  ungroup() %>%
  group_by(city_state) %>%
  top_n(10,wt = contig_count)
top_amr_color = setNames(paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(as.factor(top_amr$X6)))],levels(as.factor(top_amr$X6)))
top_amr_order = top_amr %>%
  group_by(X6) %>% summarise(contig_count = sum(contig_count)) %>% 
  arrange(desc(contig_count)) %>% ungroup %>%
  mutate(X6 = as.factor(X6))

#Total contig count per city and resolved status
city_contig_count = amr_drep_checkm2 %>% dplyr::select(name_contig, city_state, Name) %>%
  left_join(bat_bin_resolved_species) %>%
  distinct() %>%
  group_by(city_state, species) %>%
  summarise(city_contig_count = n()) %>%
  ungroup()

amr_drep_checkm2_topamr_df %>% 
  filter(X6 %in% top_amr$X6) %>%
  group_by(X6, city_state, species) %>%
  summarise(contig_count_amr = n(), .groups = "drop") %>%
  left_join(city_contig_count) %>%
  mutate(contig_count_amr_prop = contig_count_amr/city_contig_count) %>%
  left_join(top_amr_order) %>%
  arrange((contig_count)) %>%
  ggplot(aes(y = fct_inorder(X6), x = contig_count_amr_prop)) +
  geom_bar( position = "identity", stat = "identity", aes(fill = species, colour = species, linetype = species, alpha = species)) +
  scale_fill_manual(values = alpha(c("white","skyblue2"), 1)) +
  scale_alpha_manual(values = c(1,0.7)) +
  scale_color_manual(values = c("skyblue2","skyblue2")) +
  scale_linetype_manual(values = c("solid", "blank")) +
  labs(fill = "", color = "", linetype = "", alpha = "",
       x = "Proportion of contigs carrying ARG", y = "") +
  facet_wrap( ~ city_state, nrow = 1) +
  theme_bw(30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Making df with bin resolved status and MDR (if more than 2 PR). Also unique contigs for each drug class, city and bin
amr_drep_checkm2_topPR_df = amr_drep_checkm2 %>% 
  dplyr::select(X14, name_contig, city_state, Name) %>% 
  distinct() %>%
  left_join(bat_bin_resolved_species) %>%
  mutate(X14 = case_when(str_count(X14, "/") >1   ~ "MDR",
                         TRUE ~ X14))

#For each city top 10 drug classes with most number of contigs
top_pr = amr_drep_checkm2_topPR_df %>% 
  filter(!(is.na(X14))) %>%
  group_by(city_state, X14) %>%
  summarize(contig_count = n()) %>%
  ungroup() %>%
  group_by(city_state) %>%
  top_n(10,wt = contig_count)
top_pr_color = setNames(paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(as.factor(top_pr$X14)))],levels(as.factor(top_pr$X14)))
top_pr_order = top_pr %>%
  group_by(X14) %>% summarise(contig_count = sum(contig_count)) %>% 
  arrange(desc(contig_count)) %>% ungroup %>%
  mutate(X14 = as.factor(X14))

#Total contig count per city and resolved status
city_contig_count = amr_drep_checkm2 %>% dplyr::select(name_contig, city_state, Name) %>%
  left_join(bat_bin_resolved_species) %>%
  distinct() %>%
  group_by(city_state, species) %>%
  summarise(city_contig_count = n()) %>%
  ungroup()

amr_drep_checkm2_topPR_df %>% 
  filter(X14 %in% top_pr$X14) %>%
  group_by(X14, city_state, species) %>%
  summarise(contig_count_pr = n(), .groups = "drop") %>%
  left_join(city_contig_count) %>%
  mutate(contig_count_pr_prop = contig_count_pr/city_contig_count) %>%
  left_join(top_pr_order) %>%
  arrange((contig_count)) %>%
  ggplot(aes(y = fct_inorder(X14), x = contig_count_pr_prop)) +
  geom_bar( position = "identity", stat = "identity", aes(fill = species, colour = species, linetype = species, alpha = species)) +
  scale_fill_manual(values = alpha(c("white","skyblue2"), 1)) +
  scale_alpha_manual(values = c(1,0.7)) +
  scale_color_manual(values = c("skyblue2","skyblue2")) +
  scale_linetype_manual(values = c("solid", "blank")) +
  labs(fill = "", color = "", linetype = "", alpha = "",
       x = "Proportion of contigs carrying drug classes", y = "") +
  facet_wrap( ~ city_state, nrow = 1) +
  theme_bw(30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

city_contig_count_mge = amr_drep_checkm2 %>% dplyr::select(name_contig, city_state, Name) %>%
  distinct() %>%
  group_by(city_state) %>%
  summarise(city_contig_count = n()) %>%
  ungroup()


#Unique Contig count per drug class per city and resolved status
drug_class_contig_count = amr_drep_checkm2 %>% 
  mutate(X14 = case_when(str_count(X14, "/") >1   ~ "MDR",
                         TRUE ~ X14)) %>%
  dplyr::select(X14, name_contig, city_state, Name) %>%
  distinct() %>% #One drug class instance per contig.
  group_by(X14, city_state) %>% 
  summarise(dc_contig_count = n()) %>%
  ungroup()

#Unique Contig count per drug class "CONTAINING MGE" per city and resolved status
drug_class_mge_contig_count = amr_drep_checkm2 %>% 
  mutate(X14 = case_when(str_count(X14, "/") >1   ~ "MDR",
                         TRUE ~ X14)) %>%
  dplyr::select(X14, name_contig, city_state, Name) %>%
  distinct() %>%
  mutate(plasmid_transposon = case_when((name_contig %in% transposon_drep_checkm2$name_contig) &  
                                          (name_contig %in% genomad_plasmid_bin_contig$bin_contig) ~ "Yes",
                                        TRUE ~ "No")) %>%
  filter(plasmid_transposon == "Yes") %>%
  group_by(X14, city_state) %>%
  summarise(mge_dc_contig_count = n()) %>%
  ungroup()

#Unique AMR count per drug class "CONTAINING MGE" per city and resolved status
drug_class_amr_count = amr_drep_checkm2 %>% 
  mutate(X14 = case_when(str_count(X14, "/") >1   ~ "MDR",
                         TRUE ~ X14)) %>%
  dplyr::select(X6, X14, city_state, Name, name_contig) %>%
  distinct() %>%
  group_by(X14, city_state) %>%
  summarise(dc_amr_count = n()) %>%
  ungroup()


drug_class_contig_count %>% left_join(drug_class_mge_contig_count) %>% left_join(drug_class_amr_count) %>% left_join(city_contig_count_mge) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  group_by(X14) %>% mutate(contig_count_per_dc_overall = sum(dc_contig_count)) %>% ungroup() %>%
  mutate(`dc_contig_per_city_contig` = dc_contig_count/city_contig_count, #Proportion of contig for each drug class out of total contig per city and resolved status 
         mge_dc_contig_count_per_dc_contig_count  = mge_dc_contig_count/dc_contig_count) %>% #Proportion of contig with mge for each dc, for each city and resolved status
  filter(X14 %in% top_pr$X14) %>%
  arrange(desc(contig_count_per_dc_overall)) %>%
  ggplot(aes(x = dc_contig_per_city_contig, y = mge_dc_contig_count_per_dc_contig_count)) +
  geom_point(aes(fill = fct_inorder(X14), size = dc_amr_count, shape = city_state) ) +
  scale_size_continuous(range = c(2,15)) +
  #scale_alpha_continuous(guide = "none", range = c(0.5,1)) +
  scale_shape_manual(values = c(21,22,23,24)) +
  theme_bw(24) +
  labs(fill = "Drug class", size = "Unique ARG \nper Drug class", shape = "City", 
       x = "Proportion of contigs with ARG", y = "Proportion of ARG contigs with MGE") +
  guides(fill = guide_legend(override.aes = list(size=5, shape = 21)),
         shape = guide_legend(override.aes = list(size = 5, fill = "grey70")),
         size = guide_legend(override.aes = list(shape = 21, fill = "grey70"))) +
  scale_fill_manual(values = paletteer::paletteer_d("ggsci::default_igv"))


contig_amr_genomad_transposon = amr_drep_checkm2 %>% 
  dplyr::select(X6, X14, name_contig, city_state, Name) %>%
  mutate(plasmid_transposon = case_when((name_contig %in% transposon_drep_checkm2$name_contig) &  
                                          (name_contig %in% genomad_plasmid_bin_contig$bin_contig) ~ "Contigs with MGE",
                                        TRUE ~ "Contigs without MGE"))

contig_amr_genomad_transposon %>% 
  group_by(X14, name_contig, city_state, plasmid_transposon) %>% 
  summarise(pr_per_count = n(), .groups = "drop") %>%
  filter(!is.na(X14)) %>%
  filter(pr_per_count > 1) %>%
  ggplot(aes(x = city_state, y = pr_per_count)) +
  geom_jitter(aes(fill = fct_inorder(X14), shape = plasmid_transposon) , alpha = 1, size = 7, width = 0.3) +
  theme_bw(35) +
  scale_shape_manual(values = c(23,21)) +
  labs(fill = "Drug class", 
       x = "", y = "Number of ARGs per contig", shape = "") +
  guides(fill = guide_legend(override.aes = list(size=5, shape = 21)),
         shape = guide_legend(override.aes = list(size = 5, fill = "grey70")),
         size = guide_legend(override.aes = list(shape = 21, fill = "grey70"))) +
  scale_fill_manual(values = paletteer::paletteer_d("ggsci::default_igv")) 

