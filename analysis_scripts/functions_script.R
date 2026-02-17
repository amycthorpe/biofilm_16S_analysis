# biofilm analysis - 16S
# predicted functions


# load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(tidyverse)
library(glue)
library(phyloseq)
library(vegan)
library(microbiome)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(microeco)
library(file2meco)


set.seed(123)


# import data
ASV_tab <- data.frame(fread("data/ASV_counts_16S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table("data/ASV_tax_16S_processed.csv", header=T, row.names=1, check.names=F, sep=",")
meta_all <- read.table("data/seq_metadata_all_amplicons.csv", header=T, check.names=F, sep=",")

meta <- meta_all[meta_all$`16S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="16S_ID") # make seq ID rownames

meta <- meta %>% # add season column
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date

meta <- meta %>%
  mutate(season_year = paste(season, Year)) # get season x year col


# label unassigned genera
first_assigned_label <- function(Family, Order, Class, Phylum) {
  levels <- list(Family, Order, Class, Phylum)
  for (lvl in levels) {
    if (!is.na(lvl) && lvl != "unassigned") {
      return(lvl)}}
  return(NA)}

tax_tab2 <- tax_tab
tax_tab2$taxon <- rownames(tax_tab2)

tax_tab2 <- tax_tab2 %>%
  mutate(new_label = if_else(Genus != "unassigned", Genus,
                             paste("unassigned", pmap_chr(list(Family, Order, Class, Phylum), first_assigned_label))),
         taxon_label = if_else(Genus != "unassigned",
                               glue("'{taxon}: '*italic('{new_label}')"), # make genus italic (but not 'unassigned')
                               glue("'{taxon}: {new_label}'")))


# phyloseq object
TAX <- tax_table(as.matrix(tax_tab2))
ASV <- otu_table(ASV_tab, taxa_are_rows = F)
samples <- sample_data(meta)
ps <- phyloseq(ASV, TAX, samples)
ps

# convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))


# functions
meco <- phyloseq2meco(ps_rel) # phyloseq to microeco object

t1 <- trans_func$new(dataset = meco)
t1$cal_spe_func(prok_database = "FAPROTAX") # get ASV functions
functions_ASVs <- t1$res_spe_func # 1 represent trait exists, 0 represent no or cannot confirm

write.csv(functions_ASVs, "output/functions_ASVs.csv")

colsums <- colSums(functions_ASVs)
functions_to_keep <- names(colsums[colsums != 0]) # get list of functions detected
functions_ASVs <- functions_ASVs[, functions_to_keep] # keep only detected functions


# add indicators, modules and core community
titan <- read.csv("output/df_filtered.csv", row.names=1) %>% 
  select(taxon, variable, maxgrp) %>% 
  distinct(taxon, variable, maxgrp, .keep_all = TRUE) %>% 
  mutate(maxgrp = ifelse(maxgrp == 2, "Tolerant", ifelse(maxgrp == 1, "Sensitive", as.character(maxgrp)))) %>% 
  mutate(variable = recode(variable,
                           "alkalinity" = "Alkalinity",
                           "nitrate" = "Nitrate-N",
                           "temp" = "Temperature",
                           "DO" = "DO"))

modules <- read.csv("output/modules.csv", header = T, check.names = F, row.names = 1) %>% 
  select(name, taxa_roles, module)

core_df <- read.table("output/core.csv", header=T, row.names=1, check.names=F, sep=",") 

functions_merged <- merge(functions_ASVs, modules, by.x="row.names", by.y="name", all.x=T)
colnames(functions_merged)[1] <- "taxon"

functions_merged <- merge(functions_merged, titan, by="taxon", all.x=T)

functions_merged$core <- ifelse(functions_merged$taxon %in% core_df$taxon, "Core", "Not core")

functions_filtered <- functions_merged %>%
  filter(grepl("Core", core) | # keep only ASVs that are either indicators, core or are module connectors or hubs
           taxa_roles %in% c("Connectors", "Module hubs") |
           !is.na(variable))


# with abundance and occupancy
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows

abundance <- ASV_rel %>% 
  as.data.frame() %>% 
  summarise(across(everything(), mean)) %>% # calculate mean relative abundance
  t() %>% 
  as.data.frame()

colnames(abundance) <- "mean_rel_abund"

occupancy <- ASV_rel %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ sum(. != 0))) %>% # calculate occupancy (count of samples where abundance > 0)
  t() %>%
  as.data.frame()

colnames(occupancy) <- "occupancy"

n_samples <- nrow(ASV_rel)
occupancy$percent_occupancy <- (occupancy$occupancy / n_samples) * 100 # calculate % occupancy

abund_occ <- merge(abundance, occupancy, by = "row.names")

colnames(abund_occ)[1] <- "taxon"

functions_filtered <- merge(functions_filtered, abund_occ, by="taxon", all.x=T)

# get functional categories
func_groups <- t1$.__enclos_env__$private$default_prok_func_group$FAPROTAX # extract categories
mapping <- stack(func_groups)

colnames(mapping)[1] <- "Function"
colnames(mapping)[2] <- "category"

functions_long <- functions_filtered %>%
  pivot_longer(cols = -c(taxon, taxa_roles, module, core, variable, maxgrp, mean_rel_abund, occupancy, percent_occupancy), names_to = "Function", values_to = "Presence") %>%
  filter(Presence > 0) %>%
  left_join(mapping, by="Function")
  
missing_functions <- functions_long %>%
  distinct(Function) %>%
  anti_join(mapping, by = "Function")

missing_functions <- missing_functions %>%
  mutate(category = case_when(
    Function %in% c(
      "photosynthetic_cyanobacteria",
      "oxygenic_photoautotrophy",
      "phototrophy",
      "chemoheterotrophy",
      "chlorate_reducers",
      "dark_hydrogen_oxidation") ~ "Energy source",
    Function %in% c(
      "human_pathogens_pneumonia",
      "human_pathogens_all",
      "human_associated",
      "animal_parasites_or_symbionts") ~ "Host-associated",
    Function == "methanol_oxidation" ~ "C-cycle",
    TRUE ~ NA_character_))

mapping <- rbind(mapping, missing_functions)


# plot
functions_long <- functions_filtered %>%
  pivot_longer(cols = -c(taxon, taxa_roles, module, core, maxgrp, variable, mean_rel_abund, occupancy, percent_occupancy), names_to = "Function", values_to = "Presence") %>%
  filter(Presence > 0) %>%
  left_join(mapping, by="Function") %>%
  merge(tax_tab2, by="taxon") %>%
  mutate(Function = factor(Function, levels = unique(Function))) %>% 
  mutate(Function = Function %>%
           str_replace_all("_", " ") %>% # reformat text for plotting
           str_to_sentence()) %>%
  mutate(fill_colour = case_when( 
    Presence == 1 & category == "Energy source" ~ "#9E9AC8", # add colours for plot
    Presence == 1 & category == "C-cycle" ~ "#6BAED6",
    Presence == 1 & category == "N-cycle" ~ "#74C476",
    Presence == 1 & category == "Host-associated"~ "#FD8D3C",
    TRUE ~ "white")) %>% 
  filter(!Function %in% c("Chemoheterotrophy",
                          "Photosynthetic cyanobacteria", "Phototrophy", "Photoautotrophy",
                          "Human pathogens all", "Human pathogens pneumonia")) # remove redundant rows

all_taxa <- unique(functions_long$taxon)
functions_long$taxon <- factor(functions_long$taxon, levels = all_taxa)
functions_long <- functions_long %>% 
  mutate(taxon = fct_reorder(taxon, mean_rel_abund)) # set order

plot <- ggplot(functions_long, aes(x = Function, y = taxon, fill = fill_colour)) +
  geom_tile(colour = "white") +
  # scale_fill_manual(values = c("0" = "white", "1" = "steelblue"), name = "Presence") +
  scale_fill_identity() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black")) +
  scale_y_discrete(labels = function(x) parse(text = functions_long$taxon_label[match(x, functions_long$taxon)])) +
  facet_grid(. ~ category, scales = "free_x", space = "free_x")

plot


# modules
module_plot <- ggplot(functions_long, aes(x = "Network role", y = taxon, fill = taxa_roles)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_manual(values = c("Connectors" = "palevioletred",
                               "Module hubs" = "forestgreen",
                               "Peripheral nodes" = "steelblue"),
                    na.value = "white", na.translate = FALSE,
                    guide = guide_legend(ncol = 1)) +
  labs(fill = "Network role", title = "Network role") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = functions_long$taxon_label[match(x, functions_long$taxon)]))

module_plot

legend_module <- patchwork::wrap_elements(full = cowplot::get_legend(module_plot)) # extract legend


# taxonomy
phylum_colour_map <- read.table("output/phylum_colour_map.csv", check.names=F, row.names=1, header=T, sep = ",")
phylum_colours <- setNames(phylum_colour_map$colour, phylum_colour_map$Phylum)

tax_plot <- ggplot(functions_long, aes(x = "Phylum", y = taxon, fill = Phylum)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  labs(fill = "Phylum", title = "Phylum") +
  scale_x_discrete(position = "top") +
  scale_fill_manual(guide = guide_legend(nrow = 2),
                    values = phylum_colours) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = functions_long$taxon_label[match(x, functions_long$taxon)]))

tax_plot

legend_tax <- patchwork::wrap_elements(full = cowplot::get_legend(tax_plot)) # extract legend


# core
core_plot <- ggplot(functions_long, aes(x = "Core", y = taxon, colour = core)) +
  geom_point(size = 4, shape = 42) + 
  scale_colour_manual(values = c("Core" = "black", "Not core" = "white")) +
  # geom_point(data = functions_long %>% filter(core == "Core"), shape = 42, size = 4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  labs(fill = "Core", title = "Core") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = functions_long$taxon_label[match(x, functions_long$taxon)]))

core_plot


# titan indicators
plot_data <- functions_long %>%
  tidyr::complete(taxon, variable, fill = list(maxgrp = NA, taxon_label = NA))

titan_plot <- ggplot(plot_data %>% filter(!is.na(variable)),
                     aes(x = variable, y = taxon, fill = maxgrp)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_manual(values = c("Tolerant" = "#eb4a40",
                               "Sensitive" = "#2a5674"),
                    na.value = "white", na.translate = FALSE,
                    guide = guide_legend(ncol = 1)) +
  labs(fill = "TITAN indicator", title = "TITAN indicator") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = functions_long$taxon_label[match(x, functions_long$taxon)]))

titan_plot

legend_titan <- patchwork::wrap_elements(full = cowplot::get_legend(titan_plot)) # extract legend


# combined plot
plots <- core_plot + tax_plot + module_plot + titan_plot + plot +
  plot_layout(ncol = 5, widths = c(1, 1, 1, 4, 30), guides = "collect", axes = "collect")

plots

legends <- legend_tax + legend_module + legend_titan +
  plot_layout(widths=c(1,1,1))

legends

combined <- plots / legends +
  plot_layout(guides = "collect", heights = c(6, 0.8)) &
  theme(legend.position = "none",
        plot.margin = margin(t = -1, r = 1, b = -1, l = 1))

combined

ggsave("output/figures/functions.pdf", combined, width = 180, height = 190, units = "mm", dpi = 600)
ggsave("output/figures/functions.png", combined, width = 180, height = 190, units = "mm", dpi = 600)

