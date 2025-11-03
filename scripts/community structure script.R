# biofilm analysis - 16S
# community structure: composition, diversity, core community, abundance vs prevalence


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

# group at phylum level
ps_phylum<-tax_glom(ps_rel, "Phylum")

# melt
melt_rel<-psmelt(ps_phylum)


# top phyla
top_phyla <- melt_rel %>% 
  group_by(Phylum) %>% 
  summarise(mean = mean(Abundance)) %>% 
  arrange(desc(mean)) %>%
  filter(mean >= 0.005) # %>% # filter by abundance threshold
# slice_head(n=10) # top 10

melt_top <- melt_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others")) %>% 
  filter(Abundance > 0) # remove zeros

phylum_abundance <- melt_top %>% 
  select(Sample, Abundance, Phylum)

write.csv(phylum_abundance, "output/phylum_abundance.csv")

mean_abundance <- phylum_abundance %>%
  group_by(Phylum) %>%
  summarise(mean = mean(Abundance),
  SD = sd(Abundance))

phylum_abundance_seasonal <- melt_top %>% 
  select(season, Sample, Abundance, Phylum) %>% 
  group_by(season, Phylum) %>% 
  summarise(mean = mean(Abundance),
            SD = sd(Abundance))

phylum_diff <- phylum_abundance_seasonal %>%
  select(season, Phylum, mean) %>%
  pivot_wider(names_from = season,
              values_from = mean,
              names_prefix = "mean_") %>%
  mutate(diff_autumn_minus_spring = mean_Autumn - mean_Spring,
         diff_spring_minus_autumn = mean_Spring - mean_Autumn)
  
top_phylum_counts <- melt_top %>%
  group_by(Sample) %>%
  slice_max(Abundance, n = 1) %>% 
  ungroup() %>%
  count(Phylum, name = "num_top")

# set colours
unique_phyla <- sort(unique(melt_top$Phylum))

# getPalette <- colorRampPalette(brewer.pal(length(unique_phyla), "Spectral")) 
# phylum_colours <- setNames(getPalette(colourCount), unique_phyla) # spectral palette

colours <- c("#E15759FF", "#FF9D9AFF", "#59A14FFF", "#8CD17DFF", "#F28E2BFF", "#FFBE7DFF",
             "#B07AA1FF", "#D4A6C8FF", "#4E79A7FF", "#A0CBE8FF","#499894FF", "grey80")
phylum_colours <- setNames(colours, unique_phyla) # manual palette

phylum_colour_map <- data.frame(Phylum = unique_phyla,
                                colour = phylum_colours,
                                stringsAsFactors = FALSE)

write.csv(phylum_colour_map, "output/phylum_colour_map.csv")

# community composition
season_labels <- melt_top %>%
  group_by(season_year) %>%
  arrange(Sample_date) %>%
  slice(1) %>%
  ungroup() %>%
  select(Sample, season_year) # get first sample of each season x year for labelling

melt_top <- melt_top[order(melt_top$Sample_date), ] # order chronologically
melt_top$Sample <- factor(melt_top$Sample, levels=unique(melt_top$Sample[order(melt_top$Sample_date)])) # reorder samples to match

composition_plot <- ggplot(melt_top, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  labs(y = "Relative abundance") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks = season_labels$Sample,
                   labels = season_labels$season_year) +
  scale_fill_manual(values = phylum_colours,
                    guide = guide_legend(ncol = 2),
                    labels = function(x) ifelse(x == "x__Others", "Others", x)) +
  geom_vline(xintercept = which(levels(melt_top$Sample) %in% season_labels$Sample), 
             linetype = "dashed", color = "white", size = 0.2)

composition_plot

legend_composition <- ggpubr::get_legend(composition_plot) # extract legend
legend_composition <- patchwork::wrap_elements(full = legend_composition) 


# diversity
diversity <- estimate_richness(ps) # calculate diversity
rownames(diversity) <- gsub("\\.","-",rownames(diversity)) # fix rownames

write.csv(diversity, "output/diversity.csv")

meta <- rownames_to_column(meta, var="16S_ID")
diversity <- rownames_to_column(diversity, var="16S_ID")
diversity_meta <- merge(meta, diversity, by = "16S_ID") # add metadata

diversity_seasonal <- diversity_meta %>% 
  select(season, Sample_ID, Shannon, Observed) %>% 
  group_by(season) %>%
  summarise(mean_diversity = mean(Shannon),
            SD_diversity = sd(Shannon),
            mean_richness = mean(Observed),
            SD_richness = sd(Observed))

diversity_meta$season <- factor(diversity_meta$season, levels = c("Spring", "Autumn")) # set order of seasons

diversity_plot <- ggplot(diversity_meta, aes(x=season, y=Shannon, fill=season))+
  geom_boxplot(outlier.size = 0.5, alpha=0.8, size = 0.2)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none") +
  # facet_wrap(~ Year) +
  labs(y="Shannon diversity", title="Diversity")+
  scale_y_continuous(expand=c(0,0),limits=c(0,7),breaks=seq(0,7,by=2)) +
  scale_fill_manual(values=c("forestgreen", "steelblue"))

diversity_plot

richness_plot <- ggplot(diversity_meta, aes(x=season, y=Observed, fill=season))+
  geom_boxplot(outlier.size = 0.5, alpha=0.8, size = 0.2)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none") +
  # facet_wrap(~ Year) +
  labs(y="ASV richness", title="Richness")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1500),breaks=seq(0,1500,by=500)) +
  scale_fill_manual(values=c("forestgreen", "steelblue"))

richness_plot

diversity_richness_plot <- diversity_plot + richness_plot +
  plot_layout(axes = "collect", ncol=2) &
  theme(plot.title = element_blank())

diversity_richness_plot

# diversity and richness stats
diversity_meta$Year <- factor(diversity_meta$Year)
diversity_meta$season <- factor(diversity_meta$season)

anova_model <- aov(Shannon ~ Year * season, data = diversity_meta)
summary(anova_model)

anova_model <- aov(Observed ~ Year * season, data = diversity_meta)
summary(anova_model)


# core community
ps_core <- core(ps_rel, detection = 0.01, prevalence = 0.5) # core if relative abundance >detection in >prevalence % of samples
core_tax <- as.data.frame(tax_table(ps_core)) # extract core taxa
core_rel <- t(as.data.frame(otu_table(ps_core))) # extract core relative abundance

phylum_counts <- core_tax %>%
  count(Phylum, name = "count")

ps_core <- core(ps_rel, detection = 0.0001, prevalence = 0.5) # core if relative abundance >detection in >prevalence % of samples
core_tax <- as.data.frame(tax_table(ps_core)) # extract core taxa
core_rel <- t(as.data.frame(otu_table(ps_core))) # extract core relative abundance

phylum_counts <- core_tax %>%
  count(Phylum, name = "count")

detections <- c(0.0001, 0.001, 0.01)
prevalences <- seq(0.05, 1, by = 0.05) # 0.5 = 50% prevalence

core_df <- lapply(detections, function(thresh) {
  prevalence <- rowMeans(core_rel >= thresh)
  data.frame(ASV = rownames(core_rel),
             detection = thresh,
             prevalence = prevalence)
}) %>%
  bind_rows() %>% 
  merge(., core_tax, by.x = "ASV", by.y = "row.names")

write.csv(core_df, "output/core.csv")

ASV_order <- core_df %>%
  group_by(ASV) %>%
  summarise(mean_prev = mean(prevalence)) %>%
  arrange(mean_prev) %>%
  pull(ASV)

core_df$ASV <- factor(core_df$ASV, levels = ASV_order)

core_plot <- ggplot(core_df, aes(x = detection, y = ASV, fill = prevalence)) +
  geom_tile() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.8, "lines"),
        legend.position = "bottom",
        axis.text = element_text(size = 4.5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black")) +
  labs(x = "Detection threshold (relative abundance)") +
  scale_x_log10(breaks = detections) +
  scale_fill_gradientn(colours = rev(brewer.pal(9, "RdYlBu")), 
                       name = "Prevalence (%)",
                       limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                       labels = function(x) formatC(x * 100, format = "f", digits = 0)) +
  scale_y_discrete(labels = function(x) parse(text = core_df$taxon_label[match(x, core_df$taxon)]))

core_plot

legend_core <- ggpubr::get_legend(core_plot) # extract legend
legend_core <- patchwork::wrap_elements(full = legend_core) 


# abundance vs prevalence
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>% # convert to relative abundance
  t(.) # samples as rows

mean_abundance <- as.data.frame(colMeans(ASV_rel)) # calculate mean across all samples
mean_abundance <- mean_abundance %>% rownames_to_column(., var = "ASV")
colnames(mean_abundance)[2] <- "mean_rel_abund"

presence_absence <- (ASV_rel > 0) * 1 # convert coverage to binary presence/absence
occupancy <- as.data.frame(colSums(presence_absence)) # sum number of samples each ASV had abundance > 0 in
occupancy <- occupancy %>% rownames_to_column(., var = "ASV")
colnames(occupancy)[2] <- "occupancy"

total_samples <- nrow(ASV_rel) # calculate number of samples

occupancy$percentage_occupancy <- (occupancy$occupancy / total_samples) * 100 # convert to percentage of total samples

occupancy_abundance <- merge(occupancy, mean_abundance, by="ASV")
occupancy_tax <- merge(occupancy_abundance, tax_tab2, by.x="ASV", by.y="row.names")

scientific_10_labels <- function(x) {
  parse(text = paste0("10^", floor(log10(abs(x))))) # convert labels to scientific notation
}

occupancy_tax2 <- occupancy_tax %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others"))

occupancy_tax2$core <- ifelse(occupancy_tax2$ASV %in% core_tax$taxon, "Core", "Not core")

abundance_occupancy_log <- ggplot(occupancy_tax2, aes(x = percentage_occupancy, y = mean_rel_abund, colour = Phylum)) +
  geom_point(size = 0.8, alpha = 0.8, aes(shape = core)) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_blank(), # element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 4)) +
  labs(x = "Prevalence (%)", y = "Mean relative abundance") +
  scale_y_log10(labels = scientific_10_labels) +
  scale_colour_manual(values = phylum_colours,
                      labels = function(x) ifelse(x == "x__Others", "Others", x)) +
  scale_shape_manual(name = "Core ASV",
                     guide = guide_legend(ncol = 2),
                     values = c("Not core" = 16, "Core" = 17)) +
  geom_segment(aes(x = 50, xend = 50, y = 0.0001, yend = Inf),
               linetype = "dashed", color = "black", size = 0.2) +
  geom_segment(aes(x = 50, xend = Inf, y = 0.0001, yend = 0.0001),
               linetype = "dashed", color = "black", size = 0.2)
  
abundance_occupancy_log

ggsave("output/figures/ASV_abund_prev.png", abundance_occupancy_log, width = 120, height = 90, units = "mm", dpi = 600)
ggsave("output/figures/ASV_abund_prev.pdf", abundance_occupancy_log, width = 120, height = 90, units = "mm", dpi = 600)

legend_shape_plot <- abundance_occupancy_log +
  guides(colour = "none")  # remove colour legend
legend_occupancy <- ggpubr::get_legend(legend_shape_plot) # extract shape legend
legend_occupancy <- patchwork::wrap_elements(full = legend_occupancy)


# combined plot
layout <- "
AA
BE
CE
DE
"

legends <- legend_composition + (legend_occupancy / legend_core)

combined_plot <- composition_plot + abundance_occupancy_log + diversity_richness_plot +legends + core_plot +
  plot_layout(design = layout, heights = c(1.5, 0.8, 0.8, 0.5, 1)) &
  theme(legend.position = "none")

combined_plot

ggsave("output/figures/community.png", combined_plot, width = 180, height = 200, units = "mm", dpi = 800)
ggsave("output/figures/community.pdf", combined_plot, width = 180, height = 200, units = "mm", dpi = 400)

