# biofilm analysis for manuscripts - 16S
# TITAN analysis


# load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(TITAN2)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(patchwork)
library(VennDiagram)
library(ggvenn)
library(glue)


set.seed(123)


# import data
ASV_tab <- data.frame(fread("data/ASV_counts_16S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table("data/ASV_tax_16S_processed.csv", header=T, row.names=1, check.names=F, sep=",")
meta_all <- read.table("data/seq_metadata_all_amplicons.csv", header=T, check.names=F, sep=",")

meta <- meta_all[meta_all$`16S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="16S_ID") # make seq ID rownames

meta <- meta %>% # add season columm
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date


# convert to relative abundance
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows


# environmental metadata
meta_merged <- read.csv("data/meta_merged.csv", header = T, row.names = 1, check.names = F) %>% 
  select(Alkalinity_combined, Nitrate_N, DO, Water_temp)

meta_merged <- as.data.frame(na.omit(meta_merged)) # remove rows with NA

variables <- "Alkalinity_combined" # subset metadata to numeric variables of interest
meta_sub <- meta_merged[variables]

ASV_sub <- as.data.frame(ASV_rel[rownames(ASV_rel) %in% rownames(meta_sub), ]) # subset to remaining samples
meta_sub <- as.data.frame(meta_sub[intersect(rownames(meta_sub), rownames(ASV_rel)), , drop=F]) # subset to remaining samples
ASV_sub <- ASV_sub[match(rownames(meta_sub), rownames(ASV_sub)), ] # match order of rows


# filter
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 78] # increase minimum occurrence frequency to present in >5% sites to reduce dataset size and run time - 24067 to 2651 ASVs
ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 156] # increase minimum occurrence frequency to present in >10% sites to reduce dataset size and run time - 1287 ASVs
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 312] # increase minimum occurrence frequency to present in >20% sites to reduce dataset size and run time - 526 ASVs
# ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= 467] # increase minimum occurrence frequency to present in >30% sites to reduce dataset size and run time - 278 ASVs
# ASV_sub4 <- ASV_sub4[, colSums(ASV_sub4 > 0) >= 3] # minimum occurrence frequency for titan is 3


# alkalinity
# prep data
titan_alk <- readRDS("output/titan/titan_alkalinity_10perc.rds")
titan_alk$sumz.cp # check change points

tax_tab$taxonomy <- paste0( # consolidate taxonomy into one column
  "P__", tax_tab$Phylum, "; ",
  "C__", tax_tab$Class, "; ",
  "O__", tax_tab$Order, "; ",
  "F__", tax_tab$Family, "; ",
  "G__", tax_tab$Genus, "; ",
  rownames(tax_tab))


# plot with titan functions
density <- plot_sumz_density(titan_alk,
                             axis.text.y = 7,
                             xlabel="Alkalinity")
density

taxa <- plot_taxa_ridges(titan_alk,
                         grid = F,
                         n_ytaxa = 250,
                         axis.text.y = 7,
                         xlabel="Alkalinity")
taxa


# plot with.ggplot
sppmax <- as_tibble(titan_alk$sppmax, rownames = "taxon") %>% # extract summary table and clean colnames
  rename(Q5 = `5%`,
         Q10 = `10%`,
         Q50 = `50%`,
         Q90 = `90%`,
         Q95 = `95%`)

sppmax$variable <- "alkalinity"

sppmax_alk <- sppmax

n_taxa <- dim(titan_alk$metricArray)[1] # get taxa
n_boot <- dim(titan_alk$metricArray)[3] # get boot

chk_pts_df <- map_dfr(1:n_taxa, function(i) {
  tibble(
    taxon = sppmax$taxon[i],
    chk_pts = titan_alk$metricArray[i, 2, ],  # metric 2 = change points
    bootstrap = 1:n_boot)})

df_alk <- left_join(chk_pts_df, sppmax, by = "taxon") # merge

df_filtered1 <- df_alk %>%
  merge(., tax_tab, by.x = "taxon", by.y = "row.names") %>% 
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= 0.95, reliability >= 0.95) # filter to only pure and reliable taxa

threshold <- quantile(sppmax$z.median, 0.90) # calculate cut off threshold for z score - ASVs in the 90th percentile

df_filtered <- df_filtered1 %>% 
  filter(z.median > threshold) # filter based on z score median (strength and consistency of response to gradient)

df_increasers <- df_filtered %>% 
  filter(filter == 2) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                          distinct(taxonomy, Q50) %>%
                          arrange(desc(Q50)) %>%  # for increasers: high to low
                          pull(taxonomy)))

df_decreasers <- df_filtered %>% 
  filter(filter == 1) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                          distinct(taxonomy, Q50) %>%
                          arrange(Q50) %>%  # for decreases: low to high
                          pull(taxonomy)))

increasers <- ggplot(df_increasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40") + #, limits = c(30, 65), breaks=seq(30, 60, by = 10)) +
  labs(x = expression("Change point (Alkalinity, mg" ~ L^{-1} * " of CaCO"[3]*")"), fill = "Z score") +
  # geom_segment(data = df_increasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  # scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  geom_vline(xintercept = titan_alk$sumz.cp["fsumz+", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

increasers

decreasers <- ggplot(df_decreasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank()) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674") + #, limits = c(30, 65), breaks=seq(30, 60, by = 10)) +
  labs(x = expression("Change point (Alkalinity, mg" ~ L^{-1} * ")"), fill = "Z score") +
  # geom_segment(data = df_decreasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  # scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  geom_vline(xintercept = titan_alk$sumz.cp["fsumz-", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

decreasers

combined <- wrap_plots(decreasers, increasers, ncol=1, heights=c(1,3), axes="collect")
combined

phylum_cp <- df_filtered1 %>% 
  group_by(filter, Phylum) %>% 
  summarise(mean_zenv.cp = mean(zenv.cp)) # phylum level change points


# nitrate-N
# prep data
titan_nitrate <- readRDS("output/titan/titan_nitrate_10perc.rds")
titan_nitrate$sumz.cp # check change points


# plot with titan functions
density <- plot_sumz_density(titan_nitrate,
                             axis.text.y = 7,
                             xlabel="Nitrate-N")
density

taxa <- plot_taxa_ridges(titan_nitrate,
                         grid = F,
                         n_ytaxa = 250,
                         axis.text.y = 7,
                         xlabel="Nitrate-N")
taxa


# plot with.ggplot
sppmax <- as_tibble(titan_nitrate$sppmax, rownames = "taxon") %>% # extract summary table and clean colnames
  rename(Q5 = `5%`,
         Q10 = `10%`,
         Q50 = `50%`,
         Q90 = `90%`,
         Q95 = `95%`)

sppmax$variable <- "nitrate"

sppmax_nitrate <- sppmax

n_taxa <- dim(titan_nitrate$metricArray)[1] # get taxa
n_boot <- dim(titan_nitrate$metricArray)[3] # get boot

chk_pts_df <- map_dfr(1:n_taxa, function(i) {
  tibble(
    taxon = sppmax$taxon[i],
    chk_pts = titan_nitrate$metricArray[i, 2, ],  # metric 2 = change points
    bootstrap = 1:n_boot)})

df_nitrate <- left_join(chk_pts_df, sppmax, by = "taxon") # merge

df_filtered1 <- df_nitrate %>%
  merge(., tax_tab, by.x = "taxon", by.y = "row.names") %>% 
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= 0.95, reliability >= 0.95) # filter to only pure and reliable taxa

threshold <- quantile(sppmax$z.median, 0.90) # calculate cut off threshold for z score - ASVs in the 90th percentile

df_filtered <- df_filtered1 %>% 
  filter(z.median > threshold) # filter based on z score median (strength and consistency of response to gradient)

df_increasers <- df_filtered %>% 
  filter(filter == 2) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(desc(Q50)) %>%  # for increasers: high to low
                             pull(taxonomy)))

df_decreasers <- df_filtered %>% 
  filter(filter == 1) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(Q50) %>%  # for decreases: low to high
                             pull(taxonomy)))

increasers <- ggplot(df_increasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40", limits = c(19, 35), breaks=seq(20, 35, by = 5)) +
  labs(x = expression("Change point (Nitrate-N, mg" ~ L^{-1} * ")"), fill = "Z score") +
  # geom_segment(data = df_increasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1)) +
  geom_vline(xintercept = titan_nitrate$sumz.cp["fsumz+", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

increasers

decreasers <- ggplot(df_decreasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank()) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674", limits = c(19, 35), breaks=seq(20, 35, by = 5)) +
  labs(x = expression("Change point (Nitrate-N, mg" ~ L^{-1} * ")"), fill = "Z score") +
  # geom_segment(data = df_decreasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1)) +
  geom_vline(xintercept = titan_nitrate$sumz.cp["fsumz-", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

decreasers

combined <- wrap_plots(decreasers, increasers, ncol=1, heights=c(1,4), axes="collect")
combined

phylum_cp <- df_filtered1 %>% 
  group_by(filter, Phylum) %>% 
  summarise(mean_zenv.cp = mean(zenv.cp)) # phylum level change points


# DO
# prep data
titan_DO <- readRDS("output/titan/titan_DO_10perc.rds")
titan_DO$sumz.cp # check change points


# plot with titan functions
density <- plot_sumz_density(titan_DO,
                             axis.text.y = 7,
                             xlabel="Nitrate-N")
density

taxa <- plot_taxa_ridges(titan_DO,
                         grid = F,
                         n_ytaxa = 250,
                         axis.text.y = 7,
                         xlabel="DO")
taxa


# plot with.ggplot
sppmax <- as_tibble(titan_DO$sppmax, rownames = "taxon") %>% # extract summary table and clean colnames
  rename(Q5 = `5%`,
         Q10 = `10%`,
         Q50 = `50%`,
         Q90 = `90%`,
         Q95 = `95%`)

sppmax$variable <- "DO"

sppmax_DO <- sppmax

n_taxa <- dim(titan_DO$metricArray)[1] # get taxa
n_boot <- dim(titan_DO$metricArray)[3] # get boot

chk_pts_df <- map_dfr(1:n_taxa, function(i) {
  tibble(
    taxon = sppmax$taxon[i],
    chk_pts = titan_DO$metricArray[i, 2, ],  # metric 2 = change points
    bootstrap = 1:n_boot)})

df_DO <- left_join(chk_pts_df, sppmax, by = "taxon") # merge

df_filtered1 <- df_DO %>%
  merge(., tax_tab, by.x = "taxon", by.y = "row.names") %>% 
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= 0.95, reliability >= 0.95) # filter to only pure and reliable taxa

threshold <- quantile(sppmax$z.median, 0.90) # calculate cut off threshold for z score - ASVs in the 90th percentile

df_filtered <- df_filtered1 %>% 
  filter(z.median > threshold) # filter based on z score median (strength and consistency of response to gradient)

df_increasers <- df_filtered %>% 
  filter(filter == 2) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(desc(Q50)) %>%  # for increasers: high to low
                             pull(taxonomy)))

df_decreasers <- df_filtered %>% 
  filter(filter == 1) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(Q50) %>%  # for decreases: low to high
                             pull(taxonomy)))

increasers <- ggplot(df_increasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40", limits = c(19, 40), breaks=seq(20, 40, by = 5)) +
  labs(x = expression("Change point (DO, mg" ~ L^{-1} * ")"), fill = "Z score") +
  # geom_segment(data = df_increasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(8, 14), breaks = seq(8, 14, by = 1)) +
  geom_vline(xintercept = titan_DO$sumz.cp["fsumz+", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

increasers

decreasers <- ggplot(df_decreasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank()) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674", limits = c(19, 40), breaks=seq(20, 40, by = 5)) +
  labs(x = expression("Change point (DO, mg" ~ L^{-1} * ")"), fill = "Z score") +
  # geom_segment(data = df_decreasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(8, 14), breaks = seq(8, 14, by = 1)) +
  geom_vline(xintercept = titan_DO$sumz.cp["fsumz-", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

decreasers

combined <- wrap_plots(decreasers, increasers, ncol=1, heights=c(2,1), axes="collect")
combined

phylum_cp <- df_filtered1 %>% 
  group_by(filter, Phylum) %>% 
  summarise(mean_zenv.cp = mean(zenv.cp)) # phylum level change points


# temp
# prep data
titan_temp <- readRDS("output/titan/titan_temp_10perc.rds")
titan_temp$sumz.cp # check change points


# plot with titan functions
density <- plot_sumz_density(titan_temp,
                             axis.text.y = 7,
                             xlabel="Nitrate-N")
density

taxa <- plot_taxa_ridges(titan_temp,
                         grid = F,
                         n_ytaxa = 250,
                         axis.text.y = 7,
                         xlabel="DO")
taxa


# plot with.ggplot
sppmax <- as_tibble(titan_temp$sppmax, rownames = "taxon") %>% # extract summary table and clean colnames
  rename(Q5 = `5%`,
         Q10 = `10%`,
         Q50 = `50%`,
         Q90 = `90%`,
         Q95 = `95%`)

sppmax$variable <- "temp"

sppmax_temp <- sppmax

n_taxa <- dim(titan_temp$metricArray)[1] # get taxa
n_boot <- dim(titan_temp$metricArray)[3] # get boot

chk_pts_df <- map_dfr(1:n_taxa, function(i) {
  tibble(
    taxon = sppmax$taxon[i],
    chk_pts = titan_temp$metricArray[i, 2, ],  # metric 2 = change points
    bootstrap = 1:n_boot)})

df_temp <- left_join(chk_pts_df, sppmax, by = "taxon") # merge

df_filtered1 <- df_temp %>%
  merge(., tax_tab, by.x = "taxon", by.y = "row.names") %>% 
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= 0.95, reliability >= 0.95) # filter to only pure and reliable taxa

threshold <- quantile(sppmax$z.median, 0.90) # calculate cut off threshold for z score - ASVs in the 90th percentile

df_filtered <- df_filtered1 %>% 
  filter(z.median > threshold) # filter based on z score median (strength and consistency of response to gradient)

df_increasers <- df_filtered %>% 
  filter(filter == 2) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(desc(Q50)) %>%  # for increasers: high to low
                             pull(taxonomy)))

df_decreasers <- df_filtered %>% 
  filter(filter == 1) %>%
  mutate(taxonomy = factor(taxonomy, levels = df_filtered %>%
                             distinct(taxonomy, Q50) %>%
                             arrange(Q50) %>%  # for decreases: low to high
                             pull(taxonomy)))

increasers <- ggplot(df_increasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40", limits = c(18, 41), breaks=seq(20, 40, by = 5)) +
  labs(x = expression("Change point (Temperature, °C)"), fill = "Z score") +
  # geom_segment(data = df_increasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(6, 18), breaks = seq(6, 18, by = 2)) +
  geom_vline(xintercept = titan_temp$sumz.cp["fsumz+", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

increasers

decreasers <- ggplot(df_decreasers, aes(x = chk_pts, y = fct_rev(taxonomy), fill = zscore)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank()) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674", limits = c(18, 41), breaks=seq(20, 40, by = 5)) +
  labs(x = expression("Change point (Temperature, °C)"), fill = "Z score") +
  # geom_segment(data = df_decreasers %>% distinct(taxon, Q50),
  #              aes(x = Q50, xend = Q50,
  #                  y = as.numeric(fct_rev(factor(taxon))), yend = as.numeric(fct_rev(factor(taxon))) + 0.9),
  #              inherit.aes = F, colour = "black", size = 0.3) +
  scale_x_continuous(limits = c(6, 18), breaks = seq(6, 18, by = 2)) +
  geom_vline(xintercept = titan_temp$sumz.cp["fsumz-", "cp"], linetype = "dashed", color = "grey50", size = 0.5)

decreasers

combined <- wrap_plots(decreasers, increasers, ncol=1, heights=c(1,1.5), axes="collect")
combined

phylum_cp <- df_filtered1 %>% 
  group_by(filter, Phylum) %>% 
  summarise(mean_zenv.cp = mean(zenv.cp)) # phylum level change points


# combined titan results
sppmax_all <- rbind(sppmax_alk, sppmax_nitrate, sppmax_DO, sppmax_temp)
sppmax_all <- sppmax_all %>% 
  filter(purity >= 0.95, reliability >= 0.95) 

write.csv(sppmax_all, "output/titan_all.csv")

results_df <- rbind(df_alk, df_nitrate, df_DO, df_temp)

results_filtered <- results_df %>% 
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= 0.95, reliability >= 0.95) # filter to only pure and reliable taxa

# venn diagram
alk_taxa <- unique(results_filtered$taxon[results_filtered$variable == "alkalinity"])
nitrate_taxa <- unique(results_filtered$taxon[results_filtered$variable == "nitrate"])
temp_taxa <- unique(results_filtered$taxon[results_filtered$variable == "temp"])
DO_taxa <- unique(results_filtered$taxon[results_filtered$variable == "DO"])

taxon_list <- list(Alkalinity = alk_taxa,
                   'Nitrate-N' = nitrate_taxa,
                   Temperature = temp_taxa,
                   DO = DO_taxa)

plot <- ggvenn(taxon_list, fill_color = c("skyblue", "lightpink", "lightgreen", "orange"),
               stroke_size = 0.4, set_name_size = 4, text_size = 2.8)

plot

ggsave("output/figures/titan_venn.png", plot, width = 120, height = 90, units = "mm", dpi = 600)
ggsave("output/figures/titan_venn.pdf", plot, width = 120, height = 90, units = "mm", dpi = 600)


# with modules
modules <- read.csv("output/modules.csv", header = T, check.names = F, row.names = 1) %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus)

titan_modules <- merge(results_filtered, modules, by.x="taxon", by.y="name", all.x = T) # merge titan and modules
sppmax_modules <- merge(sppmax_all, modules, by.x="taxon", by.y="name", all.x = T)

# with abundance and occupancy
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

titan_modules_abund <- merge(titan_modules, abund_occ, by="taxon", all.x=T) # merge titan, modules, abundance and occupancy
sppmax_modules_abund <- merge(sppmax_modules, abund_occ, by="taxon", all.x=T)

# with taxonomy
final_df <- merge(titan_modules_abund, tax_tab, by.x="taxon", by.y="row.names", all.x=T)
sppmax_modules_abund_tax <- merge(sppmax_modules_abund, tax_tab, by.x="taxon", by.y="row.names", all.x=T)

write.csv(final_df, "output/titan_modules_abund.csv")


# with core community
core_df <- read.table("output/core.csv", header=T, row.names=1, check.names=F, sep=",") 

final_df$core <- ifelse(final_df$taxon %in% core_df$taxon, "Core", "Not core")
sppmax_modules_abund_tax$core <- ifelse(sppmax_modules_abund_tax$taxon %in% core_df$taxon, "Core", "Not core")
modules_core <- modules
modules_core$core <- ifelse(modules_core$name %in% core_df$taxon, "Core", "Not core")

# compiled plots

# names for y axis - append genus to ASV ID (if unassigned, take assignment at deepest rank after 'unassigned')
first_assigned_label <- function(Family, Order, Class, Phylum) {
  levels <- list(Family, Order, Class, Phylum)
  for (lvl in levels) {
    if (!is.na(lvl) && lvl != "unassigned") {
      return(lvl)}}
  return(NA)}

final_df <- final_df %>%
  mutate(new_label = if_else(Genus != "unassigned", Genus,
                             paste("unassigned", pmap_chr(list(Family, Order, Class, Phylum), first_assigned_label))),
         taxon_label = if_else(Genus != "unassigned",
                               glue("'{taxon}: '*italic('{new_label}')"), # make genus italic (but not 'unassigned')
                               glue("'{taxon}: {new_label}'")))

# filter ASVs
df_filtered <- final_df %>% 
  filter(z.median > 15) %>% # filter based on z median (magnitude of response)
  # filter(!is.na(taxa_roles)) %>% # filter based on modules (only those assigned to a module)
  filter(mean_rel_abund > 0.001) %>% # filter based on mean relative abundance
  # filter(percent_occupancy > 50) %>% # filter based on percent occupancy
  mutate(taxon = fct_reorder(taxon, mean_rel_abund)) # sort

# titan prep
labels <- c("alkalinity" = "Alkalinity", # set tidy labels for facets
            "DO" = "DO",
            "nitrate" = "Nitrate-N",
            "temp" = "Temperature")

change_points <- data.frame(variable = c("alkalinity", "DO", "nitrate", "temp"), # extract change points
                            fsumz_inc = c(titan_alk$sumz.cp["fsumz+", "cp"],
                                          titan_DO$sumz.cp["fsumz+", "cp"],
                                          titan_nitrate$sumz.cp["fsumz+", "cp"],
                                          titan_temp$sumz.cp["fsumz+", "cp"]),
                            fsumz_dec = c(titan_alk$sumz.cp["fsumz-", "cp"],
                                          titan_DO$sumz.cp["fsumz-", "cp"],
                                          titan_nitrate$sumz.cp["fsumz-", "cp"],
                                          titan_temp$sumz.cp["fsumz-", "cp"]))

grad_dec <- scales::col_numeric(palette = c("#d1eeea", "#2a5674"), domain = range(df_filtered$zscore)) # set colour scale for decreasers
grad_inc <- scales::col_numeric(palette = c("#fde0c5", "#eb4a40"), domain = range(df_filtered$zscore)) # set colour scale for increasers

df_filtered <- df_filtered %>%
  mutate(fill_color = ifelse(maxgrp == 1, 
                             grad_dec(zscore), # assign colour scales
                             grad_inc(zscore)))

all_taxa <- unique(df_filtered$taxon)
df_filtered$taxon <- factor(df_filtered$taxon, levels = all_taxa) # make sure individual plots have the same taxon list
df_filtered <- df_filtered %>% 
  mutate(taxon = fct_reorder(taxon, mean_rel_abund))

# titan ridge plots
titan_plot <- ggplot(df_filtered, aes(x = chk_pts, y = taxon, fill = fill_color)) +
  facet_wrap(~ variable, nrow=1, scales="free_x", labeller = labeller(variable=labels)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 3) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black")) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = "Change point", fill = "Z score") +
  geom_vline(data=change_points, aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data=change_points, aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

titan_plot

titan_alk_plot <- ggplot(data = df_filtered %>% filter(variable == "alkalinity"),
                         aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 3) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"), 
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Alkalinity CP (mg "*L^{-1}*" of CaCO"[3]*")"), fill = "Z score", title = "Alkalinity") +
  geom_vline(data = change_points %>% filter(variable == "alkalinity"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "alkalinity"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50))

titan_alk_plot

titan_DO_plot <- ggplot(data = df_filtered %>% filter(variable == "DO"),
                         aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("DO CP (mg "*L^{-1}*")"), fill = "Z score", title = "DO") +
  geom_vline(data = change_points %>% filter(variable == "DO"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "DO"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(limits = c(8.5, 12.5), breaks = seq(9, 12, by = 1))

titan_DO_plot

titan_nitrate_plot <- ggplot(data = df_filtered %>% filter(variable == "nitrate"),
                        aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Nitrate-N CP (mg "*L^{-1}*")"), fill = "Z score", title = "Nitrate-N") +
  geom_vline(data = change_points %>% filter(variable == "nitrate"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "nitrate"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2))

titan_nitrate_plot

titan_temp_plot <- ggplot(data = df_filtered %>% filter(variable == "temp"),
                        aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Temperature CP (°C)"), fill = "Z score", title = "Temperature") +
  geom_vline(data = change_points %>% filter(variable == "temp"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "temp"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(limits = c(6, 16.5), breaks = seq(6, 16, by = 2))

titan_temp_plot

legend_inc <- ggplot(data = subset(df_filtered, maxgrp == "2"), aes(x = chk_pts, y = taxon, fill = zscore)) +
  facet_wrap(~ variable) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme(legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40") +
  labs(fill = "Z score (tolerant)")
legend_inc <- patchwork::wrap_elements(full = cowplot::get_legend(legend_inc)) # extract legend

legend_dec <- ggplot(data = subset(df_filtered, maxgrp == "1"), aes(x = chk_pts, y = taxon, fill = zscore)) +
  facet_wrap(~ variable) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme(legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674") +
  labs(fill = "Z score (sensitive)")
legend_dec <- patchwork::wrap_elements(full = cowplot::get_legend(legend_dec)) # extract legend

# abundance
abundance_plot <- ggplot(df_filtered, aes(x = "Abundance", y = taxon, fill = mean_rel_abund)) +
  geom_tile(color = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(fill = "Relative abundance", title = "Abundance") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

abundance_plot

legend_abundance <- patchwork::wrap_elements(full = cowplot::get_legend(abundance_plot)) # extract legend

# occupancy
occupancy_plot <- ggplot(df_filtered, aes(x = "Prevalence", y = taxon, fill = percent_occupancy)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  labs(fill = "Prevalence (%)", title = "Prevalence") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

occupancy_plot

legend_occupancy <- patchwork::wrap_elements(full = cowplot::get_legend(occupancy_plot)) # extract legend

# modules
module_plot <- ggplot(df_filtered, aes(x = "Module role", y = taxon, fill = taxa_roles)) +
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
  labs(fill = "Module role", title = "Module role") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

module_plot

legend_module <- patchwork::wrap_elements(full = cowplot::get_legend(module_plot)) # extract legend

# taxonomy
phylum_colour_map <- read.table("output/phylum_colour_map.csv", check.names=F, row.names=1, header=T, sep = ",")
phylum_colours <- setNames(phylum_colour_map$colour, phylum_colour_map$Phylum)

tax_plot <- ggplot(df_filtered, aes(x = "Phylum", y = taxon, fill = Phylum)) +
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
  scale_fill_manual(guide = guide_legend(ncol = 2),
                    values = phylum_colours) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

tax_plot

legend_tax <- patchwork::wrap_elements(full = cowplot::get_legend(tax_plot)) # extract legend

# core
core_plot <- ggplot(df_filtered, aes(x = "Core", y = taxon)) +
  geom_point(data = df_filtered %>% filter(core == "Core"), shape = 42, size = 4) +
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
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

core_plot


# combined plot
plots <- core_plot + tax_plot + module_plot + abundance_plot + occupancy_plot + titan_alk_plot + titan_DO_plot + titan_nitrate_plot + titan_temp_plot +
  plot_layout(ncol = 9, widths = c(1, 1, 1, 1, 1, 10, 10, 10, 10), guides = "collect", axes = "collect")

legends <- legend_tax + legend_module + legend_abundance + legend_occupancy + legend_inc + legend_dec +
  plot_layout(widths=c(0.2,1,1.5,1.5,1.5,1.5))

combined <- plots / legends +
  plot_layout(guides = "collect", heights = c(6, 0.8)) &
  theme(legend.position = "none",
        plot.margin = margin(t = -1, r = 1, b = -1, l = 1))

combined

ggsave("output/figures/titan.pdf", combined, width = 180, height = 190, units = "mm", dpi = 600)
ggsave("output/figures/titan.png", combined, width = 180, height = 190, units = "mm", dpi = 600)


# stats
sppmax_modules_abund_tax %>% # unique ASVs from TITAN
  distinct(taxon) %>%     
  count()

sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are also in a module
  filter(!is.na(module)) %>%
  distinct(taxon) %>%     
  count()

sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are module hubs
  filter(taxa_roles == "Module hubs") %>%
  distinct(taxon) %>%     
  count()

sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are connectors
  filter(taxa_roles == "Connectors") %>%
  distinct(taxon) %>%     
  count()

sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are peripherals
  filter(taxa_roles == "Peripheral nodes") %>%
  distinct(taxon) %>%     
  count()

sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are core ASVs
  filter(core == "Core") %>%
  distinct(taxon) %>%     
  count()

unique <- sppmax_modules_abund_tax %>%
  group_by(taxon) %>%
  summarise(n_var = n_distinct(variable),
            variable = first(variable)) %>%
  filter(n_var == 1) %>%
  select(taxon, variable)

unique_df <- sppmax_modules_abund_tax %>%
  inner_join(unique, by = c("taxon", "variable"))

alkalinity_df <- unique_df %>% filter(variable == "alkalinity")
DO_df <- unique_df %>% filter(variable == "DO")
temp_df <- unique_df %>% filter(variable == "temp")
nitrate_df <- unique_df %>% filter(variable == "nitrate")

core_titan_modules <- sppmax_modules_abund_tax %>%
  filter(core == "Core", !is.na(module))

core_titan_modules_ASVs <- sppmax_modules_abund_tax %>% # unique ASVs from TITAN that are core ASVs and in modules
  filter(!is.na(module)) %>%
  filter(core == "Core") %>%
  distinct(taxon) %>% 
  merge(., abund_occ, by="taxon") %>% 
  merge(., tax_tab, by.x="taxon", by.y="row.names") %>% 
  count(Phylum)
