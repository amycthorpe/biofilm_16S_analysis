# biofilm analysis - 16S
# NMDS


# load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(geosphere)
library(tidyverse)
library(reshape2)
library(data.table)
library(patchwork)
library(parallel)
library(ggrepel)
library(RColorBrewer)
library(phyloseq)


set.seed(123)


# import data 
ASV_tab <- data.frame(fread("data/ASV_counts_16S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table("data/ASV_tax_16S_processed.csv", header=T, row.names=1, check.names=F, sep=",")
meta_all <- read.table("data/seq_metadata_all_amplicons.csv", header=T, check.names=F, sep=",")

meta <- meta_all[meta_all$`16S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon) - number of samples is less due to those lost during rarefaction
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="16S_ID") # make seq ID rownames

meta <- meta %>% # add season column
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta <- select(meta, -c(1:6)) # remove unwanted ID cols


# convert to relative abundance
ASV_rel <- ASV_tab %>%
  t(.) %>% 
  as.data.frame(.) %>% 
  mutate(across(everything(), ~ . / sum(.))) %>% 
  t(.) # samples as rows


# import environmental metadata
water_chem <- read.csv("data/wims_summaries_3months 2.csv", check.names=F, header=T) %>% 
  select(2, matches("(_mean|_combined)$")) # only mean columns

new_names <- c(
  "Extended_ID",
  "Alkalinity_combined",
  "Alkalinity_grans",
  "Alkalinity",
  "Ammonia_N_unionised",
  "Ammoniacal_N",
  "DOC",
  "Chloride",
  "Conductivity",
  "Nitrate_N",
  "Nitrite_N",
  "TN_oxidised",
  "TN",
  "Orthophosphate",
  "DO",
  "DO_sat",
  "TP",
  "SiO2",
  "Water_temp",
  "pH")

colnames(water_chem) <- new_names # tidy colnames

water_chem <- water_chem %>% select(-TN) # TN has lots of NAs so remove

water_chem <- water_chem %>% 
  filter(Extended_ID %in% rownames(ASV_rel)) %>% # subset to same samples
  column_to_rownames(., var="Extended_ID")

meta_merged <- merge(meta, water_chem, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

write.csv(meta_merged, "data/meta_merged.csv")

colnames(meta_merged) <- gsub("/", "_", colnames(meta_merged))


# NMDS
ASV_hel <- decostand(ASV_tab, method = "hellinger")
nmds_hel <- metaMDS(ASV_hel, autotransform = F)
nmds_hel

plot(nmds_hel, display="sites")

nmds <- nmds_hel

columns_to_exclude <- c("Sample_type", "Sample_date", "RSN_ID", "long", "lat", "Year",
                        "Ammonia_N_unionised", "TN_oxidised", "DO_sat", "TP", # all co-correlate with something else/calculated from same measurement
                        "Alkalinity_grans", "Alkalinity") # replaced by alkalinity_combined (results from both methods)
columns_to_include <- setdiff(colnames(meta_merged), columns_to_exclude)
formula <- as.formula(paste("nmds ~", paste(columns_to_include, collapse="+"))) # set envfit formula

meta_merged <- meta_merged[match(rownames(ASV_rel), rownames(meta_merged)),]

nmds.fit <- envfit(formula,
                   data=meta_merged,
                   perm=999,
                   na.rm=T)

nmds.fit


# NMDS with ggplot
df <- nmds.fit$vectors # get datafrane of vector r and p values
df_r <- as.data.frame(df[["r"]])
df_pvals <- as.data.frame(df[["pvals"]])
df_arrows <- as.data.frame(df[["arrows"]])
nmds.fit_df <- cbind(df_arrows,df_r,df_pvals)
colnames(nmds.fit_df)[3] <- "rvalue"
colnames(nmds.fit_df)[4] <- "pvalue"

nmds.fit_df <- nmds.fit_df %>% # add significance levels
  mutate(significance = case_when(
    pvalue > 0.05 ~ "ns",
    pvalue <= 0.05 & pvalue > 0.01 ~ "*",
    pvalue <= 0.01 & pvalue > 0.001 ~ "**",
    pvalue <= 0.001 ~ "***"
  ))

write.csv(nmds.fit_df, "output/envfit.csv")

# nmds.fit_sub <- nmds.fit_df %>% # filter vectors
#   filter(., pvalue <= 0.05, rvalue >= 0.1)

nmds_scores <- as.data.frame(nmds[["points"]]) # extract scores
nmds_scores_meta <- merge(nmds_scores,meta_merged,by="row.names") %>% # combine with meta
  column_to_rownames(.,"Row.names")

vectors <- as.data.frame(scores(nmds.fit, "vectors")) * ordiArrowMul(nmds.fit) # set up vector arrows
arrowhead = arrow(length = unit(0.02, "npc"))


# plot with vectors
nmds_scores_meta$season <- factor(nmds_scores_meta$season, levels = c("Spring", "Autumn"))

new_names <- c("Alkalinity_combined" = "Alkalinity",
               "Ammoniacal_N" = "Ammoniacal-N",
               "Nitrate_N" = "Nitrate-N",
               "Nitrite_N" = "Nitrite-N",
               # "SiO2" = expression(SiO[2]),
               "Water_temp" = "Temperature")

vectors2 <- vectors
matched <- intersect(rownames(vectors2), names(new_names))
rownames(vectors2)[rownames(vectors2) %in% matched] <- new_names[matched]

plot_all <- ggplot(nmds_scores_meta, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = season)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  geom_vline(xintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_hline(yintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_segment(data = vectors2, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size = 0.3, arrow = arrowhead, colour = "black") +
  # geom_text(data = vectors2, aes(x = 1.1*NMDS1, y = 1.1*NMDS2),
  #           label = row.names(vectors2), size = 2, colour="black") +
  # geom_label(data = vectors2, aes(x = 1.2*NMDS1, y = 1.2*NMDS2),
  #            label = row.names(vectors2), size = 2, colour="black") +
  geom_text_repel(data = vectors2, aes(x = NMDS1, y = NMDS2),
                  label = row.names(vectors2), size = 2, colour = "black",
                  max.overlaps = Inf, segment.colour = NA) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_colour_manual(name = "Season", values = c("forestgreen", "steelblue"))

plot_all


# plot coloured by variables
variables <- rownames_to_column(vectors, var="variable")
variable_list <- unique(variables$variable)

plot_list <- list()

for (var in variable_list) {
  plot <- ggplot(nmds_scores_meta, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = .data[[var]])) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 5, colour = "black"),
          legend.title = element_text(size = 6, colour = "black"),
          axis.text = element_text(size = 5, colour = "black"),
          axis.ticks = element_line(size = 0.2),
          axis.title = element_text(size = 6, colour = "black"),
          plot.title = element_text(size = 6, colour = "black")) +
    geom_vline(xintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
    geom_hline(yintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
    labs(x = "NMDS1", y = "NMDS2", title = var) +
    scale_colour_viridis_c(option="D", na.value="grey40")
  
  plot_list[[var]] <- plot
  
}

plots_sub <- plot_list[names(plot_list) %in% variable_list]
plots_sub <- lapply(plots_sub, \(p) p + theme(legend.position = "none"))
combined <- wrap_plots(plots_sub, ncol = 4)


# plot coloured by variables on a log scale
nmds_scores_meta_log <- nmds_scores_meta %>%
  mutate(across(10:27, ~ log(. + 1)))

variables <- rownames_to_column(vectors, var="variable")
variable_list <- unique(variables$variable)

plot_list_log <- list()

for (var in variable_list) {
  plot <- ggplot(nmds_scores_meta_log, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = .data[[var]])) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 5, colour = "black"),
          legend.title = element_text(size = 6, colour = "black"),
          axis.text = element_text(size = 5, colour = "black"),
          axis.ticks = element_line(size = 0.2),
          axis.title = element_text(size = 6, colour = "black"),
          plot.title = element_text(size = 6, colour = "black")) +
    geom_vline(xintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
    geom_hline(yintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
    labs(x = "NMDS1", y = "NMDS2", title = var) +
    scale_colour_viridis_c(option="D", na.value="grey40")
  
  plot_list_log[[var]] <- plot

}

plots_sub <- plot_list_log[names(plot_list_log) %in% variable_list]
plots_sub <- lapply(plots_sub, \(p) p + theme(legend.position = "none"))
combined <- wrap_plots(plots_sub, ncol = 4)


# arrange plot
nmds_all <- plot_all +
  theme(legend.direction = "horizontal",
        legend.title.position = "left",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_colour_manual(name = "Season", values = c("forestgreen", "steelblue")) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = bquote("Stress =" ~ .(round(nmds$stress, 2))),
           size = 2.2)
legend_season <- patchwork::wrap_elements(full = cowplot::get_legend(nmds_all)) # extract legend

nmds_alk <- plot_list_log[["Alkalinity_combined"]] +
  theme(legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_colour_viridis_c(name = "Alkalinity (log)", option="viridis", na.value="grey40") +
  # scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "Purples"), name = "Alkalinity (log)", na.value = "grey40") +
  labs(title = "Alkalinity") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = bquote(R^2 == .(round(nmds.fit_df["Alkalinity_combined", "rvalue"], 2)) * "," ~ .(nmds.fit_df["Alkalinity_combined", "significance"])),
           size = 2.2)
legend_alk <- patchwork::wrap_elements(full = cowplot::get_legend(nmds_alk)) # extract legend

nmds_DO <- plot_list[["DO"]] +
  theme(legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_colour_viridis_c(name = "DO", option="viridis", na.value="grey40") +
  # scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "Blues"), name = "DO", na.value = "grey40") +
  labs(title = "DO") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = bquote(R^2 == .(round(nmds.fit_df["DO", "rvalue"], 2)) * "," ~ .(nmds.fit_df["DO", "significance"])),
           size = 2.2)
legend_DO <- patchwork::wrap_elements(full = cowplot::get_legend(nmds_DO)) # extract legend

nmds_nitrate <- plot_list_log[["Nitrate_N"]] +
  theme(legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_colour_viridis_c(name = "Nitrate-N (log)", option="viridis", na.value="grey40") +
  # scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "Greens"), name = "Nitrate-N (log)", na.value = "grey40") +
  labs(title = "Nitrate-N") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = bquote(R^2 == .(round(nmds.fit_df["Nitrate_N", "rvalue"], 2)) * "," ~ .(nmds.fit_df["Nitrate_N", "significance"])),
           size = 2.2)
legend_nitrate <- patchwork::wrap_elements(full = cowplot::get_legend(nmds_nitrate)) # extract legend

nmds_temp <- plot_list[["Water_temp"]] +
  theme(legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_colour_viridis_c(name = "Temperature", option="viridis", na.value="grey40") +
  # scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "Oranges"), name = "Temperature", na.value = "grey40") +
  labs(title = "Temperature") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = bquote(R^2 == .(round(nmds.fit_df["Water_temp", "rvalue"], 2)) * "," ~ .(nmds.fit_df["Water_temp", "significance"])),
           size = 2.2)
legend_temp <- patchwork::wrap_elements(full = cowplot::get_legend(nmds_temp)) # extract legend

design <- c(area(t = 1, l = 1, b = 2, r = 2),  # big plot (2 rows, 2 cols)
            area(t = 1, l = 3, b = 1, r = 3),  # small_plot1
            area(t = 2, l = 3, b = 2, r = 3),  # small_plot2
            area(t = 3, l = 1, b = 3, r = 1),  # small_plot3 (under big plot col 1)
            area(t = 3, l = 2, b = 3, r = 2),  # small_plot4 (under big plot col 2)
            area(t = 3, l = 3, b = 3, r = 3))  # small_plot5 (under small 1/2 col 3)

legends1 <- legend_alk + legend_DO + legend_nitrate + legend_temp +
  plot_layout(ncol = 2)

legends <- legend_season + legends1 +
  plot_layout(ncol = 1, heights=c(1,4))

patchwork <- wrap_plots(
  nmds_all, # big plot
  nmds_alk, # small 1
  nmds_DO, # small 2
  nmds_nitrate, # small 3
  nmds_temp, # small 4
  legends, # small 5
  design = design) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

patchwork

ggsave("output/figures/NMDS_patchwork.png", patchwork, width = 180, height = 180, units = "mm", dpi = 600)
ggsave("output/figures/NMDS_patchwork.pdf", patchwork, width = 180, height = 180, units = "mm", dpi = 600)


# adonis
columns_to_exclude <- c("Sample_type", "Sample_date", "RSN_ID",
                        "Ammonia_N_unionised", "TN_oxidised", "DO_sat", "TP", # all co-correlate with something else/calculated from same measurement
                        "Alkalinity_grans", "Alkalinity") # replaced by alkalinity_combined (results from both methods)
columns_to_include <- setdiff(colnames(meta_merged), columns_to_exclude)

meta_sub <- meta_merged[ , (names(meta_merged) %in% columns_to_include)]
meta_sub <- na.omit(meta_sub)
ASV_rel2 <- as.data.frame(ASV_rel[rownames(ASV_rel) %in% rownames(meta_sub), ])
ASV_rel2 <- ASV_rel2[match(rownames(meta_sub), rownames(ASV_rel2)), ] # match order of rows

sites <- meta_merged %>% 
  select(RSN_ID)
sites <- sites %>% 
  merge(., meta_sub, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

meta_sub$Year <- as.factor(meta_sub$Year)
meta_sub$season <- as.factor(meta_sub$season)
sites$RSN_ID <- as.factor(sites$RSN_ID)

dist <- vegdist(ASV_rel2)

# formula <- as.formula(paste("dist ~", paste(columns_to_include, collapse="+"))) # set adonis formula
formula <- dist ~ Year + season
formula <- dist ~ Year * season
# formula <- dist ~ Year / season

adonis_result <- adonis2(formula,
                         data = meta_sub,
                         method = "bray",
                         strata = sites$RSN_ID,
                         permutations = 999,
                         by="margin")

adonis_result


# taxa biplot
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
ASV <- otu_table(ASV_hel, taxa_are_rows = F)
samples <- sample_data(meta)
ps <- phyloseq(ASV, TAX, samples)
ps

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# phylum level
ps_phylum<-tax_glom(ps_rel, "Phylum")

otu_tab <- as.data.frame(ps_phylum@otu_table)
taxonomy <- as.data.frame(ps_phylum@tax_table)

nmds.fit_tax <- envfit(nmds, otu_tab,
                       perm=999,
                       na.rm=T)

nmds.fit_tax

# NMDS with ggplot
df <- nmds.fit_tax$vectors # get datafrane of vector r and p values
df_r <- as.data.frame(df[["r"]])
df_pvals <- as.data.frame(df[["pvals"]])
df_arrows <- as.data.frame(df[["arrows"]])
nmds.fit_df <- cbind(df_arrows,df_r,df_pvals)
colnames(nmds.fit_df)[3] <- "rvalue"
colnames(nmds.fit_df)[4] <- "pvalue"

nmds.fit_df <- nmds.fit_df %>% # add significance levels
  mutate(significance = case_when(
    pvalue > 0.05 ~ "ns",
    pvalue <= 0.05 & pvalue > 0.01 ~ "*",
    pvalue <= 0.01 & pvalue > 0.001 ~ "**",
    pvalue <= 0.001 ~ "***"
  ))

write.csv(nmds.fit_df, "output/envfit_phylum.csv")

nmds.fit_sub <- nmds.fit_df %>% # filter vectors
  filter(., pvalue <= 0.05, rvalue >= 0.2)

nmds_scores <- as.data.frame(nmds[["points"]]) # extract scores
nmds_scores_meta <- merge(nmds_scores,meta_merged,by="row.names") %>% # combine with meta
  column_to_rownames(.,"Row.names")

vectors <- as.data.frame(scores(nmds.fit_tax, "vectors")) * ordiArrowMul(nmds.fit_tax) # set up vector arrows
arrowhead = arrow(length = unit(0.02, "npc"))
vectors_sub <- vectors[rownames(vectors) %in% rownames(nmds.fit_sub), ]
vectors_sub <- merge(vectors_sub, taxonomy, by="row.names")

nmds_scores_meta$season <- factor(nmds_scores_meta$season, levels = c("Spring", "Autumn"))

phylum <- ggplot(nmds_scores_meta, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = season)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = "none",
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  geom_vline(xintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_hline(yintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_segment(data = vectors_sub, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size = 0.3, arrow = arrowhead, colour = "black") +
  # geom_text(data = vectors, aes(x = 1.1*NMDS1, y = 1.1*NMDS2),
  #           label = row.names(vectors), size = 2, colour="black") +
  # geom_label(data = vectors, aes(x = 1.2*NMDS1, y = 1.2*NMDS2),
  #            label = row.names(vectors), size = 2, colour="black") +
  geom_text_repel(data = vectors_sub, aes(x = NMDS1, y = NMDS2),
                  label = vectors_sub$Phylum, size = 1.5, colour = "black",
                  max.overlaps = Inf, segment.colour = NA) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_colour_manual(name = "Season", values = c("forestgreen", "steelblue")) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1))

phylum


# genus level
ps_genera<-tax_glom(ps_rel, "Genus")

otu_tab <- as.data.frame(ps_genera@otu_table)
taxonomy <- as.data.frame(ps_genera@tax_table)

taxonomy_sub <- taxonomy %>%
  filter(Genus != "unassigned") 

otu_tab <- t(otu_tab)
otu_tab_sub <- otu_tab[rownames(otu_tab) %in% rownames(taxonomy_sub), ]
otu_tab_sub <- otu_tab_sub %>% 
  t() %>% 
  as.data.frame()

top_genera <- otu_tab_sub %>% 
  summarise(across(everything(), mean)) %>% 
  pivot_longer(everything(), names_to = "ASV", values_to = "mean_rel_abund") %>% 
  arrange(desc(mean_rel_abund)) %>% 
  slice_head(n=100) %>% 
  pull(ASV)

otu_tab_sub <- otu_tab_sub %>% select(all_of(top_genera))

nmds.fit_tax <- envfit(nmds, otu_tab_sub,
                       perm=999,
                       na.rm=T)

nmds.fit_tax

# NMDS with ggplot
df <- nmds.fit_tax$vectors # get datafrane of vector r and p values
df_r <- as.data.frame(df[["r"]])
df_pvals <- as.data.frame(df[["pvals"]])
df_arrows <- as.data.frame(df[["arrows"]])
nmds.fit_df <- cbind(df_arrows,df_r,df_pvals)
colnames(nmds.fit_df)[3] <- "rvalue"
colnames(nmds.fit_df)[4] <- "pvalue"

nmds.fit_df <- nmds.fit_df %>% # add significance levels
  mutate(significance = case_when(
    pvalue > 0.05 ~ "ns",
    pvalue <= 0.05 & pvalue > 0.01 ~ "*",
    pvalue <= 0.01 & pvalue > 0.001 ~ "**",
    pvalue <= 0.001 ~ "***"
  ))

write.csv(nmds.fit_df, "output/envfit_genera.csv")

nmds.fit_sub <- nmds.fit_df %>% # filter vectors
  filter(., pvalue <= 0.05, rvalue >= 0.2)

nmds_scores <- as.data.frame(nmds[["points"]]) # extract scores
nmds_scores_meta <- merge(nmds_scores,meta_merged,by="row.names") %>% # combine with meta
  column_to_rownames(.,"Row.names")

vectors <- as.data.frame(scores(nmds.fit_tax, "vectors")) * ordiArrowMul(nmds.fit_tax) # set up vector arrows
arrowhead = arrow(length = unit(0.02, "npc"))
vectors_sub <- vectors[rownames(vectors) %in% rownames(nmds.fit_sub), ]
vectors_sub <- merge(vectors_sub, taxonomy, by="row.names")

nmds_scores_meta$season <- factor(nmds_scores_meta$season, levels = c("Spring", "Autumn"))

genus <- ggplot(nmds_scores_meta, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = season)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = "none",
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  geom_vline(xintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_hline(yintercept = c(0,0), linetype = "dashed", colour = "black", size = 0.4) +
  geom_segment(data = vectors_sub, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size = 0.3, arrow = arrowhead, colour = "black") +
  # geom_text(data = vectors, aes(x = 1.1*NMDS1, y = 1.1*NMDS2),
  #           label = row.names(vectors), size = 2, colour="black") +
  # geom_label(data = vectors, aes(x = 1.2*NMDS1, y = 1.2*NMDS2),
  #            label = row.names(vectors), size = 2, colour="black") +
  geom_text_repel(data = vectors_sub, aes(x = NMDS1, y = NMDS2),
                  label = vectors_sub$Genus, size = 1.5, colour = "black",
                  max.overlaps = Inf, segment.colour = NA) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_colour_manual(name = "Season", values = c("forestgreen", "steelblue")) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1))

genus

patchwork <- wrap_plots(phylum, genus) +
  theme(legend.position = "none") +
  plot_layout(axes = "collect")

patchwork

ggsave("output/figures/NMDS_tax_patchwork.png", patchwork, width = 180, height = 90, units = "mm", dpi = 600)
ggsave("output/figures/NMDS_tax_patchwork.pdf", patchwork, width = 180, height = 90, units = "mm", dpi = 600)
