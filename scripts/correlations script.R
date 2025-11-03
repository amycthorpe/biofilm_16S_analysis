# biofilm analysis - 16S
# correlation matrix


# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(pals)
library(corrplot)
library(ggsci)
library(Hmisc)


set.seed(123)


# import data
meta_merged <- read.csv("data/meta_merged.csv", header = T, row.names = 1, check.names = F)

colnames(meta_merged) <- gsub("/", "_", colnames(meta_merged)) # replace spaces in colnames

columns_to_exclude <- c("Sample_type", "Sample_date", "RSN_ID", "season", "long", "lat", "Year",
                        # "Ammonia_N_unionised", "TN_oxidised", "DO_sat", "TP", # all co-correlate with something else/calculated from same measurement
                        "Alkalinity_grans", "Alkalinity") # replaced by alkalinity_combined (results from both methods)

meta_merged2 <- meta_merged[ , !(names(meta_merged) %in% columns_to_exclude)]


# correlation matrix - metadata
correlation_matrix <- cor(meta_merged2, use = "pairwise.complete.obs", method="spearman")
corrplot(correlation_matrix, method = "circle", tl.col = "black", tl.srt = 45)

rcors <- rcorr(as.matrix(meta_merged2), type="spearman") # get r values
rvalues <- as.matrix(rcors[["r"]])
rvalues <- reshape2::melt(rvalues) %>% 
  filter(value != 1)

pvalues <- as.matrix(rcors[["P"]]) # get p values
pvalues <- reshape2::melt(pvalues) %>% 
  filter(across(everything(), ~ !is.na(.)))

cors <- merge(rvalues, pvalues, by=c("Var1", "Var2"))
colnames(cors)[3] <- "r_value"
colnames(cors)[4] <- "p_value"

cors$significance <- ifelse(cors$p_value <= 0.001, "***", # add significance levels
                            ifelse(cors$p_value <= 0.01, "**",
                                   ifelse(cors$p_value <= 0.05, "*",
                                          "")))


# heatmap - metadata
all_vars <- sort(unique(c(cors$Var1, cors$Var2)))

cors <- cors %>%
  mutate(Var1 = factor(Var1, levels = all_vars),
         Var2 = factor(Var2, levels = all_vars))

cors_upper <- cors %>%
  filter(as.numeric(Var1) <= as.numeric(Var2))

colours <- brewer.pal(11, "RdBu")

label_data <- cors_upper %>%
  group_by(Var2) %>%
  filter(as.numeric(Var1) == max(as.numeric(Var1))) # labels at last tile in each row

heatmap <- ggplot(cors_upper, aes(x = Var1, y = Var2, fill = r_value)) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "lines")) +
  geom_text(aes(label = significance), size = 2) +
  scale_fill_gradientn(colours = rev(colours), limits = c(-1, 1), values = scales::rescale(c(-1, 0, 1))) +
  geom_text(data = label_data, aes(x = Var1, y = Var2, label = Var1),
            angle = -90, hjust = 0, vjust = 0.5, size = 5/.pt, nudge_x = 0, nudge_y = -1) +
  coord_cartesian(clip = "off", xlim = c(0.5, max(as.numeric(cors_upper$Var1)) + 1))

heatmap


# correlation matrix - water chemistry vs diversity, richness, phyla, modules
diversity <- read.csv("output/diversity.csv", check.names=F, header=T, row.names=1, sep=",") %>% 
  select('16S_ID', Shannon, Observed) %>% 
  rename(Sample = '16S_ID')

phylum_abundance <- read.csv("output/phylum_abundance.csv", check.names=F, header=T, row.names=1, sep=",") %>% 
  filter(., Phylum != "x__Others") %>% 
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = 0)
  
modules <- read.csv("output/modules.csv", check.names=F, header=T, row.names=1, sep=",") %>% 
  select(name, module)

ASV_tab <- data.frame(fread("data/ASV_counts_16S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1) %>% 
  t(.) %>% # transpose so ASVs are rows
  as.data.frame(.)

module_abundance <- ASV_tab %>% 
  rownames_to_column(., var="name") %>% 
  merge(., modules, by = "name") %>% 
  group_by(module) %>% 
  summarise(across(where(is.numeric), sum)) %>% # sum abundance by module
  column_to_rownames(., var = "module") %>%
  mutate(across(everything(), ~ . / sum(.))) %>% # convert to relative abundance
  t(.) # samples as rows

data <- diversity %>% 
  merge(., phylum_abundance, by = "Sample") %>% 
  column_to_rownames(., var = "Sample") %>% 
  merge(., module_abundance, by = "row.names") %>%
  column_to_rownames(., var = "Row.names")

meta_merged3 <- meta_merged2 %>% 
  merge(., data, by="row.names") %>% 
  column_to_rownames(., var = "Row.names")

correlation_matrix <- cor(meta_merged3, use = "pairwise.complete.obs", method="spearman")
corrplot(correlation_matrix, method = "circle", tl.col = "black", tl.srt = 45)

rcors <- rcorr(as.matrix(meta_merged3), type="spearman") # get r values
rvalues <- as.matrix(rcors[["r"]])
rvalues <- reshape2::melt(rvalues) %>% 
  filter(value != 1)

pvalues <- as.matrix(rcors[["P"]]) # get p values
pvalues <- reshape2::melt(pvalues) %>% 
  filter(across(everything(), ~ !is.na(.)))

cors <- merge(rvalues, pvalues, by=c("Var1", "Var2"))
colnames(cors)[3] <- "r_value"
colnames(cors)[4] <- "p_value"

cors$significance <- ifelse(cors$p_value <= 0.001, "***", # add significance levels
                            ifelse(cors$p_value <= 0.01, "**",
                                   ifelse(cors$p_value <= 0.05, "*",
                                          "")))

write.csv(cors, "output/correlations.csv")


# heatmap -  water chemistry vs diversity, richness, phyla, modules
cors <- cors %>% filter(Var1 %in% colnames(meta_merged2))
cors <- cors %>% filter(!Var2 %in% colnames(meta_merged2))

heatmap <- ggplot(cors, aes(x = Var1, y = Var2, fill = r_value)) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle=45, hjust=1),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "lines")) +
  geom_text(aes(label = significance), size = 2) +
  scale_fill_gradientn(colours = rev(colours), limits = c(-1, 1), values = scales::rescale(c(-1, 0, 1)))

heatmap


# co-correlation heatmap
cors <- merge(rvalues, pvalues, by=c("Var1", "Var2"))
colnames(cors)[3] <- "r_value"
colnames(cors)[4] <- "p_value"

cors$significance <- ifelse(cors$p_value <= 0.001, "***", # add significance levels
                            ifelse(cors$p_value <= 0.01, "**",
                                   ifelse(cors$p_value <= 0.05, "*",
                                          "")))

all_vars <- sort(unique(c(cors$Var1, cors$Var2)))

modules_order <- rev(paste0("M", 1:20))
phyla_order <- rev(sort(colnames(phylum_abundance)[-1]))
water_chem_order <- rev(sort(colnames(meta_merged2)))
diversity_order <- c("Shannon", "Observed")
all_vars_order <- c(water_chem_order, modules_order, diversity_order, phyla_order)

cors <- cors %>%
  mutate(Var1 = factor(Var1, levels = all_vars_order),
         Var2 = factor(Var2, levels = all_vars_order))

cors_upper <- cors %>%
  filter(as.numeric(Var1) <= as.numeric(Var2)) %>% 
  mutate(Var1 = recode(Var1,
                "Observed" = "ASV richness",
                "Shannon" = "Shannon diversity",
                "Alkalinity_combined" = "Alkalinity",
                "Chloride" = "Chloride ion",
                "Water_temp" = "Temp",
                "TN_oxidised" = "TON",
                "Nitrite_N" = "Nitrite-N",
                "Nitrate_N" = "Nitrate-N",
                "DO_sat" = "DO (saturated)",
                "Ammoniacal_N" = "Ammoniacal-N",
                "Ammonia_N_unionised" = "Ammonia-N")) %>%
  mutate(Var2 = recode(Var2,
                "Observed" = "ASV richness",
                "Shannon" = "Shannon diversity",
                "Alkalinity_combined" = "Alkalinity",
                "Chloride" = "Chloride ion",
                "Water_temp" = "Temp",
                "TN_oxidised" = "TON",
                "Nitrite_N" = "Nitrite-N",
                "Nitrate_N" = "Nitrate-N",
                "DO_sat" = "DO (saturated)",
                "Ammoniacal_N" = "Ammoniacal-N",
                "Ammonia_N_unionised" = "Ammonia-N"))

colours <- brewer.pal(11, "RdBu")

label_data <- cors_upper %>%
  group_by(Var2) %>%
  filter(as.numeric(Var1) == max(as.numeric(Var1))) # labels at last tile in each row

heatmap <- ggplot(cors_upper, aes(x = Var1, y = Var2, fill = r_value)) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black", margin = margin(r = -4)),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.size = unit(0.8, "lines"), 
        legend.position = c(0.3, 0.02),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  geom_text(aes(label = significance), size = 1.5) +
  scale_fill_gradientn(colours = rev(colours), limits = c(-1, 1), values = scales::rescale(c(-1, 0, 1))) +
  geom_text(data = label_data, aes(x = Var1, y = Var2, label = Var1),
            angle = -90, hjust = 0, vjust = 0.5, size = 5/.pt, nudge_x = 0, nudge_y = -1) +
  coord_cartesian(clip = "off", xlim = c(0.5, max(as.numeric(cors_upper$Var1)) + 1)) +
  labs(fill = "Spearman correlation")

heatmap

ggsave("output/figures/heatmap_all.png", heatmap, width=180, height=160, units="mm", dpi=600)
ggsave("output/figures/heatmap_all.pdf", heatmap, width=180, height=180, units="mm", dpi=600)



