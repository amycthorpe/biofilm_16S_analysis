# biofilm analysis - 16S
# network modules


# load libraries
library(ggplot2)
library(dplyr)
library(microeco)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(data.table)
library(patchwork)
library(file2meco)
library(rgexf)
library(RColorBrewer)
library(ggsci)
library(igraph)
library(sf)
library(glue)

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

ASV_sub <- ASV_tab[, colSums(ASV_tab > 0 ) >= 156] # minimum occurrence frequency to present in >10% sites


# import environmental metadata
meta_merged <- read.csv("data/meta_merged.csv", header = T, row.names = 1, check.names = F)


# phyloseq to microeco object
TAX <- tax_table(as.matrix(tax_tab))
ASV <- otu_table(ASV_sub, taxa_are_rows = F)
samples <- sample_data(meta_merged)
ps <- phyloseq(ASV, TAX, samples)
ps

# convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# group at phylum level
ps_phylum<-tax_glom(ps_rel, "Phylum")

# melt
melt_rel<-psmelt(ps_phylum)

# microeco object
meco <- phyloseq2meco(ps_rel)


# compute network modules
t1 <- trans_network$new(dataset = meco,
                        cor_method = "spearman",
                        filter_thres = 0.0001, # 1249 ASVs after filtering at 10% prevalence and 0.0001 relative abundance
                        taxa_level = "OTU",
                        nThreads = 10)

t1$cal_network(p_thres = 0.01, COR_cut = 0.4)
t1$cal_module(method = "cluster_fast_greedy")
t1$get_node_table(node_roles = TRUE)
modules <- t1$res_node_table

modules <- modules %>%
  group_by(module) %>%
  filter(n() > 2) %>%  # keep only modules with more than 2 taxa
  ungroup()

# export network modules
write.csv(modules, "output/modules.csv")
t1$save_network(filepath = "output/module_network.gexf")

# network properties
t1$get_edge_table()
t1$get_adjacency_matrix()

t1$cal_network_attr()
t1$res_network_attr

nodes_to_keep <- modules$name
filtered_net <- t1[["subset_network"]](node = nodes_to_keep, return_igraph = TRUE)
t1$res_network <- filtered_net
t1$cal_network_attr()

# module composition
# relative abundance
modules_rel <- modules %>% 
  group_by(module) %>%
  mutate(rel_abundance = Abundance / sum(Abundance))

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

modules_rel <- modules_rel %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>% 
  merge(., tax_tab2, by.x="name", by.y="taxon")

  
# top phyla
top_phyla <- melt_rel %>% 
  group_by(Phylum) %>% 
  summarise(mean = mean(Abundance)) %>% 
  arrange(desc(mean)) %>%
  filter(mean >= 0.005) # %>% # filter by abundance threshold
# slice_head(n=10) # top 10

modules_top <- modules_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others")) %>% 
  filter(Abundance > 0) %>% # remove zeros
  group_by(Phylum, module) %>% 
  summarise(across(where(is.numeric), sum)) # sum by phylum

# set colours
phylum_colour_map <- read.table("output/phylum_colour_map.csv", check.names=F, row.names=1, header=T, sep = ",")
phylum_colours <- setNames(phylum_colour_map$colour, phylum_colour_map$Phylum)

modules_top$module <- factor(modules_top$module, levels = paste0("M", 1:20)) # order modules
  
module_plot <- ggplot(modules_top, aes(x = module, y = rel_abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  labs(y = "Relative abundance", x = "Module") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = phylum_colours,
                    guide = guide_legend(ncol = 6),
                    labels = function(x) ifelse(x == "x__Others", "Others", x))

module_plot


# node role relative
role_rel <- modules %>% 
  group_by(taxa_roles) %>%
  mutate(rel_abundance = Abundance / sum(Abundance))

role_rel <- role_rel %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>% 
  merge(., tax_tab2, by.x="name", by.y="taxon")

role_top <- role_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others")) %>% 
  filter(Abundance > 0) %>% # remove zeros
  group_by(Phylum, taxa_roles) %>% 
  summarise(across(where(is.numeric), sum)) # sum by phylum


role_plot <- ggplot(role_top, aes(x = taxa_roles, y = rel_abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  labs(y = "Relative abundance", x = "Module role") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = phylum_colours,
                    guide = guide_legend(ncol = 6),
                    labels = function(x) ifelse(x == "x__Others", "Others", x))

role_plot

# combined
module_plots <- module_plot + role_plot +
  plot_layout(ncol=2, axes = "collect", guides = "collect", widths = c(20,3)) &
  theme(legend.position = "bottom")

module_plots

ggsave("output/figures/module_plots.png", module_plots, width = 180, height = 90, units = "mm", dpi = 600)
ggsave("output/figures/module_plots.pdf", module_plots, width = 180, height = 90, units = "mm", dpi = 600)


# stats
stats <- modules_rel %>%
  group_by(module) %>%
  summarise(ASVs = n_distinct(name))

phylum_abundance <- modules_rel %>% 
  group_by(module, Phylum) %>% 
  summarise(phylum_abund = sum(rel_abundance))

role_stats <- modules_rel %>%
  group_by(module, taxa_roles) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = taxa_roles,
                     values_from = count,
                     values_fill = 0)


# network plot
# import data (first open .gexf file in text editor and find and replace " a="1"" with nothing throughout)
gx <- rgexf::read.gexf("output/module_network_edit.gexf")


# extract nodes
nodes <- gx$nodes
id_col <- "id"
nodes <- nodes[nodes$label %in% modules$name, ] # remove ASVs only in single pair modules


# extract edges
edges <- gx$edges
el_from <- "source"
el_to   <- "target"


# normalise columns
edges$from <- trimws(as.character(edges[[el_from]])) # start node of each edge
edges$to   <- trimws(as.character(edges[[el_to]])) # end node of each edge
nodes$..key <- trimws(as.character(nodes[[id_col]])) # ID for each node

edges <- edges[edges$from %in% nodes$..key & edges$to %in% nodes$..key,] # remove edges that only connect to an ASV in a single pair module


# check for missing ends in edges not found in nodes
missing_ends <- setdiff(unique(c(edges$from, edges$to)), nodes$..key)
if (length(missing_ends)) {
  stop("GEXF endpoints not found in node IDs. Examples: ",
       paste(utils::head(missing_ends, 10), collapse=", "))
}


# build the network
el <- as.matrix(edges[, c("from","to")]) # make edge list
g  <- graph_from_edgelist(el, directed = FALSE) # convert to igraph object
nodes$name <- nodes$..key # link to node metadata
idx <- match(V(g)$name, nodes$name)


# set vertex attributes
set_vattrs <- function(g, nodes, idx, drop_cols = c("..key","name")) {
  v_keep <- setdiff(names(nodes), drop_cols)
  for (att in v_keep) {
    vals <- nodes[[att]]
    if (is.factor(vals)) vals <- as.character(vals)
    g <- set_vertex_attr(g, name = att, value = vals[idx])
  }
  g
}

g <- set_vattrs(g, nodes, idx)


# set edge attributes
e_keep <- setdiff(names(edges), c("from","to"))
for (att in e_keep) {
  vals <- edges[[att]]
  if (is.factor(vals)) vals <- as.character(vals)
  if (att == "weight") vals <- suppressWarnings(as.numeric(vals))
  g <- set_edge_attr(g, name = att, value = vals)
}


# check
cat("Vertices:", vcount(g), "  Edges:", ecount(g), "\n")
cat("Vertex attributes:", paste(vertex_attr_names(g), collapse=", "), "\n")
cat("Edge attributes:", paste(edge_attr_names(g), collapse=", "), "\n")


# add modules
modules <- read.csv("output/modules.csv")
nodemod <- left_join(nodes, modules, by = c("label" = "name"))
idx <- match(V(g)$name, nodemod$id)


# assign modules to vertices
V(g)$label <- nodemod$label[idx]
V(g)$module <- nodemod$module[idx]
V(g)$taxa_roles <- nodemod$taxa_roles[idx]
V(g)$Phylum <- nodemod$Phylum[idx]
V(g)$Abundance <- nodemod$Abundance[idx]
V(g)$deg_mod <- nodemod$degree[idx]
V(g)$bc <- nodemod$betweenness_centrality[idx]
V(g)$cc <- nodemod$closeness_centrality[idx]
V(g)$ec <- nodemod$eigenvector_centrality[idx]


# colour map
mods <- factor(V(g)$module)
num <- as.numeric(sub("M", "", levels(mods))) # make sure modules are in numeric order
mod_levels <- levels(mods)[order(num)]
mods <- factor(V(g)$module, levels = mod_levels)

# pal_vec <- rev(ggsci::pal_d3("category20b")(length(mod_levels)))
pal_vec <- c("#9C9EDEFF", "#DE9ED6FF", "#CEDB9CFF", "#E7969CFF", "#E7CB94FF",
             "#6B6ECFFF", "#CE6DBDFF", "#B5CF6BFF", "#D6616BFF", "#E7BA52FF",
             "#5254A3FF", "#A55194FF", "#8CA252FF", "#AD494AFF", "#BD9E39FF",
             "#393B79FF", "#7B4173FF", "#637939FF", "#843C39FF", "#8C6D31FF")

pal <- setNames(pal_vec, mod_levels)
V(g)$color <- pal[as.character(V(g)$module)]
V(g)$color[is.na(V(g)$color)] <- "#B3B3B3"


# plot
set.seed(42)

lay <- layout_with_fr(g) # set layout

V(g)$size <- 1 + sqrt(degree(g)) * 1.2  # scale nodes by degree

roles_raw <- V(g)$taxa_roles
roles_clean <- ifelse(is.na(roles_raw) | roles_raw == "" | roles_raw == "NA", "Unassigned", roles_raw) # NAs are unassigned roles
V(g)$taxa_roles <- roles_clean

add_triangle_shape <- function() { # add additional shape
  add.vertex.shape(
    "triangle",
    clip = igraph.shape.noclip,
    plot = function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) vertex.color <- vertex.color[v]
      
      vertex.frame.color <- params("vertex", "frame.color")
      if (length(vertex.frame.color) != 1 && !is.null(v)) vertex.frame.color <- vertex.frame.color[v]
      
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) vertex.size <- vertex.size[v]
      
      for (i in seq_along(coords[,1])) {
        x <- coords[i,1]
        y <- coords[i,2]
        s <- vertex.size[i]
        
        x_tri <- c(0, -0.866, 0.866) * s + x
        y_tri <- c(1, -0.5, -0.5) * s + y
        
        polygon(x_tri, y_tri, col=vertex.color[i], border=vertex.frame.color[i])
      }
    }
  )
}

add_triangle_shape()

shape_map <- c("Peripheral nodes" = "circle", # identify hubs and connectors with shapes
               "Connectors" = "triangle",
               "Module hubs" = "square",
               "Unassigned" = "circle")

V(g)$shape <- unname(shape_map[V(g)$taxa_roles])

# V(g)$frame.color <- scales::alpha("black", 0.6) # add outline
V(g)$frame.width <- 0.6

V(g)$frame.width[V(g)$taxa_roles == "Module hubs"] <- 1.2 # make hub outline thicker

# V(g)$color <- alpha(V(g)$color, 0.8) # set transparency

pdf("output/figures/network_plot.pdf", width = 7.09, height = 7.09)

plot(g, layout = lay,
     vertex.label = NA,
     vertex.size = V(g)$size,
     vertex.color = V(g)$color,
     vertex.shape = V(g)$shape,
     edge.width = 0.4, edge.curved = 0.4)

# module colour legend
legend("topleft", #x=-2, y=1,
       legend = mod_levels,
       col = unname(pal[mod_levels]),
       pch = 16, pt.cex = 0.8, bty = "n", title = "Module", ncol = 2, cex = 0.8)

# role shape legend
role_levels <- c("Peripheral nodes", "Connectors", "Module hubs")
role_shapes <- c("circle", "triangle", "square")
shape_map_pch <- c(circle = 21, square = 22, triangle = 24)
legend("bottomleft", #x=-2, y=-0.5,
       legend = role_levels,
       pch = shape_map_pch[role_shapes],
       pt.cex = 1.5, pt.bg = "gray", bty = "n", title = "Role", cex = 0.8)

# degree size legend
deg_vals <- c(1, 50, 100)
size_vals <- 1 + sqrt(deg_vals) * 1.2

legend_cex <- size_vals / max(V(g)$size)  # scale relative to max plotted size
legend_cex <- legend_cex * 3

legend("topright", #x=-2, y=0,
       legend = paste(deg_vals),
       pt.cex = legend_cex, pch = 21, pt.bg = "gray", bty = "n", title = "Degree", cex = 0.8)

dev.off()
