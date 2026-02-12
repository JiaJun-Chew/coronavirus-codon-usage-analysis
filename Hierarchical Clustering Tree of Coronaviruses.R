#Hierarchical Clustering Tree of Coronaviruses

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ape)
library(dendextend)
library(pvclust)
library(RColorBrewer)
library(scales)

# Load data
data <- read_excel("C:/Users/User/OneDrive/Documents/CUB_Frequency.xlsx")

data <- data %>%
  mutate(
    Organism = as.character(Organism),
    Genus = factor(Genus),
    Host  = factor(Host)
  )

# Prepare RSCU matrix
rscu_codons <- grep("^[AUCG]{3}$", names(data), value = TRUE)
exclude_codons <- c("AUG","UGG","UAA","UAG","UGA")
rscu_filtered <- setdiff(rscu_codons, exclude_codons)

rscu_data <- data %>% select(all_of(rscu_filtered))
rscu_data[] <- lapply(rscu_data, as.numeric)
rownames(rscu_data) <- data$Organism

# pvclust bootstrap clustering
set.seed(123)
pv_fit <- pvclust(
  t(rscu_data),
  method.hclust = "ward.D2",
  method.dist   = "euclidean",
  nboot         = 1000
)

hc   <- as.hclust(pv_fit$hclust)
tree <- as.phylo(hc)
tree$tip.label <- data$Organism

# Metadata 
tip_metadata <- data.frame(
  label = data$Organism,
  Host  = data$Host,
  Genus = data$Genus
)

# Define colors & shapes
host_list   <- levels(tip_metadata$Host)
shape_values <- c(16,17,15,18,8,3,4,7,9,10,11,12) 
shape_map    <- setNames(shape_values[1:length(host_list)], host_list)

genus_list  <- levels(tip_metadata$Genus)
genus_colors <- setNames(brewer.pal(max(3, length(genus_list)), "Dark2")[1:length(genus_list)], genus_list)

# Extract AU p-values for internal nodes 
au_pvals <- pv_fit$edges[, "au"]
names(au_pvals) <- 1:length(au_pvals)

# Build ggtree plot
p <- ggtree(tree, layout = "rectangular", size = 0.8) %<+% tip_metadata +
  # Host shapes + Genus colors
  geom_tippoint(aes(color = Genus, shape = Host), size = 3, stroke = 0.9) +
  geom_tiplab(aes(label = label, color = Genus), size = 2.6,
              align = TRUE, offset = 0.02, hjust = 0, fontface = "bold") +
  # AU node labels: red (<95%), black (≥95%)
  geom_nodelab(aes(
    label = ifelse(!is.na(as.numeric(node)),
                   round(au_pvals[as.character(node)], 1), "")),
    hjust = -0.3, vjust = -0.3, size = 2.8,
    color = ifelse(au_pvals < 95, "red", "black")) +
  scale_color_manual(values = genus_colors, name = "Genus") +
  scale_shape_manual(values = shape_map, name = "Host") +
  ggtitle(" ") +
  xlim(0, max(dist(rscu_data)) * 1.0) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 1),
    legend.position = c(0.03, 0.95),
    legend.justification = c(0, 1),
    legend.box = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 5, 10, 5)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3), ncol = 1),
    shape = guide_legend(order = 2, override.aes = list(size = 3), ncol = 1)
  )

# Save image
ggsave(
  filename = "Fig 4. Hierarchical clustering of coronaviruses spike gene based on relative synonymous codon usage profiles.tif",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

print(p)
