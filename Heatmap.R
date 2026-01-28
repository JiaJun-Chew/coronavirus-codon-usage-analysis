#Heatmap

# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Load data
df <- read_excel("C:/Users/User/OneDrive/Documents/CUB_Frequency.xlsx")

# Select RSCU codon columns, exclude AUG, UGG, stop codons
rscu_cols <- grep("^[AUCG]{3}$", names(df), value = TRUE)
exclude_codons <- c("AUG", "UGG", "UAA", "UAG", "UGA")
rscu_cols <- setdiff(rscu_cols, exclude_codons)

rscu_data <- df %>% select(all_of(rscu_cols))
rscu_data[] <- lapply(rscu_data, as.numeric)
rownames(rscu_data) <- df$Organism

# Define codon → amino acid mapping (standard genetic code, single-letter AA code)
codon2aa <- c(
  "GCU"="A","GCC"="A","GCA"="A","GCG"="A",
  "UGU"="C","UGC"="C",
  "GAU"="D","GAC"="D",
  "GAA"="E","GAG"="E",
  "UUU"="F","UUC"="F",
  "GGU"="G","GGC"="G","GGA"="G","GGG"="G",
  "CAU"="H","CAC"="H",
  "AUA"="I","AUU"="I","AUC"="I",
  "AAA"="K","AAG"="K",
  "UUA"="L","UUG"="L","CUU"="L","CUC"="L","CUA"="L","CUG"="L",
  "AUG"="M", # excluded
  "AAU"="N","AAC"="N",
  "CCU"="P","CCC"="P","CCA"="P","CCG"="P",
  "CAA"="Q","CAG"="Q",
  "CGU"="R","CGC"="R","CGA"="R","CGG"="R","AGA"="R","AGG"="R",
  "UCU"="S","UCC"="S","UCA"="S","UCG"="S","AGU"="S","AGC"="S",
  "ACU"="T","ACC"="T","ACA"="T","ACG"="T",
  "GUU"="V","GUC"="V","GUA"="V","GUG"="V",
  "UGG"="W", # excluded
  "UAU"="Y","UAC"="Y"
)

# Rename codons as "Codon(AA)" for heatmap columns
colnames(rscu_data) <- paste0(rscu_cols, "(", codon2aa[rscu_cols], ")")

# Plot heatmap with hierarchical clustering
heatmap_plot <- pheatmap(
  as.matrix(rscu_data),
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  main = " ",
  fontsize_col = 8,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  border_color = NA
)

# Save image
ggsave(
  filename = "Figure 5. Heatmap of relative synonymous codon usage values across coronavirus genomes.tif",
  plot = heatmap_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)
