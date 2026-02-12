#Correspondence analysis

# Load required libraries
library(readxl)
library(dplyr)
library(FactoMineR)
library(ggplot2)

# 1. Load  data
df <- read_excel("C:/Users/User/OneDrive/Documents/CUB_Frequency.xlsx")

# 2. Select RSCU columns and exclude AUG, UGG, and stop codons
rscu_cols <- grep("^[AUCG]{3}$", names(df), value = TRUE)
exclude_codons <- c("AUG", "UGG", "UAA", "UAG", "UGA")
rscu_cols <- setdiff(rscu_cols, exclude_codons)

rscu_data <- df %>% select(all_of(rscu_cols))
rscu_data[] <- lapply(rscu_data, as.numeric)

# 3. Correspondence Analysis (COA)
coa_res <- CA(rscu_data, graph = FALSE)

# 4. Extract individual coordinates (first 2 dimensions)
ind_coords <- as.data.frame(coa_res$row$coord[, 1:2])
colnames(ind_coords) <- c("Dim1", "Dim2")
ind_coords$Organism <- df$Organism
ind_coords$Genus <- df$Genus
ind_coords$Host <- df$Host

# 5. Explained variance
var1 <- round(coa_res$eig[1, 2], 1)
var2 <- round(coa_res$eig[2, 2], 1)

# 6. Plot COA results
p_coa <- ggplot(ind_coords, aes(x = Dim1, y = Dim2, color = Genus)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = Genus, color = Genus),
    linetype = 2,      # Dotted
    level = 0.95,
    size = 1.1,
    alpha = 1
  ) +
  labs(
    title = " ",
    x = paste0("Dimension 1 (", var1, "% variance)"),
    y = paste0("Dimension 2 (", var2, "% variance)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save figure
ggsave(
  filename = "Fig 3. Two-dimensional correspondence analysis (CA) of relative synonymous codon usage (RSCU) values across coronavirus genomes.tif",
  plot = p_coa,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

# Print to screen
print(p_coa)
