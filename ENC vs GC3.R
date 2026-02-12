#ENC vs GC3  

# Load libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)  # for better labels

# Load data
data <- read_excel("C:/Users/User/OneDrive/Documents/CUB_Frequency.xlsx")

# Calculate theoretical ENC curve
gc3s_theory <- seq(0, 1, by=0.01)
enc_theory <- 2 + gc3s_theory + (29 / (gc3s_theory^2 + (1 - gc3s_theory)^2))

# Create a data frame for shading
curve_df <- data.frame(
  GC3 = gc3s_theory,
  ENC = enc_theory,
  ENC_upper = enc_theory + 2,  # Neutral mutation zone
  ENC_lower = enc_theory - 2   # Selection zone
)

# Basic scatter plot with shapes
p <- ggplot() +
  # Shaded neutral mutation region
  geom_ribbon(data = curve_df, aes(x = GC3, ymin = ENC_lower, ymax = ENC_upper),
              fill = "lightgrey", alpha = 0.5) +
  # Theoretical ENC curve
  geom_line(data = curve_df, aes(x = GC3, y = ENC), 
            color = "black", linetype = "dashed", size=1) +
  # Points with both color & shape
  geom_point(data = data, aes(x = GC3, y = ENC, color = Genus, shape = Genus), size=3) +
  # Label outliers
  geom_text_repel(data = data, aes(x = GC3, y = ENC, 
                                   label = ifelse(ENC < 40 | ENC > 53, Organism, "")),
                  size=2.5, max.overlaps = 20) +
  # Themes and labels
  labs(
    title = " ",
    subtitle = " ",
    x = "GC3s (GC content at 3rd codon position)",
    y = "ENC (Effective Number of Codons)",
    color = "Genus", shape = "Genus"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.major = element_blank(),   # remove gridlines
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15, 18))  # Circle, Triangle, Square, Diamond

# Save image
ggsave(
  filename = "Fig 6. ENC-GC3s plot of coronavirus genomes.tif",
  plot = p,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

# Show plot
print(p)
