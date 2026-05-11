#Comparison across host

# Load libraries
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(broom)
library(ggh4x) 

# Load data
data <- read_excel("D:/Manuscript/Submission/S1 Table.xlsx") %>%
  mutate(
    Genus = factor(Genus, levels = c("Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus", "Deltacoronavirus")),
  )

# Function for clean p-values
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
}


# SUMMARY STATISTICS
summary_stats <- data %>%
  pivot_longer(cols = c(ENC, CAI, GC3s), names_to = "Metric", values_to = "Value") %>%
  group_by(Metric) %>%
  summarise(
    Min = round(min(Value, na.rm = TRUE), 3),
    Max = round(max(Value, na.rm = TRUE), 3),
    Median = round(median(Value, na.rm = TRUE), 3),
    Mean = round(mean(Value, na.rm = TRUE), 3)
  )

print(summary_stats)
write.csv(summary_stats, "Summary_statistics.csv", row.names = FALSE)


# Function for clean p-values
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.6f", p))
}

# SUMMARY STATISTICS (PER GENUS)
summary_stats <- data %>%
  pivot_longer(cols = c(ENC, CAI, GC3s),
               names_to = "Metric",
               values_to = "Value") %>%
  group_by(Genus, Metric) %>%
  summarise(
    n = n(),
    Min = round(min(Value, na.rm = TRUE), 3),
    Max = round(max(Value, na.rm = TRUE), 3),
    Mean = round(mean(Value, na.rm = TRUE), 3),
    Median = round(median(Value, na.rm = TRUE), 3),
    SD = round(sd(Value, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_stats)

write.csv(summary_stats,
          "Summary_statistics_by_genus.csv",
          row.names = FALSE)

# KRUSKAL-WALLIS TEST
kruskal_genus <- map_dfr(c("ENC", "CAI", "GC3s"), ~{
  kruskal.test(as.formula(paste(.x, "~ Genus")), data = data) %>%
    tidy() %>%
    mutate(
      Metric = .x,
      p_label = format_p(p.value)
    )
})

print(kruskal_genus)

# DUNN TEST
dunn_genus <- map_dfr(c("ENC", "CAI", "GC3s"), ~{
  dunn_test(data, as.formula(paste(.x, "~ Genus")), p.adjust.method = "bonferroni") %>%
    mutate(
      Metric = .x,
      p_label = format_p(p.adj)
    )
})

write.csv(dunn_genus, "Dunn_genus_results_clean.csv", row.names = FALSE)

# GENUS PLOT
genus_plot <- data %>%
  pivot_longer(cols = c(ENC, CAI, GC3s), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Genus, y = Value, fill = Genus)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_npg() + # Restored your original colors
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 4) + # Restored p-values
  labs(x = "Genus", y = "Metric Value") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  
  # ADJUST Y-AXIS
  facetted_pos_scales(
    y = list(
      Metric == "GC3s" ~ scale_y_continuous(limits = c(0.1, 0.5)),
      Metric == "CAI" ~ scale_y_continuous(), # Keep original free scale
      Metric == "ENC" ~ scale_y_continuous()  # Keep original free scale
    )
  )

# Display result
print(genus_plot)

# SAVE FIGURES
ggsave("Fig1_Genus.tiff", genus_plot,
       width = 14, height = 7, dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")

