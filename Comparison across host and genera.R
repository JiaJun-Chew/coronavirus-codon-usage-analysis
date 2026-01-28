#Comparison across host and genera

library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(ggsci)
library(ggthemes)

# Read the Excel file
data <- read_excel("C:/Users/User/OneDrive/Documents/CUB_Frequency.xlsx") %>%
  mutate(
    Genus = factor(Genus, levels = c("Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus", "Deltacoronavirus")),
    Host = factor(Host, levels = unique(Host)) # Automatically include all hosts
  )

# Kruskal-Wallis for each metric by Genus
kruskal_genus <- map_dfr(c("ENC", "CAI", "GC3"), ~{
  kruskal.test(as.formula(paste(.x, "~ Genus")), data = data) %>%
    broom::tidy() %>%
    mutate(Metric = .x)
})

# Dunn's post-hoc for each metric by Genus
dunn_genus <- map_dfr(c("ENC", "CAI", "GC3"), ~{
  dunn_test(data, as.formula(paste(.x, "~ Genus")), p.adjust.method = "bonferroni") %>%
    mutate(Metric = .x)
})

# Kruskal-Wallis for each metric by Host
kruskal_host <- map_dfr(c("ENC", "CAI", "GC3"), ~{
  kruskal.test(as.formula(paste(.x, "~ Host")), data = data) %>%
    broom::tidy() %>%
    mutate(Metric = .x)
})

# Dunn's post-hoc for each metric by Host
dunn_host <- map_dfr(c("ENC", "CAI", "GC3"), ~{
  dunn_test(data, as.formula(paste(.x, "~ Host")), p.adjust.method = "bonferroni") %>%
    mutate(Metric = .x)
})

genus_plot <- data %>%
  pivot_longer(cols = c(ENC, CAI, GC3), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Genus, y = Value, fill = Genus)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black", lwd = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, shape = 21, color = "black") +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_npg() +
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.y.npc = "top", vjust = -3.0, size = 4) +
  labs(title = " ",
       subtitle = " ",
       x = "Genus", y = "Metric Value") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

host_plot <- data %>%
  pivot_longer(cols = c(ENC, CAI, GC3), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Host, y = Value)) +
  # Neutral boxplots for cleaner Host comparison
  geom_boxplot(fill = "grey90", color = "black", lwd = 0.5, outlier.shape = NA, alpha = 0.7) +
  # Genus as point color (main biological signal)
  geom_jitter(aes(color = Genus), width = 0.2, size = 2.5, alpha = 0.8) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_color_npg() +  # color by Genus only
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.y.npc = "top", vjust = -3.5, size = 4) +
  labs(title = " ",
       subtitle = " ",
       x = "Host Species", y = "Metric Value", color = "Genus") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Save refined host plot
ggsave(
  filename = "Figure 1. Boxplots of codon usage bias across coronavirus genera.tif",
  plot = genus_plot,
  width = 14,
  height = 7,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = "Figure 2. Comparison of codon usage bias metrics across coronavirus host species.tif",
  plot = host_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


print(kruskal_genus)
print(dunn_genus)
print(kruskal_host)
print(dunn_host)

