library(tidyverse)
library(ggbeeswarm)

soma_data <- read_tsv(file = "Documents/Data/Somascan/data_phase2_somalogic_curated_20251021.tsv") 
meta_manifest <- read_tsv("Documents/Data/Olink/metadata_olink_soma_phase2.tsv") # Static predefined metadata (for plots mostly)

#source("Documents/Projects/Somalogic-data-exploration/scripts/functions_utility.R")


long_data <- 
  soma_data |> 
  pivot_longer(names_to = "aptamer", values_to = "rfu", cols = -DAid) 

outlier_thresholds <- 
  long_data |>
  group_by(aptamer) |> 
  mutate(outlier_threshold = 6 * mad(rfu),
        median = median(rfu),
      upper_limit = median + outlier_threshold,
      lower_limint = median - outlier_threshold,
    measurements = n()) |> 
  arrange(- outlier_threshold)

outlier_measurements <- 
  outlier_thresholds |> 
  filter(!(rfu > lower_limint & rfu < upper_limit)) |> 
  nrow()

total_measurements <- outlier_thresholds |> nrow()

(outlier_measurements/total_measurements) * 100


outliers_per_aptamer <- 
  outlier_thresholds |> 
  mutate(outlier = ifelse(!(rfu > lower_limint & rfu < upper_limit), "Yes", "No")) |> 
  count(aptamer, outlier) |> 
  filter(outlier == "Yes") |> 
  arrange(-n) 

outliers_per_aptamer |> 
  ggplot(aes(n)) +
  geom_histogram()

outlier_thresholds |> 
  distinct(aptamer, outlier_threshold) |>
  ggplot(aes(outlier_threshold)) +
  geom_histogram() +
    geom_vline(
    aes(xintercept = median(outlier_threshold, na.rm = TRUE)),
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  ggtitle("Distribution of outlier thresholds (6 * MAD)")

outlier_thresholds |> 
  distinct(aptamer, outlier_threshold)

current_aptamers <- "seq.6405.74"

current_outlier_threshold <- 
  outlier_thresholds |>
  distinct(aptamer, median, upper_limit, lower_limint) |> 
  filter(aptamer == current_aptamers)

outlier_thresholds |> 
  mutate(outlier = ifelse(!(rfu > lower_limint & rfu < upper_limit), "Yes", "No")) |> 
  left_join(meta_manifest, by = "DAid") |> 
  filter(aptamer %in% current_aptamers) |> # ""seq.6402.8
 # mutate(Disease = factor(Disease, levels = unique(disease_n$Disease))) |> 
  ggplot(aes(Disease, rfu)) +
  geom_quasirandom(size = 0.5, aes(color = outlier)) +
  geom_boxplot(alpha = 0.2, color = "grey") +
  geom_hline(yintercept = current_outlier_threshold$lower_limint, lty = "dashed", color = "darkgrey") +
  geom_hline(yintercept = current_outlier_threshold$upper_limit, lty = "dashed", color = "darkgrey") +

  #  scale_fill_manual(values = pal_class) +
#  scale_color_manual(values = pal_class) +
#  facet_wrap(~aptamer, scales = "free_y", ncol = 1) +
  theme_hpa(angled = T) #+
#  ggtitle(unique(top_protein_aptamers$TargetFullName))


outlier_thresholds |> 
  mutate(outlier = ifelse(!(rfu > lower_limint & rfu < upper_limit), "Yes", "No")) |> 
  left_join(meta_manifest, by = "DAid") |> 
  filter(aptamer %in% current_aptamers) |> 
# mutate(Disease = factor(Disease, levels = unique(disease_n$Disease))) |> 
  ggplot(aes(Disease, rfu)) +
  
  # shaded outlier zones
  annotate(
    "rect", xmin = -Inf, xmax = Inf, 
    ymin = -Inf, ymax = current_outlier_threshold$lower_limint,
    fill = "lightgrey", alpha = 0.2
  ) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, 
    ymin = current_outlier_threshold$upper_limit, ymax = Inf,
    fill = "lightgrey", alpha = 0.2
  ) +

  geom_quasirandom(size = 0.5, alpha = 0.5, aes(color = outlier)) +
  geom_boxplot(alpha = 0.2, color = "grey") +
  geom_hline(yintercept = current_outlier_threshold$lower_limint, 
             lty = "dashed", color = "darkgrey") +
  geom_hline(yintercept = current_outlier_threshold$upper_limit, 
             lty = "dashed", color = "darkgrey") +
  
  # optional: add fill legend for shaded zones
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.2))
  ) +
  facet_wrap(~aptamer, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("grey20", "darkred")) +
  scale_fill_manual(
    name = "", 
    values = c("lightgrey" = "lightgrey"),
    labels = c("lightgrey" = "Outlier zone")
  ) + 
  ylab("log2(RFU)")
  theme_hpa(angled = TRUE)








theme_hpa <-
  function(angled = F, axis_x = T, axis_x_title = F, axis_y = T, facet_title = T) {
    t <-
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )

    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }

    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    }

    if(axis_x_title == T) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }

    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }