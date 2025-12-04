# Script: plot_atac_p53_overlap.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("workflow/scripts/plot_atac_p53_overlap.rda"))
message("Saved Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING PACKAGES ==========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})


### LOADING FILES =============================================================

# atac peaks
atac_peaks <- read.table(snakemake@input$atac, header = FALSE)
colnames(atac_peaks) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
atac_peaks <- makeGRangesFromDataFrame(atac_peaks, keep.extra.columns = TRUE)

# p53 full peaks
p53_full_original <- read.table(snakemake@input$p53_full, header = FALSE)
colnames(p53_full_original) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
p53_full_original <- makeGRangesFromDataFrame(p53_full, keep.extra.columns = TRUE)

# p53 IDR peaks
p53_idr_original <- read.table(snakemake@input$p53_idr, header = FALSE)
colnames(p53_idr_original) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
p53_idr_original <- makeGRangesFromDataFrame(p53_idr, keep.extra.columns = TRUE)


### PARAMETERS ================================================================

slop_distance <- 500  # Modify this to extend peaks


### OVERLAP ANALYSIS ==========================================================
# Extend peaks by slop distance if needed
extend_peaks <- function(peaks, slop) {
  if (slop > 0) {
    peaks <- resize(peaks, width = width(peaks) + 2*slop, fix = "center")
  }
  return(peaks)
}

p53_full <- extend_peaks(p53_full_original, slop_distance)
p53_idr <- extend_peaks(p53_idr_original, slop_distance)

# Find overlaps with ATAC peaks
overlaps_full <- suppressWarnings(findOverlaps(p53_full, atac_peaks))
overlaps_idr <- suppressWarnings(findOverlaps(p53_idr, atac_peaks))

# Calculate overlap statistics
p53_full_with_overlap <- length(unique(queryHits(overlaps_full)))
p53_idr_with_overlap <- length(unique(queryHits(overlaps_idr)))

overlap_stats <- data.frame(
  dataset = c("p53 Full", "p53 IDR"),
  total_peaks = c(length(p53_full), length(p53_idr)),
  overlapping = c(p53_full_with_overlap, p53_idr_with_overlap),
  not_overlapping = c(length(p53_full) - p53_full_with_overlap, 
                      length(p53_idr) - p53_idr_with_overlap)
) %>%
  mutate(
    overlap_percent = (overlapping / total_peaks) * 100,
    no_overlap_percent = (not_overlapping / total_peaks) * 100,
    overlap_prop = overlapping / total_peaks,
    no_overlap_prop = not_overlapping / total_peaks
  )

### DISTANCE ANALYSIS =========================================================
# Calculate distance to nearest ATAC peak for all p53 peaks
dist_full <- suppressWarnings(distanceToNearest(p53_full, atac_peaks))
dist_idr <- suppressWarnings(distanceToNearest(p53_idr, atac_peaks))

dist_data <- bind_rows(
  data.frame(
    dataset = "p53 Full",
    distance = mcols(dist_full)$distance,
    overlaps = mcols(dist_full)$distance == 0
  ),
  data.frame(
    dataset = "p53 IDR",
    distance = mcols(dist_idr)$distance,
    overlaps = mcols(dist_idr)$distance == 0
  )
)

### PLOTTING ==================================================================
# Color palette - more distinct blues
colors_blue <- c("#08519c", "#3182bd", "#6baed6", "#c6dbef")

# 1. Overlap proportion plot
plot_overlap <- overlap_stats %>%
  pivot_longer(cols = c(overlap_prop, no_overlap_prop),
               names_to = "status", values_to = "proportion") %>%
  mutate(status = factor(status, levels = c("overlap_prop", "no_overlap_prop"),
                         labels = c("Overlap", "No Overlap")),
         dataset_label = paste0(dataset, "\n(n=", total_peaks, ")")) %>%
  ggplot(aes(x = dataset_label, y = proportion, fill = status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Overlap" = colors_blue[2], "No Overlap" = colors_blue[4])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "p53 Peak Overlap with ATAC-seq",
       subtitle = if(slop_distance > 0) paste0("Slop: ", slop_distance, " bp") else "No extension",
       x = NULL, y = "Proportion of p53 peaks", fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "bold"))

# 2. Distance violin plot - non-overlapping peaks only
summary_stats <- dist_data %>%
  filter(!overlaps) %>%
  group_by(dataset) %>%
  summarize(
    median_dist = median(distance),
    n = n(),
    .groups = "drop"
  )

plot_distance_violin <- dist_data %>%
  filter(!overlaps) %>%
  mutate(dataset_label = case_when(
    dataset == "p53 Full" ~ paste0("p53 Full\n(n=", summary_stats$n[summary_stats$dataset == "p53 Full"], ")"),
    dataset == "p53 IDR" ~ paste0("p53 IDR\n(n=", summary_stats$n[summary_stats$dataset == "p53 IDR"], ")")
  )) %>%
  ggplot(aes(x = dataset_label, y = distance, fill = dataset)) +
  geom_violin(alpha = 0.6, color = colors_blue[1], size = 0.8) +
  geom_segment(data = summary_stats %>%
                 mutate(dataset_label = case_when(
                   dataset == "p53 Full" ~ paste0("p53 Full\n(n=", n, ")"),
                   dataset == "p53 IDR" ~ paste0("p53 IDR\n(n=", n, ")")
                 ),
                 x_start = as.numeric(factor(dataset_label)) - 0.4,
                 x_end = as.numeric(factor(dataset_label)) + 0.4),
               aes(x = x_start, xend = x_end, y = median_dist, yend = median_dist, color = dataset),
               linetype = "dashed", size = 0.5, inherit.aes = FALSE) +
  geom_text(data = summary_stats %>%
              mutate(dataset_label = case_when(
                dataset == "p53 Full" ~ paste0("p53 Full\n(n=", n, ")"),
                dataset == "p53 IDR" ~ paste0("p53 IDR\n(n=", n, ")")
              )),
            aes(x = dataset_label, y = median_dist * 1.5, 
                label = paste0("Median: ", scales::comma(median_dist), " bp")),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_y_log10(labels = scales::comma, 
                breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
  scale_fill_manual(values = c("p53 Full" = colors_blue[2], "p53 IDR" = colors_blue[3])) +
  scale_color_manual(values = c("p53 Full" = colors_blue[1], "p53 IDR" = colors_blue[1])) +
  labs(title = "Distance to Nearest ATAC Peak",
       subtitle = "Non-overlapping peaks only",
       x = NULL, y = "Distance (bp, log scale)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(face = "bold")) +
  guides(color = "none")

# 3. Distance histogram for non-overlapping peaks
non_overlap_stats <- dist_data %>%
  filter(!overlaps) %>%
  group_by(dataset) %>%
  summarize(
    n = n(),
    median_dist = median(distance),
    mean_dist = mean(distance),
    .groups = "drop"
  )

plot_distance_hist <- dist_data %>%
  filter(!overlaps) %>%
  ggplot(aes(x = distance, fill = dataset)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity", color = "white", size = 0.3) +
  geom_vline(data = non_overlap_stats,
             aes(xintercept = median_dist, color = dataset),
             linetype = "dashed", size = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = c("p53 Full" = colors_blue[1], "p53 IDR" = colors_blue[3])) +
  scale_color_manual(values = c("p53 Full" = colors_blue[1], "p53 IDR" = colors_blue[3])) +
  labs(title = "Distance Distribution (Non-overlapping Peaks)",
       subtitle = "Dashed lines show median distances",
       x = "Distance to nearest ATAC peak (bp, log scale)", 
       y = "Count", fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black")) +
  guides(color = "none")

# Combine plots
combined_plot <- (plot_overlap | plot_distance_violin) / plot_distance_hist +
  plot_annotation(
    title = "ATAC-seq and p53 ChIP-seq Peak Overlap Analysis",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  )

### SAVE OUTPUT ===============================================================

# message("Saving output files")
# 
# # Save the combined plot
# ggsave(snakemake@output$plot, plot = combined_plot, width = 14, height = 16, device = "pdf")
# https://ondemand.sherlock.stanford.edu/rnode/sh02-10n04.int/57438/graphics/d0bf28a7-7beb-4df9-a34d-f84935598f48.png
# # Save overlap statistics
# write_csv(overlap_stats, snakemake@output$overlap_stats)
# 
# # Save distance data
# write_csv(dist_data, snakemake@output$distance_data)
# 
# # Save summary statistics
# summary_stats <- dist_data %>%
#   group_by(dataset) %>%
#   summarize(
#     total_peaks = n(),
#     overlapping = sum(overlaps),
#     not_overlapping = sum(!overlaps),
#     overlap_percent = (overlapping / total_peaks) * 100,
#     median_distance = median(distance),
#     mean_distance = mean(distance),
#     median_distance_nonoverlap = median(distance[!overlaps]),
#     mean_distance_nonoverlap = mean(distance[!overlaps])
#   )
# 
# write_csv(summary_stats, snakemake@output$summary_stats)
# 
# message("Analysis complete!")
# message("Output files saved:")
# message("  - Plot: ", snakemake@output$plot)
# message("  - Overlap stats: ", snakemake@output$overlap_stats)
# message("  - Distance data: ", snakemake@output$distance_data)
# message("  - Summary stats: ", snakemake@output$summary_stats)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)