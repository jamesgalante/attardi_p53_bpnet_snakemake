# Script: plot_p53_strand_bias.R
# Analyze strand bias in p53 ChIP-seq peaks with ATAC-seq signal

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("workflow/scripts/plot_p53_strand_bias.rda"))
message("Saved Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING PACKAGES ==========================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})


### LOADING FILES =============================================================

message("Loading input files")

# Load peaks
message("Loading p53 full peaks")
p53_full <- read.table(snakemake@input$p53_full, header = FALSE)
colnames(p53_full) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
p53_full <- makeGRangesFromDataFrame(p53_full, keep.extra.columns = TRUE)

message("Loading p53 IDR peaks")
p53_idr <- read.table(snakemake@input$p53_idr, header = FALSE)
colnames(p53_idr) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
p53_idr <- makeGRangesFromDataFrame(p53_idr, keep.extra.columns = TRUE)

# Load bigWig files
message("Loading bigWig files")
bw_p53_plus <- import(snakemake@input$p53_plus, as = "RleList")
bw_p53_minus <- import(snakemake@input$p53_minus, as = "RleList")
bw_atac <- import(snakemake@input$atac, as = "RleList")
bw_control_plus <- import(snakemake@input$control_plus, as = "RleList")
bw_control_minus <- import(snakemake@input$control_minus, as = "RleList")


### EXTRACT SIGNAL ============================================================

message("Extracting signals from bigWigs")

# Function to extract mean signal from bigWig
extract_signal <- function(bw, peaks, func) {
  scores <- numeric(length(peaks))
  for (i in seq_along(peaks)) {
    chr <- as.character(seqnames(peaks[i]))
    if (chr %in% names(bw)) {
      region <- bw[[chr]][start(peaks[i]):end(peaks[i])]
      scores[i] <- func(region, na.rm = TRUE)
    } else {
      scores[i] <- 0
    }
  }
  return(scores)
}

# Extract signals for p53 Full peaks
message("Extracting signals for p53 Full peaks")
signal_full <- data.frame(
  peak_id = paste0(seqnames(p53_full), ":", start(p53_full), "-", end(p53_full)),
  chr = as.character(seqnames(p53_full)),
  start = start(p53_full),
  end = end(p53_full),
  p53_plus = extract_signal(bw_p53_plus, p53_full, sum),
  p53_minus = extract_signal(bw_p53_minus, p53_full, sum),
  atac = extract_signal(bw_atac, p53_full, mean),
  control_plus = extract_signal(bw_control_plus, p53_full, sum),
  control_minus = extract_signal(bw_control_minus, p53_full, sum)
) %>%
  mutate(
    plus_ratio = (p53_plus + 1) / (control_plus + 1),
    minus_ratio = (p53_minus + 1) / (control_minus + 1),
    total_p53 = p53_plus + p53_minus,
    total_control = control_plus + control_minus,
    p53_enrichment = (total_p53 + 1) / (total_control + 1)
  )

# Extract signals for p53 IDR peaks
message("Extracting signals for p53 IDR peaks")
signal_idr <- data.frame(
  peak_id = paste0(seqnames(p53_idr), ":", start(p53_idr), "-", end(p53_idr)),
  chr = as.character(seqnames(p53_idr)),
  start = start(p53_idr),
  end = end(p53_idr),
  p53_plus = extract_signal(bw_p53_plus, p53_idr, sum),
  p53_minus = extract_signal(bw_p53_minus, p53_idr, sum),
  atac = extract_signal(bw_atac, p53_idr, mean),
  control_plus = extract_signal(bw_control_plus, p53_idr, sum),
  control_minus = extract_signal(bw_control_minus, p53_idr, sum)
) %>%
  mutate(
    plus_ratio = (p53_plus + 1) / (control_plus + 1),
    minus_ratio = (p53_minus + 1) / (control_minus + 1),
    total_p53 = p53_plus + p53_minus,
    total_control = control_plus + control_minus,
    p53_enrichment = (total_p53 + 1) / (total_control + 1)
)


### PLOTTING ==================================================================

message("Creating plots")

# Plot p53 signal + by - strand
plot_full <- signal_full %>%
  ggplot(aes(x = plus_ratio, y = minus_ratio, color = atac)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_viridis_c(option = "plasma", trans = "log10") +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray40", size = 0.5) +
  labs(title = "Strand Bias Analysis - p53 Full",
       subtitle = paste0("n = ", nrow(signal_full), " peaks"),
       x = "Exp Counts in + strand / Control counts in + strand",
       y = "Exp Counts in - strand / Control counts in - strand",
       color = "ATAC signal\n(log10)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"),
        legend.position = "right")

plot_idr <- signal_idr %>%
  ggplot(aes(x = plus_ratio, y = minus_ratio, color = atac)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_viridis_c(option = "plasma", trans = "log10") +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray40", size = 0.5) +
  labs(title = "Strand Bias Analysis - p53 IDR",
       subtitle = paste0("n = ", nrow(signal_idr), " peaks"),
       x = "Exp Counts in + strand / Control counts in + strand",
       y = "Exp Counts in - strand / Control counts in - strand",
       color = "ATAC signal\n(log10)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"),
        legend.position = "right")

# Plot p53 vs ATAC
plot_p53_atac_full <- signal_full %>%
  ggplot(aes(x = atac, y = total_p53)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#1976D2") +
  geom_smooth(method = "lm", color = "#D32F2F", se = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "p53 Signal vs ATAC Signal - p53 Full",
       subtitle = paste0("n = ", nrow(signal_full), " peaks"),
       x = "Mean ATAC signal (log10)",
       y = "Total p53 signal (log10)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"))

plot_p53_atac_idr <- signal_idr %>%
  ggplot(aes(x = atac, y = total_p53)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#1976D2") +
  geom_smooth(method = "lm", color = "#D32F2F", se = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "p53 Signal vs ATAC Signal - p53 IDR",
       subtitle = paste0("n = ", nrow(signal_idr), " peaks"),
       x = "Mean ATAC signal (log10)",
       y = "Total p53 signal (log10)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text = element_text(color = "black"))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")

ggsave(snakemake@output$plot_full, plot = plot_full, 
       width = 10, height = 8, dpi = 300)

ggsave(snakemake@output$plot_idr, plot = plot_idr, 
       width = 10, height = 8, dpi = 300)

write_csv(signal_full, snakemake@output$signal_data_full)
write_csv(signal_idr, snakemake@output$signal_data_idr)

# Summary statistics
summary_full <- signal_full %>%
  summarize(
    dataset = "p53 Full",
    total_peaks = n(),
    mean_p53_plus = mean(p53_plus),
    mean_p53_minus = mean(p53_minus),
    mean_atac = mean(atac),
    mean_control_plus = mean(control_plus),
    mean_control_minus = mean(control_minus),
    median_exp_ratio = median(exp_ratio),
    median_control_ratio = median(control_ratio)
  )

summary_idr <- signal_idr %>%
  summarize(
    dataset = "p53 IDR",
    total_peaks = n(),
    mean_p53_plus = mean(p53_plus),
    mean_p53_minus = mean(p53_minus),
    mean_atac = mean(atac),
    mean_control_plus = mean(control_plus),
    mean_control_minus = mean(control_minus),
    median_exp_ratio = median(exp_ratio),
    median_control_ratio = median(control_ratio)
  )

summary_stats <- bind_rows(summary_full, summary_idr)
write_csv(summary_stats, snakemake@output$summary_stats)

message("Output files saved:")
message("  - Plot (Full): ", snakemake@output$plot_full)
message("  - Plot (IDR): ", snakemake@output$plot_idr)
message("  - Signal data (Full): ", snakemake@output$signal_data_full)
message("  - Signal data (IDR): ", snakemake@output$signal_data_idr)
message("  - Summary stats: ", snakemake@output$summary_stats)


### CLEAN UP ==================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)