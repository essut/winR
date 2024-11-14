#!/usr/bin/env Rscript
library(readxl)
library(ggplot2)

## Change FIXME lines to the appropriate parameters

# FIXME: change to path to metadata
metadata_file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata_file)


# FIXME: change to path of COI summary
coi_summary_file <- "location/to/mhap_COI_summary.tsv"

coi_summary <- read.delim(coi_summary_file)

# FIXME: change to path of effective COI summary
effective_coi_summary_file <- "location/to/mhap_effective_COI_summary.tsv"

effective_coi_summary <- read.delim(effective_coi_summary_file)


# FIXME: change according sample and group column in metadata
metadata_sample_column <- "ID"
metadata_group_column <- "Year"

coi_summary <-
  merge(
    coi_summary,
    metadata,
    by.x = "sample_id",
    by.y = metadata_sample_column
  )

effective_coi_summary <-
  merge(
    effective_coi_summary,
    metadata,
    by.x = "sample_id",
    by.y = metadata_sample_column
  )

## STOP and check the number of samples remaining after the previous steps


# FIXME: adjust name and size (in inches) of COI plot
coi_plot_file <- "location/to/mhap_COI.pdf"

pdf(file = coi_plot_file, width = 7, height = 7)

# FIXME: COI plot framework, adjust metadata as necessary
ggplot(
  coi_summary,
  aes(x = as.factor(.data[[metadata_group_column]]), y = post_coi_mean)
) +
  geom_boxplot() +
  xlab(metadata_group_column) +
  ylab("COI") +
  theme_classic() +
  scale_y_continuous(breaks = 1:ceiling(max(coi_summary[["post_coi_mean"]])))

dev.off()


# FIXME: adjust name and size (in inches) of polyclonal prevalence plot
polyclonal_prevalence_plot_file <- "location/to/mhap_polyclonal_prevalence.pdf"

pdf(file = polyclonal_prevalence_plot_file, width = 7, height = 7)

# FIXME: polyclonal prevalence plot framework, adjust metadata as necessary
ggplot(
  coi_summary,
  aes(x = as.factor(.data[[metadata_group_column]]), fill = prob_polyclonal > 0.5)
) +
  geom_bar(position = "fill") +
  xlab(metadata_group_column) +
  ylab("Proportion") +
  labs(fill = "Polyclonal") +
  theme_classic() +
  theme(legend.position = "bottom")

dev.off()


# FIXME: adjust name and size (in inches) of COI plot
effective_coi_plot_file <- "location/to/mhap_effective_COI.pdf"

pdf(file = effective_coi_plot_file, width = 7, height = 7)

# FIXME: COI plot framework, adjust metadata as necessary
ggplot(
  effective_coi_summary,
  aes(x = as.factor(.data[[metadata_group_column]]), y = post_effective_coi_mean)
) +
  geom_boxplot() +
  xlab(metadata_group_column) +
  ylab("eMOI") +
  theme_classic() +
  scale_y_continuous(breaks = 1:ceiling(max(effective_coi_summary[["post_effective_coi_mean"]])))

dev.off()
