#!/usr/bin/env Rscript
library(readxl)
library(ggplot2)

## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output_dir <- "location/to/MOIRE"
dir.create(output_dir, recursive = TRUE)

# FIXME: change to path to metadata
metadata_file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata_file)


# FIXME: change to path of COI summary
coi_summary_file <- "location/to/MOIRE/mhap_COI_summary.tsv"

coi_summary <- read.delim(coi_summary_file)

# FIXME: change to path of effective COI summary
effective_coi_summary_file <- "location/to/MOIRE/mhap_effective_COI_summary.tsv"

effective_coi_summary <- read.delim(effective_coi_summary_file)


# FIXME: change according sample and group column in metadata
metadata_sample_column <- "ID"
metadata_group_column <- "Year"

coi_summary <- merge(
  coi_summary,
  metadata,
  by.x = "sample_id",
  by.y = metadata_sample_column
)

effective_coi_summary <- merge(
  effective_coi_summary,
  metadata,
  by.x = "sample_id",
  by.y = metadata_sample_column
)

## STOP and check the number of samples remaining after the previous steps

## setup colours here
z <- as.factor(metadata[[metadata_group_column]])

# FIXME: select the colour palette according to the number of groups
# The "Tableau 10" colour palette can accommodate up to 10 groups (default)
# There are also "ggplot2" and other palettes from `palette.pals()`
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
nlevels(z)

cols <- palette.colors(palette = "Tableau 10")

palette <- setNames(cols[1:nlevels(z)], levels(z))


coi_plot_file <- paste0(output_dir, "/", "mhap_COI.pdf")

# FIXME: adjust plotting template if necessary
pdf(file = coi_plot_file, width = 4, height = 4)

ggplot(
  coi_summary,
  aes(x = as.factor(.data[[metadata_group_column]]), y = post_coi_mean)
) +
  geom_jitter(
    aes(colour = as.factor(.data[[metadata_group_column]])),
    width = 0.3,
    height = 0,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  xlab(metadata_group_column) +
  ylab("COI") +
  theme_classic() +
  scale_y_continuous(breaks = 1:ceiling(max(coi_summary[["post_coi_mean"]]))) +
  scale_colour_manual(breaks = levels(z), values = palette)

dev.off()


polyclonal_proportion_file <- paste0(
  output_dir,
  "/",
  "mhap_polyclonal_proportion.tsv"
)

# estimate polyclonal proportion by taking the mean of polyclonal probabilities
polyclonal_proportion <- aggregate(
  coi_summary["prob_polyclonal"],
  coi_summary[metadata_group_column],
  mean
)
names(polyclonal_proportion)[2] <- "polyclonal_proportion"

write.table(
  polyclonal_proportion,
  polyclonal_proportion_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

polyclonal_proportion_plot_file <- paste0(
  output_dir,
  "/",
  "mhap_polyclonal_proportion.pdf"
)

# FIXME: adjust plotting template if necessary
pdf(file = polyclonal_proportion_plot_file, width = 4, height = 4)

ggplot(
  polyclonal_proportion,
  aes(x = as.factor(.data[[metadata_group_column]]), y = polyclonal_proportion)
) +
  geom_col(
    aes(fill = as.factor(.data[[metadata_group_column]])),
    show.legend = FALSE
  ) +
  ylim(0, 1) +
  xlab(metadata_group_column) +
  ylab("Estimated polyclonal proportion") +
  theme_classic() +
  scale_fill_manual(breaks = levels(z), values = palette)

dev.off()


effective_coi_plot_file <- paste0(output_dir, "/", "mhap_effective_COI.pdf")

# FIXME: adjust plotting template if necessary
pdf(file = effective_coi_plot_file, width = 4, height = 4)

ggplot(
  effective_coi_summary,
  aes(
    x = as.factor(.data[[metadata_group_column]]),
    y = post_effective_coi_mean
  )
) +
  geom_jitter(
    aes(colour = as.factor(.data[[metadata_group_column]])),
    width = 0.3,
    height = 0,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  xlab(metadata_group_column) +
  ylab("eMOI") +
  theme_classic() +
  scale_y_continuous(
    breaks = 1:ceiling(max(effective_coi_summary[["post_effective_coi_mean"]]))
  ) +
  scale_colour_manual(breaks = levels(z), values = palette)

dev.off()
