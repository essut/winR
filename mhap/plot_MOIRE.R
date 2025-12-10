#!/usr/bin/env Rscript
library(readxl)
library(ggplot2)

## Load all functions by running everything from line 6-48
calculate_polyclonal_prevalence <- function(
  polyclonal_status,
  metadata_group_column
) {
  n_polyclonal <-
    aggregate(
      polyclonal_status["is_polyclonal"],
      polyclonal_status[metadata_group_column],
      sum
    )
  n_total <-
    aggregate(
      polyclonal_status["is_polyclonal"],
      polyclonal_status[metadata_group_column],
      length
    )

  polyclonal_prevalence <-
    merge(n_polyclonal, n_total, by = metadata_group_column, sort = FALSE)
  names(polyclonal_prevalence)[2:3] <- c("n_polyclonal", "n_total")

  htests <-
    tapply(
      polyclonal_prevalence,
      polyclonal_prevalence[[metadata_group_column]],
      function(x) binom.test(x[["n_polyclonal"]], x[["n_total"]])
    )
  
  polyclonal_prevalence <-
    data.frame(
      metadata_group_column = names(htests),
      n_polyclonal = vapply(htests, function(x) x[["statistic"]], numeric(1)),
      n_total = vapply(htests, function(x) x[["parameter"]], numeric(1)),
      pc_polyclonal = 100 * vapply(htests, function(x) x[["estimate"]], numeric(1)),
      CI_lower_95 = 100 * vapply(htests, function(x) x[["conf.int"]], numeric(2))[1, ],
      CI_upper_95 = 100 * vapply(htests, function(x) x[["conf.int"]], numeric(2))[2, ]
    )
  names(polyclonal_prevalence)[
    names(polyclonal_prevalence) %in% "metadata_group_column"
  ] <- metadata_group_column

  return(polyclonal_prevalence)
}


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

# FIXME: change to path of polyclonal status
polyclonal_status_file <- "location/to/MOIRE/mhap_polyclonal_status.tsv"

polyclonal_status <- read.delim(polyclonal_status_file)

# FIXME: change to path of effective COI summary
effective_coi_summary_file <- "location/to/MOIRE/mhap_effective_COI_summary.tsv"

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

polyclonal_status <-
  merge(
    polyclonal_status,
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

coi_plot_file <- paste0(output_dir, "/", "mhap_COI.pdf")

# FIXME: adjust plotting template if necessary
pdf(file = coi_plot_file, width = 4, height = 4)

ggplot(
  coi_summary,
  aes(x = as.factor(.data[[metadata_group_column]]), y = post_coi_mean)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  xlab(metadata_group_column) +
  ylab("COI") +
  theme_classic() +
  scale_y_continuous(breaks = 1:ceiling(max(coi_summary[["post_coi_mean"]])))

dev.off()


polyclonal_prevalence_file <-
  paste0(output_dir, "/", "mhap_polyclonal_prevalence.tsv")

polyclonal_prevalence <-
  calculate_polyclonal_prevalence(polyclonal_status, metadata_group_column)
write.table(
  polyclonal_prevalence,
  polyclonal_prevalence_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

polyclonal_prevalence_plot_file <-
  paste0(output_dir, "/", "mhap_polyclonal_prevalence.pdf")

# FIXME: adjust plotting template if necessary
pdf(file = polyclonal_prevalence_plot_file, width = 4, height = 4)

ggplot(
  polyclonal_prevalence,
  aes(
    x = as.factor(.data[[metadata_group_column]]),
    y = pc_polyclonal,
    ymin = CI_lower_95,
    ymax = CI_upper_95
  )
) +
  geom_col() +
  geom_errorbar() +
  ylim(0, 100) +
  xlab(metadata_group_column) +
  ylab("% polyclonal") +
  theme_classic()

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
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  xlab(metadata_group_column) +
  ylab("eMOI") +
  theme_classic() +
  scale_y_continuous(
    breaks = 1:ceiling(max(effective_coi_summary[["post_effective_coi_mean"]]))
  )

dev.off()
