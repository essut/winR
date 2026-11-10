#!/usr/bin/env Rscript
library(readxl)
library(ggplot2)

## Load all functions by running everything from line 6-88
decompose.loci <- function(dlong) {
  chrom.pos <- strsplit(dlong[["locus"]], ":")
  chrom <- vapply(chrom.pos, "[[", character(1), 1)
  pos <- vapply(chrom.pos, "[[", character(1), 2)
  marker <- vapply(chrom.pos, "[[", character(1), 3)

  start.end <- strsplit(pos, "-")
  start <- as.numeric(vapply(start.end, "[[", character(1), 1))
  end <- as.numeric(vapply(start.end, "[[", character(1), 2))

  data.frame(
    chromosome = chrom,
    start.pos = start,
    end.pos = end,
    marker = marker,
    dlong
  )
}


## assumes chromosomes have been given order
sort.dlong <- function(dlong) {
  sorted.dlong <- dlong[
    do.call(order, dlong[, c("chromosome", "start.pos", "end.pos", "allele")]),
  ]

  locus.order <- unique(sorted.dlong[["locus"]])
  sorted.dlong[["locus"]] <- factor(
    sorted.dlong[["locus"]],
    levels = locus.order
  )

  marker.order <- unique(sorted.dlong[["marker"]])
  sorted.dlong[["marker"]] <- factor(
    sorted.dlong[["marker"]],
    levels = marker.order
  )

  sorted.dlong
}


calculate.prop <- function(sorted.dlong) {
  dlong.count <- aggregate(count ~ sample_id + locus, sorted.dlong, sum)
  names(dlong.count)[length(dlong.count)] <- "total"

  sorted.dlong.wprop <- merge(sorted.dlong, dlong.count, sort = FALSE)
  sorted.dlong.wprop[["prop"]] <- sorted.dlong.wprop[["count"]] /
    sorted.dlong.wprop[["total"]]

  sorted.dlong.wprop
}


## calculating heterozygosity based on https://www.uwyo.edu/dbmcd/molmark/lect04/lect4.html
calculate.expected.heterozygosity <- function(dlong.wprop) {
  dlong.wprop[["prop2"]] <- dlong.wprop[["prop"]]^2

  # first, sum up the squared allele proportions per locus per sample
  prop2.per.allele.per.sample <- aggregate(
    prop2 ~ sample_id + locus,
    dlong.wprop,
    sum
  )

  # second, sum up the squared allele proportions per sample, normalised by number of locus per sample
  normalised.prop2.per.sample <- aggregate(
    prop2 ~ sample_id,
    prop2.per.allele.per.sample,
    function(x) sum(x) / length(x)
  )

  # third, substract one from the calculated results to get heterozygosity
  normalised.prop2.per.sample[["expected_heterozygosity"]] <- 1 -
    normalised.prop2.per.sample[["prop2"]]

  heterozygosity.per.sample <- normalised.prop2.per.sample[, c(
    "sample_id",
    "expected_heterozygosity"
  )]

  return(heterozygosity.per.sample)
}


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/expected_heterozygosity"
dir.create(output.dir, recursive = TRUE)

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/data/mhap_filtered.tsv"

dlong <- read.delim(sfile)
dlong <- decompose.loci(dlong)

# FIXME: adjust chromosome name scheme
chrom.order <- paste0("PvP01_", sprintf("%02d", 1:14), "_v2")
dlong[["chromosome"]] <- factor(dlong[["chromosome"]], levels = chrom.order)

sorted.dlong <- sort.dlong(dlong)
dlong.wprop <- calculate.prop(sorted.dlong)
heterozygosity.per.sample <- calculate.expected.heterozygosity(dlong.wprop)

heterozygosity.per.sample.file <- paste0(
  output.dir,
  "/",
  "mhap_expected_heterozygosity.tsv"
)

write.table(
  heterozygosity.per.sample,
  file = heterozygosity.per.sample.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path to metadata
metadata.file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata.file)

# FIXME: change according sample and group column in metadata
metadata.sample.column <- "ID"
metadata.group.column <- "Year"

heterozygosity.per.sample.meta <- merge(
  metadata,
  heterozygosity.per.sample,
  by.x = metadata.sample.column,
  by.y = "sample_id"
)


z <- as.factor(metadata[[metadata.group.column]])

# FIXME: select the colour palette according to the number of groups
# The "Okabe-Ito" colour palette can accommodate up to 9 groups (default)
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
nlevels(z)

cols <- palette.colors(palette = "Tableau 10")

palette <- setNames(cols[1:nlevels(z)], levels(z))


heterozygosity.per.sample.plot.file <- paste0(
  output.dir,
  "/",
  "mhap_sample_heterozygosity.pdf"
)

# FIXME: adjust size (in inches) of expected heterozygosity plot
pdf(file = heterozygosity.per.sample.plot.file, width = 4, height = 4)

ggplot(
  heterozygosity.per.sample.meta,
  aes(
    x = as.factor(.data[[metadata.group.column]]),
    y = expected_heterozygosity
  )
) +
  geom_jitter(
    aes(colour = as.factor(.data[[metadata.group.column]])),
    width = 0.35,
    height = 0,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  xlab(metadata.group.column) +
  ylab("Expected heterozygosity") +
  theme_classic() +
  scale_colour_manual(values = palette)

dev.off()
