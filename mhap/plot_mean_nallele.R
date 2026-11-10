#!/usr/bin/env Rscript
library(readxl)
library(ggplot2)

## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/mean_nallele"
dir.create(output.dir, recursive = TRUE)

# FIXME: change to path of filtered microhaplotype data
mhap.filtered <- read.delim("location/to/data/mhap_filtered.tsv")

mhap.filtered.sum.unique.alleles.per.marker <- aggregate(
  allele ~ sample_id + locus,
  mhap.filtered,
  function(x) length(unique(x))
)

mhap.filtered.mean.nallele <- aggregate(
  allele ~ sample_id,
  mhap.filtered.sum.unique.alleles.per.marker,
  mean
)
names(mhap.filtered.mean.nallele)[2] <- "mean_nallele"

mhap.filtered.mean.nallele.file <- paste0(
  output.dir,
  "/",
  "mhap_mean_nallele.tsv"
)

write.table(
  mhap.filtered.mean.nallele,
  file = mhap.filtered.mean.nallele.file,
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

mhap.filtered.mean.nallele.meta <- merge(
  mhap.filtered.mean.nallele,
  metadata,
  by.x = "sample_id",
  by.y = metadata.sample.column
)

mhap.filtered.mean.nallele.meta.summary <- aggregate(
  as.formula(paste("mean_nallele ~", metadata.group.column)),
  mhap.filtered.mean.nallele.meta,
  summary
)
mhap.filtered.mean.nallele.meta.summary <- cbind(
  mhap.filtered.mean.nallele.meta.summary[[metadata.group.column]],
  as.data.frame(mhap.filtered.mean.nallele.meta.summary[["mean_nallele"]])
)
names(mhap.filtered.mean.nallele.meta.summary)[1] <- metadata.group.column

mhap.filtered.mean.nallele.meta.summary.file <- paste0(
  output.dir,
  "/",
  "mhap_mean_nallele_summary.csv"
)

write.csv(
  mhap.filtered.mean.nallele.meta.summary,
  file = mhap.filtered.mean.nallele.meta.summary.file,
  row.names = FALSE
)


## setup colours here
z <- as.factor(metadata[[metadata.group.column]])

# FIXME: select the colour palette according to the number of groups
# The "Tableau 10" colour palette can accommodate up to 10 groups (default)
# There are also "ggplot2" and other palettes from `palette.pals()`
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
nlevels(z)

cols <- palette.colors(palette = "Tableau 10")

palette <- setNames(cols[1:nlevels(z)], levels(z))


mhap.filtered.mean.nallele.meta.plot.file <- paste0(
  output.dir,
  "/",
  "mhap_mean_nallele_summary.pdf"
)

# FIXME: adjust size (in inches) of mean number of alleles plot
pdf(file = mhap.filtered.mean.nallele.meta.plot.file, width = 4, height = 4)

ggplot(
  mhap.filtered.mean.nallele.meta,
  aes(x = as.factor(.data[[metadata.group.column]]), y = mean_nallele)
) +
  geom_jitter(
    aes(colour = as.factor(.data[[metadata.group.column]])),
    width = 0.3,
    height = 0,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  xlab(metadata.group.column) +
  ylab("Mean number of alleles") +
  theme_classic() +
  scale_colour_manual(breaks = levels(z), values = palette)

dev.off()


# FIXME: determine cutoff of high complexity infection based on mean number of alleles
high.complexity.mean.nallele.cutoff <- 1.05

mhap.filtered.mean.nallele.sorted.file <- paste0(
  output.dir,
  "/",
  "mhap_mean_nallele_sorted.pdf"
)

# FIXME: adjust size (in inches) of mean number of alleles plot
pdf(file = mhap.filtered.mean.nallele.sorted.file, width = 4, height = 4)
par(mar = c(1, 4, 1, 1) + 0.1)
plot(
  sort(mhap.filtered.mean.nallele[["mean_nallele"]]),
  ylab = "Mean number of alleles",
  xaxt = "n",
  xlab = NA
)
abline(h = high.complexity.mean.nallele.cutoff, col = "red", lt = "dashed")
dev.off()

mhap.filtered.mean.nallele[["complexity"]] <- ifelse(
  mhap.filtered.mean.nallele[["mean_nallele"]] >
    high.complexity.mean.nallele.cutoff,
  "high",
  "low"
)

mhap.filtered.mean.nallele.complexity.file <- paste0(
  output.dir,
  "/",
  "mhap_complexity_status.tsv"
)

write.table(
  mhap.filtered.mean.nallele[, c("sample_id", "complexity")],
  file = mhap.filtered.mean.nallele.complexity.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
