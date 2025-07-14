#!/usr/bin/env Rscript
library(ggplot2)

# FIXME: change to folder to store outputs
output.dir <- "location/to/paneljudge"
dir.create(output.dir, recursive = TRUE)

# FIXME: change according to group column in metadata
metadata.group.column <- "Year"


# FIXME: change to path of diversities
diversities.file <- "location/to/paneljudge/mhap_paneljudge_diversities.tsv"

diversities <- read.delim(diversities.file)

diversities.plot.file <- paste0(output.dir, "/", "mhap_paneljudge_diversities.pdf")

# FIXME: adjust size (in inches) of paneljudge diversities plot
pdf(file = diversities.plot.file, width = 4, height = 4)

ggplot(diversities, aes(x = as.factor(.data[[metadata.group.column]]), y = diversity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.5) +
  ylim(0, 1) +
  xlab(metadata.group.column) +
  ylab("Diversity") +
  theme_classic()

dev.off()


# FIXME: change to path of effective cardinalities
eff.cardinalities.file <- "location/to/paneljudge/mhap_paneljudge_eff_cardinalities.tsv"

eff.cardinalities <- read.delim(eff.cardinalities.file)

eff.cardinalities.plot.file <-
  paste0(output.dir, "/", "mhap_paneljudge_eff_cardinalities.pdf")

# FIXME: adjust size (in inches) of paneljudge effective cardinalities plot
pdf(file = eff.cardinalities.plot.file, width = 4, height = 4)

ggplot(eff.cardinalities, aes(x = as.factor(.data[[metadata.group.column]]), y = eff_cardinality)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.5) +
  xlab(metadata.group.column) +
  ylab("Effective Cardinality") +
  theme_classic()

dev.off()


# FIXME: change to path of estimates of k and r
krhats.file <- "location/to/paneljudge/mhap_paneljudge_k_r_estimates.tsv"

krhats <- read.delim(krhats.file)

rmse.krhats <-
  aggregate(
    as.formula(paste("(rhat - r)^2 ~", metadata.group.column, "+ k + r")),
    krhats,
    function(x) sqrt(mean(x))
  )
names(rmse.krhats)[length(rmse.krhats)] <- "RMSE"

krhats.plot.file <- paste0(output.dir, "/", "mhap_paneljudge_k_r_estimates.pdf")

# FIXME: adjust size (in inches) of paneljudge k and r estimates plot
pdf(file = krhats.plot.file, width = 6, height = 4)

ggplot(rmse.krhats, aes(x = r, y = RMSE)) +
  facet_wrap(metadata.group.column, scales = "free_y") +
  geom_point() +
  geom_line() +
  xlab(expression(paste("data-generating ", italic(r)))) +
  ylab(expression(paste("RMSE of ", italic(hat(r)), " around data-generating ", italic(r)))) +
  theme_classic()

dev.off()
