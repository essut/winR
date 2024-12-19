#!/usr/bin/env Rscript
library(ggplot2)

# FIXME: change according to group column in metadata
metadata.group.column <- "Year"


# FIXME: change to path of diversities
diversities.file <- "location/to/mhap_paneljudge_diversities.tsv"

diversities <- read.delim(diversities.file)

# FIXME: adjust name and size (in inches) of paneljudge diversities plot
diversities.plot.file <- "location/to/mhap_paneljudge_diversities.pdf"

pdf(file = diversities.plot.file, width = 7, height = 7)

ggplot(diversities, aes(x = as.factor(.data[[metadata.group.column]]), y = diversity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  ylim(0, 1) +
  xlab(metadata.group.column) +
  ylab("Diversity") +
  theme_classic()

dev.off()


# FIXME: change to path of effective cardinalities
eff.cardinalities.file <- "location/to/mhap_paneljudge_eff_cardinalities.tsv"

eff.cardinalities <- read.delim(eff.cardinalities.file)

# FIXME: adjust name and size (in inches) of paneljudge effective cardinalities plot
eff.cardinalities.plot.file <- "location/to/mhap_paneljudge_eff_cardinalities.pdf"

pdf(file = eff.cardinalities.plot.file, width = 7, height = 7)

ggplot(eff.cardinalities, aes(x = as.factor(.data[[metadata.group.column]]), y = eff_cardinality)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  xlab(metadata.group.column) +
  ylab("Effective Cardinality") +
  theme_classic()

dev.off()


# FIXME: change to path of estimates of k and r
krhats.file <- "location/to/mhap_paneljudge_k_r_estimates.tsv"

krhats <- read.delim(krhats.file)

rmse.krhats <-
  aggregate(
    as.formula(paste("(rhat - r)^2 ~", metadata.group.column, "+ k + r")),
    krhats,
    function(x) sqrt(mean(x))
  )
names(rmse.krhats)[length(rmse.krhats)] <- "RMSE"

# FIXME: adjust name and size (in inches) of paneljudge k and r estimates plot
krhats.plot.file <- "location/to/mhap_paneljudge_k_r_estimates.pdf"

pdf(file = krhats.plot.file, width = 7, height = 7)

ggplot(rmse.krhats, aes(x = r, y = RMSE)) +
  facet_wrap(metadata.group.column, scales = "free_y") +
  geom_point() +
  geom_line() +
  xlab(expression(paste("data-generating ", italic(r)))) +
  ylab(expression(paste("RMSE of ", italic(hat(r)), " around data-generating ", italic(r)))) +
  theme_classic()

dev.off()
