#!/usr/bin/env Rscript
library(ggplot2)

# FIXME: change to folder to store outputs
output.dir <- "location/to/paneljudge"
dir.create(output.dir, recursive = TRUE)

# FIXME: change according to group column in metadata
metadata.group.column <- "Year"


## setup colours here
# FIXME: select the colour palette according to the number of groups
# The "Tableau 10" colour palette can accommodate up to 10 groups (default)
# There are also "ggplot2" and other palettes from `palette.pals()`
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
cols <- palette.colors(palette = "Tableau 10")


# FIXME: change to path of diversities
diversities.file <- "location/to/paneljudge/mhap_paneljudge_diversities.tsv"

diversities <- read.delim(diversities.file)

diversities.plot.file <- paste0(
  output.dir,
  "/",
  "mhap_paneljudge_diversities.pdf"
)

diversities[[metadata.group.column]] <-
  factor(diversities[[metadata.group.column]])
diversities.palette <-
  setNames(
    cols[1:nlevels(diversities[[metadata.group.column]])],
    levels(diversities[[metadata.group.column]])
  )

# FIXME: adjust size (in inches) of paneljudge diversities plot
pdf(file = diversities.plot.file, width = 4, height = 4)

ggplot(
  diversities,
  aes(x = as.factor(.data[[metadata.group.column]]), y = diversity)
) +
  geom_jitter(
    aes(colour = as.factor(.data[[metadata.group.column]])),
    width = 0.3,
    height = 0,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  ylim(0, 1) +
  xlab(metadata.group.column) +
  ylab("Marker diversity") +
  theme_classic() +
  scale_colour_manual(
    breaks = names(diversities.palette),
    values = diversities.palette
  )

dev.off()


# FIXME: change to path of effective cardinalities
eff.cardinalities.file <- "location/to/paneljudge/mhap_paneljudge_eff_cardinalities.tsv"

eff.cardinalities <- read.delim(eff.cardinalities.file)

eff.cardinalities.plot.file <-
  paste0(output.dir, "/", "mhap_paneljudge_eff_cardinalities.pdf")

eff.cardinalities[[metadata.group.column]] <-
  factor(eff.cardinalities[[metadata.group.column]])
eff.cardinalities.palette <-
  setNames(
    cols[1:nlevels(eff.cardinalities[[metadata.group.column]])],
    levels(eff.cardinalities[[metadata.group.column]])
  )

# FIXME: adjust size (in inches) of paneljudge effective cardinalities plot
pdf(file = eff.cardinalities.plot.file, width = 4, height = 4)

ggplot(
  eff.cardinalities,
  aes(x = as.factor(.data[[metadata.group.column]]), y = eff_cardinality)
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
  ylab("Marker effective cardinality") +
  theme_classic() +
  scale_colour_manual(
    breaks = names(eff.cardinalities.palette),
    values = eff.cardinalities.palette
  )

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

rmse.krhats[[metadata.group.column]] <-
  factor(rmse.krhats[[metadata.group.column]])
rmse.krhats.palette <-
  setNames(
    cols[1:nlevels(rmse.krhats[[metadata.group.column]])],
    levels(rmse.krhats[[metadata.group.column]])
  )

# FIXME: adjust size (in inches) of paneljudge k and r estimates plot
pdf(file = krhats.plot.file, width = 6, height = 4)

ggplot(rmse.krhats, aes(x = r, y = RMSE)) +
  facet_wrap(metadata.group.column, scales = "free_y") +
  geom_line(
    aes(colour = as.factor(.data[[metadata.group.column]])),
    show.legend = FALSE
  ) +
  geom_point() +
  xlab(expression(paste("data-generating ", italic(r)))) +
  ylab(expression(paste(
    "RMSE of ",
    italic(hat(r)),
    " around data-generating ",
    italic(r)
  ))) +
  theme_classic() +
  scale_colour_manual(
    breaks = names(rmse.krhats.palette),
    values = rmse.krhats.palette
  )

dev.off()
