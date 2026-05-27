#!/usr/bin/env Rscript
library(readxl)
library(network)
library(ggnetwork)

## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/Dcifer"
dir.create(output.dir, recursive = TRUE)

# FIXME: change to path to metadata
metadata.file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata.file)


# FIXME: change to path of between-infection relatedness estimate
mall.estimate.file <- "location/to/Dcifer/mhap_between_relatedness_estimate.tsv"

mall.estimate <- read.delim(mall.estimate.file)

## STOP and make sure the metadata has the ID in the first column
## the ID must be the ones used in Dcifer analysis

# FIXME: change according sample and group column in metadata
metadata.sample.column <- "ID"
metadata.group.column <- "Year"

mall.estimate.meta <-
  merge(
    mall.estimate,
    metadata,
    by.x = "sample_id1",
    by.y = metadata.sample.column
  )
mall.estimate.meta <-
  merge(
    mall.estimate.meta,
    metadata,
    by.x = "sample_id2",
    by.y = metadata.sample.column
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


## setup for relatedness between-infections within groups
mall.estimate.meta.within <-
  mall.estimate.meta[
    mall.estimate.meta[[paste0(metadata.group.column, ".x")]] ==
      mall.estimate.meta[[paste0(metadata.group.column, ".y")]],
  ]

mall.estimate.meta.within[[metadata.group.column]] <-
  mall.estimate.meta.within[[paste0(metadata.group.column, ".x")]]

between.relatedness.plot.file <-
  paste0(output.dir, "/", "mhap_between_relatedness.pdf")

# FIXME: adjust size (in inches) of between-infection relatedness plot
pdf(file = between.relatedness.plot.file, width = 4, height = 4)

ggplot(
  mall.estimate.meta.within,
  aes(x = as.factor(.data[[metadata.group.column]]), y = scaled_r)
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
  ylab("relatedness between-infections (IBD)") +
  theme_classic() +
  scale_colour_manual(breaks = levels(z), values = palette)

dev.off()


prefix.filename <- paste0(output.dir, "/", "mhap_network")

# FIXME: adjust the IBD thresholds as needed
IBD.thresholds <- c(1, 1 / 2, 1 / 4, 1 / 8, 1 / 16) * 0.95

# make sure the IBD thresholds are in decreasing order
IBD.thresholds <- sort(IBD.thresholds, decreasing = TRUE)

colour.palette <- setNames(
  colorRampPalette(c("#252525", "#F0F0F0"))(length(IBD.thresholds)),
  IBD.thresholds
)
linewidth.palette <- setNames(
  seq(from = 4, to = 1, length.out = length(IBD.thresholds)),
  IBD.thresholds
)

g <-
  network(
    mall.estimate.meta[, c("sample_id1", "sample_id2", "scaled_r")],
    directed = FALSE
  )

for (i in seq_along(IBD.thresholds)) {
  IBD.threshold <- IBD.thresholds[i]

  net <- network.copy(g)
  delete.edges(net, which(get.edge.value(net, "scaled_r") < IBD.threshold))

  # make plots reproducible
  set.seed(1)
  net <- ggnetwork(net)

  net <-
    merge(
      net,
      metadata,
      by.x = "vertex.names",
      by.y = metadata.sample.column,
      sort = FALSE
    )

  # manually bin IBD values for finer control of plots
  net[["scaled_r"]] <- cut(
    net[["scaled_r"]],
    breaks = c(Inf, IBD.thresholds),
    labels = rev(IBD.thresholds),
    right = FALSE
  )

  # sort by relatedness to emphasize close relationships
  net <- net[order(net[["scaled_r"]]), ]

  # FIXME: adjust size (in inches) of network plot
  pdf(
    paste0(prefix.filename, "_", IBD.threshold, ".pdf"),
    width = 9,
    height = 8
  )

  # FIXME: relatedness network framework, adjust metadata as necessary
  print(
    ggplot(
      net,
      aes(
        x,
        y,
        xend = xend,
        yend = yend,
        fill = as.factor(.data[[metadata.group.column]])
      )
    ) +
      geom_edges(aes(colour = scaled_r, linewidth = scaled_r)) +
      geom_nodes(size = 5, shape = 21) +
      theme_blank() +
      labs(fill = metadata.group.column) +
      # ensure legend order is consistent for all plots
      guides(fill = guide_legend(order = 1)) +
      scale_fill_manual(breaks = levels(z), values = palette) +
      scale_colour_manual(
        name = "IBD",
        values = colour.palette,
        breaks = IBD.thresholds[1:i],
        labels = paste0(">=", IBD.thresholds[1:i] * 100, "%")
      ) +
      scale_linewidth_manual(
        name = "IBD",
        values = linewidth.palette,
        breaks = IBD.thresholds[1:i],
        labels = paste0(">=", IBD.thresholds[1:i] * 100, "%")
      )
  )

  dev.off()
}
