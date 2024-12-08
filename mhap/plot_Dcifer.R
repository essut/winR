#!/usr/bin/env Rscript
library(readxl)
library(network)
library(sna)
library(GGally)

## Change FIXME lines to the appropriate parameters

# FIXME: change to path to metadata
metadata.file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata.file)


# FIXME: change to path of between-infection relatedness estimate
mall.estimate.file <- "location/to/mhap_between_relatedness_estimate.tsv"

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


mall.estimate.meta <-
  mall.estimate.meta[
    mall.estimate.meta[[paste0(metadata.group.column, ".x")]] ==
      mall.estimate.meta[[paste0(metadata.group.column, ".y")]],
    
  ]


# FIXME: adjust name and size (in inches) of between-infection relatedness plot
between.relatedness.plot.file <- "location/to/mhap_between_relatedness.pdf"

pdf(file = between.relatedness.plot.file, width = 7, height = 7)

ggplot(
  mall.estimate.meta,
  aes(x = as.factor(.data[[paste0(metadata.group.column, ".x")]]), y = relatedness)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  ylim(0, 1) +
  xlab(metadata.group.column) +
  theme_classic()

dev.off()


# FIXME: adjust the IBD thresholds as needed
IBD.thresholds <- c(1, 1/2, 1/4, 1/8) * 0.95

# FIXME: determine path to save plots
prefix.filename <- "location/to/mhap_network"

g <- network(mall.estimate, directed = FALSE, vertices = metadata)

for (IBD.threshold in IBD.thresholds) {
  net <- network.copy(g)
  net <-
    delete.edges(net, which(get.edge.value(net, "relatedness") < IBD.threshold))

  # FIXME: adjust size (in inches) of network plot
  pdf(paste0(prefix.filename, "_", IBD.threshold, ".pdf"), width = 7, height = 7.5)
  
  # FIXME: adjust metadata as necessary
  print(
    ggnet2(
      net,
      color = metadata.group.column,
      size = 6,
      palette = "Paired",
      legend.position = "bottom"
    ) +
      geom_point(aes(color = color), size = 6, color = "black") +
      geom_point(aes(color = color), size = 5) +
      ggtitle(IBD.threshold)
  )
  
  dev.off()
}
