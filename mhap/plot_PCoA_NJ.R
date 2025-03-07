#!/usr/bin/env Rscript
library(readxl)
library(ape)
library(ggplot2)

## Load all functions by running everything from line 7-151
create.allele.matrix <- function(dlong) {
  # calls the major allele in each locus for each sample
  dlong.major <-
    merge(aggregate(count ~ sample_id + locus, data = dlong, max), dlong)
  
  # convert as a matrix, randomly picks one allele if multiple exists
  dlong.major.wide <-
    reshape(
      dlong.major[, c("sample_id", "locus", "allele")],
      direction = "wide",
      idvar = "sample_id",
      timevar = "locus"
    )
  
  row.names(dlong.major.wide) <- dlong.major.wide[["sample_id"]]
  dlong.major.wide[["sample_id"]] <- NULL
  
  return(dlong.major.wide)
}


calculate.missingness <- function(dlong.major.wide) {
  missing.alleles <- is.na(dlong.major.wide)
  missing.locus.per.sample <- rowSums(missing.alleles)
  missing.sample.per.locus <- colSums(missing.alleles)
  
  return(
    list(
      missing.locus.per.sample = missing.locus.per.sample,
      missing.sample.per.locus = missing.sample.per.locus
    )
  )
}


get.usable.markers <- function(dlong.major.wide.missingness) {
  return(sum(dlong.major.wide.missingness[["missing.sample.per.locus"]] == 0))
}


plot.missing.sample.per.locus <-
  function(dlong.major.wide.missingness, missing.sample.per.locus.cutoff) {
    hist(
      dlong.major.wide.missingness[["missing.sample.per.locus"]],
      breaks = "FD",
      right = FALSE,
      main = "Number of missing sample per locus",
      xlab = NULL
    )
    
    abline(v = missing.sample.per.locus.cutoff + 1, col = "red", lty = "dashed")
  }


plot.missing.locus.per.sample <-
  function(dlong.major.wide.missingness, missing.locus.per.sample.cutoff) {
    hist(
      dlong.major.wide.missingness[["missing.locus.per.sample"]],
      breaks = "FD",
      right = FALSE,
      main = "Number of missing locus per sample",
      xlab = NULL
    )
    
    abline(v = missing.locus.per.sample.cutoff + 1, col = "red", lty = "dashed")
  }


.compute.parents.groups <-
  function(edge.groups, edges, mixed.group) {
    # determine what the parents' group would be based on
    # the groups passed on from the children
    tapply(edge.groups[, 1], edges[, 1], function(x) {
      # retrieve children's groups
      unique.groups <- na.omit(unique(x))
      
      # if the children's groups are mixed
      if (mixed.group %in% unique.groups) {
        return(mixed.group)
      }
      
      unique.length <- length(unique.groups)
      
      # if no information about the children's groups
      # we'll return back to this parent later
      if (unique.length == 0) {
        return(NA)
      }
      
      # if the children only have one group
      if (unique.length == 1) {
        return(unique.groups)
      }
      
      # if the children have multiple groups
      return(mixed.group)
    })
  }

.update.parents.groups <-
  function(edge.groups, edges, mixed.group) {
    parents.groups <-
      .compute.parents.groups(edge.groups, edges, mixed.group)
    
    parent.as.child <- match(edges[, 2], names(parents.groups))
    
    # if the parents are children of other parents,
    # carry the group of parents as the child's group
    edge.groups[, 2] <-
      mapply(
        function(x, y) {
          if (is.na(y)) {
            return(x)
          } else {
            return(y)
          }
        },
        x = edge.groups[, 2],
        y = parents.groups[parent.as.child]
      )
    
    edge.groups[, 1] <- edge.groups[, 2]
    
    edge.groups
  }

group.NJ <- function(NJ, tip.groups, mixed.group) {
  edges <- NJ[["edge"]]
  
  # copy available tip groups as edge groups
  edge.groups <- tip.groups[edges]
  
  # copy groups on the tip to its direct parent
  dim(edge.groups) <- dim(edges)
  edge.groups[, 1] <- edge.groups[, 2]
  
  # trace back parent groups until they meet different groups
  while (sum(is.na(edge.groups))) {
    edge.groups <-
      .update.parents.groups(edge.groups, edges, mixed.group)
  }
  
  # finalise edge groups
  .update.parents.groups(edge.groups, edges, mixed.group)
}


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to path to metadata
metadata.file <- "location/to/metadata.xlsx"

# FIXME: adjust based on file format
metadata <- read_excel(metadata.file)


# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/mhap_filtered.tsv"

dlong <- read.delim(sfile)

# FIXME: change to path of polyclonal status
polyclonal.status.file <- "location/to/mhap_polyclonal_status.tsv"

polyclonal.status <- read.delim(polyclonal.status.file)

## filter to monoclonal samples
dlong <-
  dlong[
    dlong[["sample_id"]] %in%
      polyclonal.status[!polyclonal.status[["is_polyclonal"]], "sample_id"],
    
    ]

dlong.major.wide <- create.allele.matrix(dlong)
dlong.major.wide.missingness <- calculate.missingness(dlong.major.wide)

## check number of samples and markers currently available
dim(dlong.major.wide)

## check how many markers will be used in downstream analyses
get.usable.markers(dlong.major.wide.missingness)


# FIXME: set a cutoff to remove as many missing locus as possible
missing.sample.per.locus.cutoff <- 5

# FIXME: change to path for missing sample per locus
missing.sample.per.locus.plot.file <- "location/to/mhap_missing_sample_per_locus.pdf"

pdf(file = missing.sample.per.locus.plot.file)
plot.missing.sample.per.locus(
  dlong.major.wide.missingness,
  missing.sample.per.locus.cutoff
)
dev.off()


## check if happy with the cutoff for missing sample per locus
## if not, go back and change the cutoff
dlong.major.wide.locfilt <-
  dlong.major.wide[
    ,
    dlong.major.wide.missingness[["missing.sample.per.locus"]] <=
      missing.sample.per.locus.cutoff
  ]


dlong.major.wide.locfilt.missingness <-
  calculate.missingness(dlong.major.wide.locfilt)

## check number of samples and markers currently available
dim(dlong.major.wide.locfilt)

## check how many markers will be used in downstream analyses
get.usable.markers(dlong.major.wide.locfilt.missingness)


# FIXME: set a cutoff to remove as many missing sample as possible
missing.locus.per.sample.cutoff <- 0

# FIXME: change to path for missing locus per sample
missing.locus.per.sample.plot.file <- "location/to/mhap_missing_locus_per_sample.pdf"

pdf(file = missing.locus.per.sample.plot.file)
plot.missing.locus.per.sample(
  dlong.major.wide.locfilt.missingness,
  missing.locus.per.sample.cutoff
)
dev.off()


## check if happy with the cutoff for missing sample per locus
## if not, go back and change the cutoff
dlong.major.wide.locfilt.samfilt <-
  dlong.major.wide.locfilt[
    dlong.major.wide.locfilt.missingness[["missing.locus.per.sample"]] <=
      missing.locus.per.sample.cutoff
    ,
  ]


dlong.major.wide.locfilt.samfilt.missingness <-
  calculate.missingness(dlong.major.wide.locfilt.samfilt)

## check number of samples and markers currently available
dim(dlong.major.wide.locfilt.samfilt)

## check how many markers will be used in downstream analyses
get.usable.markers(dlong.major.wide.locfilt.samfilt.missingness)


dist <- dist.gene(dlong.major.wide.locfilt.samfilt)

# FIXME: change according sample and group column in metadata
metadata.sample.column <- "ID"
metadata.group.column <- "Year"

z <- as.factor(metadata[[metadata.group.column]])

# FIXME: select the colour palette according to the number of groups
# The "Okabe-Ito" colour palette can accommodate up to 10 groups (default)
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
nlevels(z)
cols <- palette.colors(palette = "Okabe-Ito")

palette <- setNames(cols[1:nlevels(z)], levels(z))


## lines below create PCoA plots

res <- pcoa(dist)
vectors <- res[["vectors"]]
Broken_stick <- res[["values"]][["Broken_stick"]] * 100


data <-
  merge(
    vectors,
    metadata,
    by.x = "row.names",
    by.y = metadata.sample.column
  )


# FIXME: determine path to save PCoA plots
PCoA.prefix.filename <- "location/to/mhap_PCoA"

# FIXME: number of PCoA axes to plot
n.axis <- 3

m <- combn(n.axis, 2)
for (j in seq_len(ncol(m))) {
  x <- m[1, j]
  y <- m[2, j]
  
  # FIXME: adjust size (in inches) of PCoA plot
  pdf(paste0(PCoA.prefix.filename, "_", x, "_", y, ".pdf"), width = 5, height = 5)
  
  # FIXME: adjust PCoA framework as required
  print(
    ggplot(
      data,
      aes(
        x = .data[[paste0("Axis.", x)]],
        y = .data[[paste0("Axis.", y)]],
        colour = as.factor(.data[[metadata.group.column]])
      )
    ) +
      geom_point(alpha = 0.7, size = 4) +
      theme_classic() +
      labs(colour = metadata.group.column) +
      theme(legend.position = "bottom") +
      xlab(paste0("Coordinate ", x, " (", round(Broken_stick[x], digits = 2), "%)")) +
      ylab(paste0("Coordinate ", y, " (", round(Broken_stick[y], digits = 2), "%)")) +
      scale_colour_manual(values = palette) +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))
  )
  
  dev.off()
}


## lines below create NJ plots
tr <- nj(dist)
cls <- setNames(metadata[[metadata.group.column]], metadata[[metadata.sample.column]])
tip.groups <- cls[tr[["tip.label"]]]

tip.colours <- palette[as.character(tip.groups)]
edge.colours <- group.NJ(tr, tip.colours, "#E5E4E2")
## if the script hangs at edge.colours,
## there is probably a sample with no available metadata


# FIXME: determine path to save NJ plots
NJ.prefix.filename <- "location/to/mhap_NJ"


## rooted, labelled
# FIXME: adjust size (in inches) of NJ plot
pdf(paste0(NJ.prefix.filename, "_rooted_labelled.pdf"), width = 7, height = 7.5)
par(xpd = TRUE)

plot(
  tr,
  type = "fan",
  underscore = TRUE,
  lab4ut = "axial",
  show.tip.label = TRUE,
  edge.color = edge.colours,
  tip.color = tip.colours
)

legend(
  "bottom",
  levels(z),
  fill = palette,
  ncol = 4,
  title = metadata.group.column,
  inset = -0.17
)

dev.off()

## unrooted, unlabelled
# FIXME: adjust size (in inches) of NJ plot
pdf(paste0(NJ.prefix.filename, "_unrooted_unlabelled.pdf"), width = 7, height = 7.5)
par(xpd = TRUE)

plot(
  tr,
  type = "unrooted",
  underscore = TRUE,
  lab4ut = "axial",
  show.tip.label = FALSE,
  edge.color = edge.colours,
  tip.color = tip.colours
)

legend(
  "bottom",
  levels(z),
  fill = palette,
  ncol = 4,
  title = metadata.group.column,
  inset = -0.17
)

dev.off()
