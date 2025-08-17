#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-111
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
  sorted.dlong <-
    dlong[do.call(order, dlong[, c("chromosome", "start.pos", "end.pos", "allele")]), ]
  
  locus.order <- unique(sorted.dlong[["locus"]])
  sorted.dlong[["locus"]] <- factor(sorted.dlong[["locus"]], levels = locus.order)
  
  marker.order <- unique(sorted.dlong[["marker"]])
  sorted.dlong[["marker"]] <- factor(sorted.dlong[["marker"]], levels = marker.order)

  sorted.dlong
}


calculate.prop <- function(sorted.dlong) {
  dlong.count <- aggregate(count ~ sample_id + locus, sorted.dlong, sum)
  names(dlong.count)[length(dlong.count)] <- "total"
  
  sorted.dlong.wprop <- merge(sorted.dlong, dlong.count, sort = FALSE)
  sorted.dlong.wprop[["prop"]] <-
    sorted.dlong.wprop[["count"]] / sorted.dlong.wprop[["total"]]
  
  sorted.dlong.wprop
}


## assumes loci are ordered
plot.rainbow <- function(dlong.wprop, sample_id, add.marker.label = TRUE, sort.by.prop = FALSE) {
  if (length(sample_id) != 1) {
    stop("Please select a single sample to plot.")
  }
  
  index <- dlong.wprop[["sample_id"]] %in% sample_id
  if (any(index) == FALSE) {
    stop("Sample not found in the dataset: ", sample_id)
  }
  
  sample.dlong <- dlong.wprop[index, ]

  # maximum number of alleles per locus is dependent on the dataset
  dataset.alleles <- unique(dlong.wprop[, c("marker", "allele")])
  
  sample.dlong <- merge(dataset.alleles, sample.dlong, all = TRUE)
  
  allele.prop <-
    tapply(sample.dlong[["prop"]], sample.dlong[["marker"]], unlist)
  allele.prop <-
    lapply(allele.prop, function(x) {
      x[is.na(x)] <- 0
      return(x)
    })
  
  allele.length <- lengths(allele.prop)
  allele.pad <-
    lapply(max(allele.length) - allele.length, function(n) rep(NA, n))
  
  allele.complete <-
    mapply(function(x, y) c(x, y), allele.prop, allele.pad)
  
  # plot each bar according to number of alleles per marker
  first.plot <- TRUE
  for (marker in colnames(allele.complete)) {
    plot.device <- allele.complete
    
    # empty out all non-plot
    plot.device[, !colnames(allele.complete) %in% marker] <- NA
    
    col <- scales::pal_hue()(allele.length[marker])
    
    if (sort.by.prop) {
      unsorted <- plot.device[, colnames(allele.complete) %in% marker]
      names(unsorted) <- col
      
      sorted <- sort(unsorted, decreasing = TRUE, na.last = TRUE)
      
      col <- names(sorted)[seq_len(allele.length[marker])]
      plot.device[, colnames(allele.complete) %in% marker] <- sorted
    }
    
    if (first.plot) {
      barplot(plot.device, col = col, axisnames = add.marker.label, las = 2)
      first.plot <- FALSE
    } else {
      barplot(plot.device, col = col, axes = FALSE, axisnames = FALSE, add = TRUE)
    }
  }
}


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/rainbow"
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

# FIXME: by default, plot all samples
sample_ids <- sort(unique(dlong.wprop[["sample_id"]]))

for (sample_id in sample_ids) {
  rainbow.plot.file <- paste0(output.dir, "/", sample_id, "_rainbow_plot.pdf")
  
  # FIXME: adjust size of "rainbow" plot
  pdf(file = rainbow.plot.file, width = 20, height = 5)
  
  plot.rainbow(dlong.wprop, sample_id)
  title(sample_id)
  
  dev.off()
}
