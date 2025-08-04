#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-146
output.to.long <- function(output, keep.unused.alleles = FALSE) {
  long <-
    reshape(
      output,
      direction = "long",
      varying = names(output)[-1],
      v.names = "count",
      timevar = "pseudoCIGAR",
      times = names(output)[-1],
      idvar = "sample_id",
      ids = output[["sample"]]
    )
  row.names(long) <- NULL
  
  if (!keep.unused.alleles) {
    long <- long[long[["count"]] > 0, ]
  }
  
  # format data to required columns
  locus.allele <- strsplit(long[["pseudoCIGAR"]], ",")
  long[["locus"]] <- vapply(locus.allele, "[[", character(1), 1)
  long[["allele"]] <- vapply(locus.allele, "[[", character(1), 2)
  
  long[, c("sample_id", "locus", "allele", "count")]
}


merge.outputs <- function(output.files) {
  longs <- list()
  
  for (i in seq_along(output.files)) {
    output.file <- output.files[i]
    output <- read.delim(output.file, check.names = FALSE)
    longs[[i]] <- output.to.long(output)
  }
  
  long <- do.call(rbind, longs)
  
  # consolidate allele counts from different runs
  aggregate(count ~ sample_id + locus + allele, long, sum)
}


.rmindel.allele.outputCIGAR <- function(long) {
  # insertion pattern: {position}I={bases}
  # deletion pattern: {position}D={bases}
  long[["allele"]] <- gsub("([0-9]+)[DI]=([ACGT]+)", "", long[["allele"]])
  
  # default to wild type if all variants were removed
  long[long[["allele"]] %in% "", "allele"] <- "."
  
  long
}


rmindel.allele <- function(long) {
  long <- .rmindel.allele.outputCIGAR(long)
  
  # consolidate allele counts after indel removal
  aggregate(count ~ sample_id + locus + allele, long, sum)
}


select.markers <- function(long, keep.marker = "microhaplotype") {
  switch (
    keep.marker,
    microhaplotype = {
      long[!grepl("MIT|DHPS|MDR1", long[["locus"]]), ]
    },
    drugR = {
      long[grepl("DHPS|MDR1", long[["locus"]]), ]
    },
    mitochondria = {
      long[grepl("MIT", long[["locus"]]), ]
    },
    stop("Valid options are 'microhaplotype', 'drugR', 'mitochondria'")
  )
}


.calculate.prop <- function(long) {
  dlong.count <- aggregate(count ~ sample_id + locus, long, sum)
  names(dlong.count)[length(dlong.count)] <- "total"
  
  long.wprop <- merge(long, dlong.count, sort = FALSE)
  long.wprop[["prop"]] <-
    long.wprop[["count"]] / long.wprop[["total"]]
  
  long.wprop
}


allele.proportion.filter <- function(long, allele.proportion.cutoff) {
  long.wprop <- .calculate.prop(long)
  
  long.wprop[
    long.wprop[["prop"]] >= allele.proportion.cutoff,
    c("sample_id", "locus", "allele", "count")
  ]
}


allele.count.filter <- function(long, allele.count.cutoff) {
  long[long[["count"]] >= allele.count.cutoff, ]
}


minimum.total.filter <- function(long, minimum.total.cutoff) {
  long.counts.per.locus <- aggregate(count ~ sample_id + locus, long, sum)
  
  pass.filter <-
    long.counts.per.locus[long.counts.per.locus[["count"]] >= minimum.total.cutoff, ]
  
  merge(long, pass.filter[, c("sample_id", "locus")], sort = FALSE)
}


calculate.remaining.nloci <- function(long.filtered) {
  long.filtered.nloci.per.sample <-
    aggregate(locus ~ sample_id, long.filtered, function(x) length(unique(x)))
  
  names(long.filtered.nloci.per.sample)[2] <- "nloci"
  
  long.filtered.nloci.per.sample
}


plot.remaining.nloci <- function(long.filtered.nloci.per.sample, minimum.nloci) {
  plot(
    sort(long.filtered.nloci.per.sample[["nloci"]]),
    main = NA,
    xlab = NA,
    ylab = "Remaining number of locus",
    xaxt = "n"
  )
  
  abline(h = minimum.nloci, col = "red", lty = "dashed")
}


filter.failed.samples <-
  function(long.filtered, long.filtered.nloci.per.sample, minimum.nloci) {
    pass.samples <-
      long.filtered.nloci.per.sample[
        long.filtered.nloci.per.sample[["nloci"]] >= minimum.nloci,
        "sample_id"
      ]
    
    long.filtered[long.filtered[["sample_id"]] %in% pass.samples, ]
  }


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/data"
dir.create(output.dir, recursive = TRUE)

# FIXME: change to outputCIGAR.tsv / outputHaplotypes.tsv file paths
output.files <-
  c(
    "location/to/outputCIGAR.tsv",
    "another/outputCIGAR.tsv"
  )

# use outputs for complete sample list
unfiltered.samples <- character()
for (output.file in output.files) {
  unfiltered.samples <-
    c(unfiltered.samples, read.delim(output.file, check.names = FALSE)[["sample"]])
}
unfiltered.samples <- sort(unique(unfiltered.samples))
unfiltered.samples <- data.frame(sample_id = unfiltered.samples)

print(unfiltered.samples)
## STOP and check the sample list, correct the outputs as necessary

long <- merge.outputs(output.files)


# FIXME: remove samples if needed (e.g. samples from a different cohort)
removed.samples <- c("sWGA_1_10", "sWGA_1_20", "sWGA_1_30")

unfiltered.samples <-
  unfiltered.samples[
    !unfiltered.samples[["sample_id"]] %in% removed.samples,
    ,
    drop = FALSE
  ]

long <- long[!long[["sample_id"]] %in% removed.samples, ]

  
# FIXME: comment the line below if you want to keep insertions and deletions in alleles
long <- rmindel.allele(long)

long.file <- paste0(output.dir, "/", "long_unfiltered.tsv")

write.table(long, long.file, quote = FALSE, sep = "\t", row.names = FALSE)


# FIXME: change read-pair threshold for failed genotype
allele.proportion.cutoff <- 0.02 # allele proportion is at least 2% per locus per sample
allele.count.cutoff <- 5 # at least 5 read-pairs needed for each allele
# minimum.total.cutoff <- 10 # at least 10 read-pairs in total per locus per sample

# FIXME: comment any lines below to disable certain filters
mhap <- select.markers(long, keep.marker = "microhaplotype")
mhap.filtered <- allele.proportion.filter(mhap, allele.proportion.cutoff)
mhap.filtered <- allele.count.filter(mhap.filtered, allele.count.cutoff)
# mhap.filtered <- minimum.total.filter(mhap.filtered, minimum.total.cutoff)
mhap.filtered.nloci.per.sample <- calculate.remaining.nloci(mhap.filtered)

mhap.filtered.nloci.per.sample <-
  merge(mhap.filtered.nloci.per.sample, unfiltered.samples, all.y = TRUE)
mhap.filtered.nloci.per.sample[
  is.na(mhap.filtered.nloci.per.sample[["nloci"]]),
  "nloci"
] <- 0

mhap.filtered.nloci.per.sample.text <-
  paste0(output.dir, "/", "mhap_filtered_nloci.tsv")

write.table(
  mhap.filtered.nloci.per.sample,
  mhap.filtered.nloci.per.sample.text,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change remaining locus threshold for failed samples
minimum.nloci <- 75 # ~80% markers

mhap.filtered.nloci.per.sample.figure <-
  paste0(output.dir, "/", "mhap_filtered_nloci.pdf")

pdf(file = mhap.filtered.nloci.per.sample.figure)
plot.remaining.nloci(mhap.filtered.nloci.per.sample, minimum.nloci)
dev.off()

mhap.filtered.subset <-
  filter.failed.samples(
    mhap.filtered,
    mhap.filtered.nloci.per.sample,
    minimum.nloci
  )

## STOP and look at the generated text and figure statistics
## Was the remaining locus threshold appropriate for this dataset?

print(data.frame(sample_id = sort(unique(mhap.filtered.subset[["sample_id"]]))))

# FIXME: remove any unused samples (e.g. controls)
unused.samples <- c("CQE98-Pv")

mhap.filtered.subset <-
  mhap.filtered.subset[!mhap.filtered.subset[["sample_id"]] %in% unused.samples, ]

mhap.filtered.subset.file <- paste0(output.dir, "/", "mhap_filtered.tsv")

write.table(
  mhap.filtered.subset,
  mhap.filtered.subset.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
