#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-101
outputCIGAR.to.long <- function(outputCIGAR, keep.unused.alleles = FALSE) {
  long <-
    reshape(
      outputCIGAR,
      direction = "long",
      varying = names(outputCIGAR)[-1],
      v.names = "count",
      timevar = "pseudoCIGAR",
      times = names(outputCIGAR)[-1],
      idvar = "sample_id",
      ids = outputCIGAR[["sample"]]
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


rmindel.allele <- function(long) {
  # insertion pattern: {position}I={bases}
  # deletion pattern: {position}D={bases}
  long[["allele"]] <- gsub("([0-9]+)[DI]=([ACGT]+)", "", long[["allele"]])
  
  # default to wild type if all variants were removed
  long[long[["allele"]] %in% "", "allele"] <- "."
  
  # consolidate allele counts after indel removal
  aggregate(count ~ sample_id + locus + allele, long, sum)
}


filter.long <- function(long, minimum.count, keep.marker = "microhaplotype") {
  long <-
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
  
  long.counts.per.locus <- aggregate(count ~ sample_id + locus, long, sum)
  
  pass.filter <-
    long.counts.per.locus[long.counts.per.locus[["count"]] >= minimum.count, ]
  
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

# FIXME: change to outputCIGAR.tsv file path
outputCIGAR.file <- "location/to/outputCIGAR.tsv"

outputCIGAR <- read.delim(outputCIGAR.file, check.names = FALSE)
long <- outputCIGAR.to.long(outputCIGAR)

# FIXME: comment the line below if you want to keep insertions and deletions in alleles
long <- rmindel.allele(long)

# FIXME: change to path of unfiltered data
long.file <- "location/to/long_unfiltered.tsv"

write.table(long, long.file, quote = FALSE, sep = "\t", row.names = FALSE)


# FIXME: change read-pair threshold for failed genotype
minimum.count <- 10

mhap.filtered <- filter.long(long, minimum.count, keep.marker = "microhaplotype")
mhap.filtered.nloci.per.sample <- calculate.remaining.nloci(mhap.filtered)


# FIXME: change to path for statistics on remaining locus (text)
mhap.filtered.nloci.per.sample.text <- "location/to/mhap_filtered_nloci.tsv"

write.table(
  mhap.filtered.nloci.per.sample,
  mhap.filtered.nloci.per.sample.text,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change remaining locus threshold for failed samples
minimum.nloci <- 75 # ~80% markers

# FIXME: change to path for statistics on remaining locus (figure)
mhap.filtered.nloci.per.sample.figure <- "location/to/mhap_filtered_nloci.pdf"

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


# FIXME: change to path of filtered microhaplotype data
mhap.filtered.subset.file <- "location/to/mhap_filtered.tsv"

write.table(
  mhap.filtered.subset,
  mhap.filtered.subset.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
