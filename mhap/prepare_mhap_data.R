#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-122
outputCIGAR.to.long <- function(outputCIGAR, keep.unused.alleles = FALSE) {
  outputCIGAR.long <-
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
  row.names(outputCIGAR.long) <- NULL
  
  if (!keep.unused.alleles) {
    outputCIGAR.long <- outputCIGAR.long[outputCIGAR.long[["count"]] > 0, ]
  }
  
  # format data to required columns
  locus.allele <- strsplit(outputCIGAR.long[["pseudoCIGAR"]], ",")
  outputCIGAR.long[["locus"]] <-
    vapply(locus.allele, "[[", character(1), 1)
  outputCIGAR.long[["allele"]] <-
    vapply(locus.allele, "[[", character(1), 2)
  outputCIGAR.long <-
    outputCIGAR.long[, c("sample_id", "locus", "allele", "count")]
  
  outputCIGAR.long
}


filter.outputCIGAR.long <-
  function(
    outputCIGAR.long,
    minimum.count,
    keep.marker = "microhaplotype"
  ) {
    
    outputCIGAR.long <-
      switch (
        keep.marker,
        microhaplotype = {
          outputCIGAR.long[!grepl("MIT|DHPS|MDR1", outputCIGAR.long[["locus"]]), ]
        },
        drugR = {
          outputCIGAR.long[grepl("DHPS|MDR1", outputCIGAR.long[["locus"]]), ]
        },
        mitochondria = {
          outputCIGAR.long[grepl("MIT", outputCIGAR.long[["locus"]]), ]
        },
        stop("Valid options are 'microhaplotype', 'drugR', 'mitochondria'")
      )
    
    outputCIGAR.long.counts.per.locus <-
      aggregate(count ~ sample_id + locus, outputCIGAR.long, sum)
    
    pass.filter <-
      outputCIGAR.long.counts.per.locus[
        outputCIGAR.long.counts.per.locus[["count"]] >= minimum.count,
        
      ]
    
    outputCIGAR.long.filtered <-
      merge(
        outputCIGAR.long,
        pass.filter[, c("sample_id", "locus")],
        sort = FALSE
      )
    
    outputCIGAR.long.filtered
}


calculate.remaining.nloci <- function(outputCIGAR.long.filtered) {
  outputCIGAR.long.filtered.nloci.per.sample <-
    aggregate(
      locus ~ sample_id,
      outputCIGAR.long.filtered,
      function(x) length(unique(x))
    )
  names(outputCIGAR.long.filtered.nloci.per.sample)[2] <- "nloci"
  
  outputCIGAR.long.filtered.nloci.per.sample
}


plot.remaining.nloci <-
  function(outputCIGAR.long.filtered.nloci.per.sample, minimum.nloci) {
    plot(
      sort(outputCIGAR.long.filtered.nloci.per.sample[["nloci"]]),
      main = NA,
      xlab = NA,
      ylab = "Remaining number of locus",
      xaxt = "n"
    )
    abline(h = minimum.nloci, col = "red", lty = "dashed")
}


filter.failed.samples <-
  function(
    outputCIGAR.long.filtered,
    outputCIGAR.long.filtered.nloci.per.sample,
    minimum.nloci
  ) {
    pass.samples <-
      outputCIGAR.long.filtered.nloci.per.sample[
        outputCIGAR.long.filtered.nloci.per.sample[["nloci"]] >= minimum.nloci,
        "sample_id"
      ]
    
    outputCIGAR.long.filtered.subset <-
      outputCIGAR.long.filtered[
        outputCIGAR.long.filtered[["sample_id"]] %in% pass.samples,
        
      ]
    
    outputCIGAR.long.filtered.subset
  }


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to outputCIGAR.tsv file path
outputCIGAR.file <- "location/to/outputCIGAR.tsv"

outputCIGAR <- read.delim(outputCIGAR.file, check.names = FALSE)
outputCIGAR.long <- outputCIGAR.to.long(outputCIGAR)


# FIXME: change to path of unfiltered microhaplotype data
outputCIGAR.long.file <- "location/to/mhap_unfiltered.tsv"

write.table(
  outputCIGAR.long,
  outputCIGAR.long.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change read-pair threshold for failed genotype
minimum.count <- 10

outputCIGAR.long.filtered <-
  filter.outputCIGAR.long(outputCIGAR.long, minimum.count)

outputCIGAR.long.filtered.nloci.per.sample <-
  calculate.remaining.nloci(outputCIGAR.long.filtered)


# FIXME: change to path for statistics on remaining locus (text)
outputCIGAR.long.filtered.nloci.per.sample.text <- "location/to/mhap_filtered_nloci.tsv"

write.table(
  outputCIGAR.long.filtered.nloci.per.sample,
  outputCIGAR.long.filtered.nloci.per.sample.text,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change remaining locus threshold for failed samples
minimum.nloci <- 80

# FIXME: change to path for statistics on remaining locus (figure)
outputCIGAR.long.filtered.nloci.per.sample.figure <- "location/to/mhap_filtered_nloci.pdf"

pdf(file = outputCIGAR.long.filtered.nloci.per.sample.figure)
plot.remaining.nloci(outputCIGAR.long.filtered.nloci.per.sample, minimum.nloci)
dev.off()

outputCIGAR.long.filtered.subset <-
  filter.failed.samples(
    outputCIGAR.long.filtered,
    outputCIGAR.long.filtered.nloci.per.sample,
    minimum.nloci
  )

## STOP and look at the generated text and figure statistics
## Was the remaining locus threshold appropriate for this dataset?


# FIXME: change to path of filtered microhaplotype data
outputCIGAR.long.filtered.subset.file <- "location/to/mhap_filtered.tsv"

write.table(
  outputCIGAR.long.filtered.subset,
  outputCIGAR.long.filtered.subset.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
