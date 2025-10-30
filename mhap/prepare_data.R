#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-513
load.outputs <- function(output.files) {
  outputs <- list()
  for (output.file in output.files) {
    outputs[[output.file]] <- read.delim(output.file, check.names = FALSE)
  }
  outputs
}


get.sample.list <- function(outputs) {
  sample.list <- data.frame()

  for (i in names(outputs)) {
    sample.list <-
      rbind(
        sample.list,
        data.frame(
          sample_id = outputs[[i]][["sample"]],
          total_read_pairs = rowSums(outputs[[i]][, 2:ncol(outputs[[i]])]),
          file = i
        )
      )
  }

  sample.list[["first"]] <- !duplicated(sample.list[["sample_id"]])

  # make sure the row names are valid
  row.names(sample.list) <-
    make.unique(paste0(
      sample.list[["sample_id"]],
      sample.list[["total_read_pairs"]]
    ))

  sample.list[["maximum"]] <- FALSE
  sample.list[
    row.names(sample.list) %in%
      do.call(
        paste0,
        aggregate(total_read_pairs ~ sample_id, sample.list, max)
      ),
    "maximum"
  ] <- TRUE

  sample.list <- sample.list[order(sample.list[["sample_id"]]), ]
  row.names(sample.list) <- NULL

  sample.list
}


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


merge.outputs <- function(outputs, sample.list, how) {
  sample.list <-
    switch(
      how,
      sum = {
        sample.list
      },
      maximum = {
        sample.list[sample.list[["maximum"]], ]
      },
      first = {
        sample.list[sample.list[["first"]], ]
      },
      stop("Valid options are 'sum', 'maximum', 'first'")
    )
  sample.list <- sample.list[, c("sample_id", "total_read_pairs", "file")]

  longs <- list()

  for (output.file in names(outputs)) {
    working.list <- sample.list[sample.list[["file"]] %in% output.file, ]

    if (nrow(working.list) == 0) {
      next
    }

    output <- outputs[[output.file]]
    output <- output[output[["sample"]] %in% working.list[["sample_id"]], ]

    longs[[output.file]] <- output.to.long(output)
  }

  long <- do.call(rbind, longs)

  # check if user supplied different allele formats together
  is.outputCIGAR <- any(grepl("(\\.|[0-9][ACGT])", long[["allele"]]))
  is.outputHaplotypes <- any(grepl("[:*]", long[["allele"]]))

  if (is.outputCIGAR & is.outputHaplotypes) {
    stop(
      "Detected both pseudoCIGAR and cs tag, please do not mix outputCIGAR.tsv and outputHaplotypes.tsv together"
    )
  }

  # consolidate allele counts from different runs
  list(
    long = aggregate(count ~ sample_id + locus + allele, long, sum),
    sample.list = merge(
      aggregate(total_read_pairs ~ sample_id, sample.list, sum),
      aggregate(file ~ sample_id, sample.list, function(x) {
        paste(x, collapse = ",")
      })
    )
  )
}


.rmindel.allele.outputCIGAR <- function(long) {
  # insertion pattern: {position}I={bases}
  # deletion pattern: {position}D={bases}
  long[["allele"]] <- gsub("([0-9]+)[DI]=([ACGT]+)", "", long[["allele"]])

  # default to wild type if all variants were removed
  long[long[["allele"]] %in% "", "allele"] <- "."

  long
}


.rmindel.allele.outputHaplotypes <- function(long) {
  # outputHaplotypes.tsv have a suffix format to determine
  # whether to keep or remove soft clips at each end
  trim.start <- grepl("_s-.*$", long[["allele"]])
  trim.end <- grepl("-e$", long[["allele"]])

  # trim the suffix
  long[["allele"]] <- sub("_.*$", "", long[["allele"]])

  # trim the soft clips, pattern: [+-][ACGT]+
  long[trim.start, "allele"] <- sub(
    "^[+-][ACGT]+",
    "",
    long[trim.start, "allele"]
  )
  long[trim.end, "allele"] <- sub("[+-][ACGT]+$", "", long[trim.end, "allele"])

  # insertion pattern: +[acgt]+
  long[["allele"]] <- gsub("\\+[acgt]+", "", long[["allele"]])

  ## removing deletion is more involved in cs tag
  # deletion pattern: -[acgt]+
  del.match <- gregexpr("-[acgt]+", long[["allele"]])

  # find how long the deletions are
  del.sub <- lapply(del.match, function(x) {
    paste0(":", attr(x, "match.length") - 1)
  })
  del.sub <- lapply(del.sub, function(x) sub(":-2", "", x))

  # mark and replace the deletions with the correct reference length
  long[["allele"]] <- gsub("-[acgt]+", "|", long[["allele"]])
  long[["allele"]] <- strsplit(long[["allele"]], "\\|")
  del.sub <-
    mapply(
      function(x, y) c(x, character(y)),
      x = del.sub,
      y = lengths(long[["allele"]]) - lengths(del.sub)
    )
  long[["allele"]] <-
    mapply(
      function(x, y) paste0(c(rbind(x, y)), collapse = ""),
      x = long[["allele"]],
      y = del.sub
    )

  # scan and sum up the split references
  cs.modified <- character(length(long[["allele"]]))
  i <- 1

  for (cs.original in strsplit(long[["allele"]], character())) {
    numeric.store <- numeric()
    potential.numbers <- character()
    cs <- c()

    for (char in cs.original) {
      # if not numeric
      if (is.na(as.numeric(char))) {
        if (length(potential.numbers) != 0) {
          numeric.store <-
            c(
              numeric.store,
              as.numeric(paste0(potential.numbers, collapse = ""))
            )

          numeric.store <- na.omit(numeric.store)
        }

        # if not reference
        if (char != ":") {
          if (length(numeric.store) > 0) {
            cs <- c(cs, sum(numeric.store))
          }
          numeric.store <- numeric()
          potential.numbers <- character()
        }

        # store digits if there are other digits in storage
        if (length(numeric.store) == 0) {
          cs <- c(cs, char)
        } else {
          potential.numbers <- character()
        }
      } else {
        potential.numbers <- c(potential.numbers, char)
      }
    }

    # end of loop, make sure everything is written
    numeric.store <-
      sum(
        numeric.store,
        as.numeric(paste0(potential.numbers, collapse = "")),
        na.rm = TRUE
      )

    if (numeric.store != 0) {
      cs <- c(cs, as.character(numeric.store))
    }

    cs.modified[i] <- paste0(cs, collapse = "")
    i <- i + 1
  }

  long[["allele"]] <- cs.modified

  long
}


rmindel.allele <- function(long) {
  is.outputCIGAR <- any(grepl("[DI]=", long[["allele"]]))
  # also remove the indel information at the end
  is.outputHaplotypes <- any(grepl("([+-][acgt]+|_)", long[["allele"]]))

  if (is.outputCIGAR & is.outputHaplotypes) {
    stop(
      "Detected both pseudoCIGAR and cs tag, please do not mix outputCIGAR.tsv and outputHaplotypes.tsv together"
    )
  }
  if (is.outputCIGAR) {
    long <- .rmindel.allele.outputCIGAR(long)
  }
  if (is.outputHaplotypes) {
    long <- .rmindel.allele.outputHaplotypes(long)
  }

  # consolidate allele counts after indel removal
  aggregate(count ~ sample_id + locus + allele, long, sum)
}


select.markers <- function(long, keep.marker = "microhaplotype") {
  switch(
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


.calculate.prop.major <- function(long) {
  dlong.count <- aggregate(count ~ sample_id + locus, long, max)
  names(dlong.count)[length(dlong.count)] <- "major"

  long.wprop <- merge(long, dlong.count, sort = FALSE)
  long.wprop[["prop"]] <-
    long.wprop[["count"]] / long.wprop[["major"]]

  long.wprop
}


.calculate.prop.total <- function(long) {
  dlong.count <- aggregate(count ~ sample_id + locus, long, sum)
  names(dlong.count)[length(dlong.count)] <- "total"

  long.wprop <- merge(long, dlong.count, sort = FALSE)
  long.wprop[["prop"]] <-
    long.wprop[["count"]] / long.wprop[["total"]]

  long.wprop
}


allele.proportion.filter <-
  function(long, allele.proportion.cutoff, on = "total", how = "naive") {
    if (on %in% "total") {
      long.wprop <- .calculate.prop.total(long)
    } else if (on %in% "major") {
      long.wprop <- .calculate.prop.major(long)
    } else {
      stop("Valid options are 'total' or 'major' allele")
    }

    long.wprop.max <- aggregate(prop ~ locus + allele, long.wprop, max)
    names(long.wprop.max)[3] <- "prop.max"
    long.wprop <- merge(long.wprop, long.wprop.max)

    if (how %in% "naive") {
      how.col <- "prop"
    } else if (how %in% "population-aware") {
      how.col <- "prop.max"
    } else {
      stop("Valid options are 'naive' or 'population-aware'")
    }

    long.wprop[
      long.wprop[[how.col]] >= allele.proportion.cutoff,
      c("sample_id", "locus", "allele", "count")
    ]
  }


allele.count.filter <- function(long, allele.count.cutoff) {
  long[long[["count"]] >= allele.count.cutoff, ]
}


minimum.total.filter <- function(long, minimum.total.cutoff) {
  long.counts.per.locus <- aggregate(count ~ sample_id + locus, long, sum)

  pass.filter <-
    long.counts.per.locus[
      long.counts.per.locus[["count"]] >= minimum.total.cutoff,
    ]

  merge(long, pass.filter[, c("sample_id", "locus")], sort = FALSE)
}


calculate.allele.statistics <- function(long) {
  chrom.pos <- strsplit(long[["locus"]], ":")
  chrom <- vapply(chrom.pos, "[[", character(1), 1)
  pos <- vapply(chrom.pos, "[[", character(1), 2)
  marker <- vapply(chrom.pos, "[[", character(1), 3)

  start.end <- strsplit(pos, "-")
  start <- as.numeric(vapply(start.end, "[[", character(1), 1))
  end <- as.numeric(vapply(start.end, "[[", character(1), 2))

  long <-
    data.frame(
      chromosome = chrom,
      start.pos = start,
      end.pos = end,
      marker = marker,
      long
    )

  long.allele.statistics <-
    aggregate(
      sample_id ~ chromosome + start.pos + end.pos + marker + locus + allele,
      long,
      length
    )
  names(long.allele.statistics)[7] <- "nsample"

  long.allele.count <- as.data.frame(table(long[["locus"]]))
  names(long.allele.count) <- c("locus", "total.allele")

  long.allele.statistics <-
    merge(long.allele.statistics, long.allele.count, sort = FALSE)

  long.allele.statistics["prop.allele"] <-
    long.allele.statistics[["nsample"]] /
    long.allele.statistics[["total.allele"]]

  long.allele.statistics
}


plot.prop.allele <-
  function(long.allele.statistics, prop.allele.cutoff = -Inf) {
    hist(
      long.allele.statistics[["prop.allele"]],
      breaks = "FD",
      right = FALSE,
      main = "Proportion of alleles",
      xlab = NULL
    )

    abline(v = prop.allele.cutoff, col = "red", lty = "dashed")
  }


plot.nsample.per.allele <-
  function(long.allele.statistics, nsample.per.allele.cutoff = -Inf) {
    hist(
      long.allele.statistics[["nsample"]],
      breaks = "FD",
      right = FALSE,
      main = "Number of sample per allele",
      xlab = NULL
    )

    abline(v = nsample.per.allele.cutoff + 1, col = "red", lty = "dashed")
  }


prop.allele.filter <-
  function(long, long.allele.statistics, prop.allele.cutoff) {
    pass.prop.allele <-
      long.allele.statistics[
        long.allele.statistics[["prop.allele"]] >= prop.allele.cutoff,
        c("locus", "allele")
      ]

    long <- merge(long, pass.prop.allele, sort = FALSE)
    long[, c("sample_id", "locus", "allele", "count")]
  }


nsample.per.allele.filter <-
  function(long, long.allele.statistics, nsample.per.allele.cutoff) {
    pass.nsample.per.allele <-
      long.allele.statistics[
        long.allele.statistics[["nsample"]] >= nsample.per.allele.cutoff,
        c("locus", "allele")
      ]

    long <- merge(long, pass.nsample.per.allele, sort = FALSE)
    long[, c("sample_id", "locus", "allele", "count")]
  }


minimum.total.filter <- function(long, minimum.total.cutoff) {
  long.counts.per.locus <- aggregate(count ~ sample_id + locus, long, sum)

  pass.filter <-
    long.counts.per.locus[
      long.counts.per.locus[["count"]] >= minimum.total.cutoff,
    ]

  merge(long, pass.filter[, c("sample_id", "locus")], sort = FALSE)
}


calculate.remaining.nloci <- function(long) {
  long.nloci.per.sample <-
    aggregate(locus ~ sample_id, long, function(x) length(unique(x)))

  names(long.nloci.per.sample)[2] <- "nloci"

  long.count.per.sample <-
    aggregate(count ~ sample_id, long, sum)

  names(long.count.per.sample)[2] <- "filtered_read_pairs"

  merge(long.nloci.per.sample, long.count.per.sample)
}


plot.remaining.nloci <- function(long.nloci.per.sample, minimum.nloci) {
  plot(
    sort(long.nloci.per.sample[["nloci"]]),
    main = NA,
    xlab = NA,
    ylab = "Remaining number of locus",
    xaxt = "n"
  )

  abline(h = minimum.nloci, col = "red", lty = "dashed")
}


filter.failed.samples <-
  function(long, long.nloci.per.sample, minimum.nloci) {
    pass.samples <-
      long.nloci.per.sample[
        long.nloci.per.sample[["nloci"]] >= minimum.nloci,
        "sample_id"
      ]

    long[long[["sample_id"]] %in% pass.samples, ]
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


outputs <- load.outputs(output.files)

# use outputs for complete sample list
sample.list <- get.sample.list(outputs)

print(sample.list)
## STOP and check the sample list, take notes on what needs to be changed

# FIXME: choose an action for samples with the same name in different runs
# - "sum" the read pairs from different runs
# - pick the sample with "maximum" read pairs
# - pick the "first" available sample option (dependent on output.files order)
merged <- merge.outputs(outputs, sample.list, how = "sum")
long <- merged[["long"]]
sample.list <- merged[["sample.list"]]


# FIXME: example on how to change sample IDs based on above notes
# long[long[["sample_id"]] %in% "Human-AT2", "sample_id"] <- "NC"
# long[long[["sample_id"]] %in% "NTC-1xTE", "sample_id"] <- "EMPTY"
# long[long[["sample_id"]] %in% "Pf-K1", "sample_id"] <- "NC"
# long <- aggregate(count ~ sample_id + locus + allele, long, sum)

# FIXME: remove samples if needed (e.g. samples from a different cohort)
removed.samples <- c("sWGA_1_10", "sWGA_1_20", "sWGA_1_30")

sample.list <-
  sample.list[
    !sample.list[["sample_id"]] %in% removed.samples,
    ,
    drop = FALSE
  ]

long <- long[!long[["sample_id"]] %in% removed.samples, ]

long.file <- paste0(output.dir, "/", "long_unfiltered.tsv")

write.table(long, long.file, quote = FALSE, sep = "\t", row.names = FALSE)


mhap.filtered <- select.markers(long, keep.marker = "microhaplotype")

## FIXME: comment/uncomment any lines below to disable/enable filters
## FIXME: also, feel free to reorder filters as needed

# FIXME: remove insertions and deletions in alleles
mhap.filtered <- rmindel.allele(mhap.filtered)

# FIXME: filter locus with less than X read-pairs per sample
# minimum.total.cutoff <- 10 # at least 10 read-pairs in total per locus per sample
# mhap.filtered <- minimum.total.filter(mhap.filtered, minimum.total.cutoff)

# FIXME: filter alleles with less than X% proportion per sample
# FIXME: in relation to total allele (on = "total") or major allele (on = "major")
allele.proportion.cutoff <- 0.02 # at least 2% allele proportion per sample
allele.proportion.on <- "total" # per locus instead of per major allele
mhap.filtered <-
  allele.proportion.filter(
    mhap.filtered,
    allele.proportion.cutoff,
    on = allele.proportion.on,
    how = "naive"
  )

# FIXME: filter alleles with less than X read-pairs per sample
allele.count.cutoff <- 5 # at least 5 read-pairs needed for each allele per sample
mhap.filtered <- allele.count.filter(mhap.filtered, allele.count.cutoff)


mhap.filtered.allele.statistics <- calculate.allele.statistics(mhap.filtered)

mhap.filtered.allele.statistics.text <-
  paste0(output.dir, "/", "mhap_filtered_allele_statistics.tsv")

write.table(
  mhap.filtered.allele.statistics,
  mhap.filtered.allele.statistics.text,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

mhap.filtered.prop.allele.figure <-
  paste0(output.dir, "/", "mhap_filtered_prop_allele.pdf")

pdf(file = mhap.filtered.prop.allele.figure)
plot.prop.allele(mhap.filtered.allele.statistics)
dev.off()

mhap.filtered.nsample.per.allele.figure <-
  paste0(output.dir, "/", "mhap_filtered_nsample_per_allele.pdf")

pdf(file = mhap.filtered.nsample.per.allele.figure)
plot.nsample.per.allele(mhap.filtered.allele.statistics)
dev.off()

# FIXME: filter alleles with less than X% proportion per marker
# prop.allele.cutoff <- 0.01 # alleles contribute to at least 1% of total alleles
# mhap.filtered <-
#   prop.allele.filter(mhap.filtered, mhap.filtered.allele.statistics, prop.allele.cutoff)

# FIXME: filter alleles with less than X samples
# nsample.per.allele.cutoff <- 2 # alleles observed in at least 2 samples
# mhap.filtered <-
#   nsample.per.allele.filter(
#     mhap.filtered,
#     mhap.filtered.allele.statistics,
#     nsample.per.allele.cutoff
#   )

## END of filtering functions

mhap.filtered.nloci.per.sample <- calculate.remaining.nloci(mhap.filtered)

mhap.filtered.nloci.per.sample <-
  merge(mhap.filtered.nloci.per.sample, sample.list, all = TRUE)
mhap.filtered.nloci.per.sample[
  is.na(mhap.filtered.nloci.per.sample[["nloci"]]),
  "nloci"
] <- 0
mhap.filtered.nloci.per.sample[
  is.na(mhap.filtered.nloci.per.sample[["filtered_read_pairs"]]),
  "filtered_read_pairs"
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
  mhap.filtered.subset[
    !mhap.filtered.subset[["sample_id"]] %in% unused.samples,
  ]

mhap.filtered.subset.file <- paste0(output.dir, "/", "mhap_filtered.tsv")

write.table(
  mhap.filtered.subset,
  mhap.filtered.subset.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
