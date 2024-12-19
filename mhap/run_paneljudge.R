#!/usr/bin/env Rscript
library(readxl)
library(paneljudge)

## Load all functions by running everything from line 6-239
create.markers <- function(dlong) {
  loci <- unique(dlong[["locus"]])
  Chr.Pos.Amplicon_name <- strsplit(loci, ":")
  
  Chr <- vapply(Chr.Pos.Amplicon_name, "[[", character(1), 1)
  Pos <- vapply(Chr.Pos.Amplicon_name, "[[", character(1), 2)
  Amplicon_name <- vapply(Chr.Pos.Amplicon_name, "[[", character(1), 3)
  
  chrom <- as.numeric(vapply(strsplit(Chr, "_"), "[[", character(1), 2))
  
  Start.Stop <- strsplit(Pos, "-")
  Start <- as.numeric(vapply(Start.Stop, "[[", character(1), 1))
  Stop <- as.numeric(vapply(Start.Stop, "[[", character(1), 2))
  
  pos <- rowMeans(cbind(Start, Stop))
  
  markers <- 
    data.frame(
      Amplicon_name,
      Chr,
      Start,
      Stop,
      length = Stop - Start,
      pos,
      chrom
    )
  
  markers <- markers[do.call(order, markers[, c("chrom", "Start", "Stop")]), ]
  row.names(markers) <- markers[["Amplicon_name"]]
  
  markers[["distances"]] <- c(diff(markers[["pos"]]), Inf)
  markers[c(diff(markers[["chrom"]]), 1) > 0, "distances"] <- Inf
  
  return(markers)
}


create.fs <- function(dlong, markers) {
  markers[["locus"]] <-
    paste0(
      markers[["Chr"]],
      ":",
      markers[["Start"]],
      "-",
      markers[["Stop"]],
      ":",
      markers[["Amplicon_name"]]
    )
  
  dlong.markers <- merge(dlong, markers, sort = FALSE)
  dlong.markers[["Amplicon_name"]] <-
    factor(
      dlong.markers[["Amplicon_name"]],
      levels = markers[["Amplicon_name"]]
    )
  
  allele.count.per.locus <-
    tapply(
      dlong.markers[["allele"]],
      dlong.markers[["Amplicon_name"]],
      table
    )
  
  allele.frequency.per.locus <- lapply(allele.count.per.locus, proportions)
  allele.frequency.per.locus <-
    lapply(allele.frequency.per.locus, sort, decreasing = TRUE)
  
  # pad allele frequencies to have the same length
  allele.frequency.per.locus <-
    lapply(
      allele.frequency.per.locus,
      function(x) {
        c(x, numeric(max(lengths(allele.frequency.per.locus)) - length(x)))
      }
    )
  
  fs <- do.call(rbind, allele.frequency.per.locus)
  colnames(fs) <- paste0("Allele.", seq_len(ncol(fs)))
  
  missing <- rownames(fs)[rowSums(fs) == 0]
  if (length(missing) > 0) {
    warning("Data unavailable for marker ", paste(missing, collapse = ", "))
    fs <- fs[!rownames(fs) %in% missing, ]
  }
  
  return(fs)
}


create.frequencies <- function(dlong, markers, metadata.group.column) {
  groups <- unique(dlong[[metadata.group.column]])
  
  frequencies <-
    lapply(
      groups,
      function(x) {
        create.fs(dlong[dlong[[metadata.group.column]] %in% x, ], markers)
      }
    )
  
  names(frequencies) <- groups
  
  return(frequencies)
}


# use the same markers across different groups
synchronise.used.markers <- function(frequencies) {
  markers <- NULL
  
  # find markers that overlap between all groups
  for (fs in frequencies) {
    if (is.null(markers)) {
      markers <- rownames(fs)
    } else {
      markers <- intersect(markers, rownames(fs))
    }
  }
  
  if (length(markers) == 0) {
    stop("No markers overlap between all different groups.")
  }
  
  # subset to all overlap markers
  for (i in seq_along(frequencies)) {
    fs <- frequencies[[i]]
    frequencies[[i]] <- fs[markers, ]
  }
  
  return(frequencies)
}


compute.diversities <- function(frequencies, metadata.group.column) {
  diversities <- list()
  
  for (group in names(frequencies)) {
    diversity <- compute_diversities(fs = frequencies[[group]])
    diversity <-
      data.frame(group, Amplicon_name = names(diversity), diversity = diversity)
    names(diversity)[1] <- metadata.group.column
    row.names(diversity) <- NULL
    
    diversities[[group]] <- diversity
  }
  diversities <- do.call(rbind, diversities)
  row.names(diversities) <- NULL
  
  return(diversities)
}


compute.eff.cardinalities <- function(frequencies, metadata.group.column) {
  eff.cardinalities <- list()
  
  for (group in names(frequencies)) {
    eff.cardinality <- compute_eff_cardinalities(fs = frequencies[[group]])
    eff.cardinality <-
      data.frame(group, Amplicon_name = names(eff.cardinality), eff_cardinality = eff.cardinality)
    names(eff.cardinality)[1] <- metadata.group.column
    row.names(eff.cardinality) <- NULL
    
    eff.cardinalities[[group]] <- eff.cardinality
  }
  eff.cardinalities <- do.call(rbind, eff.cardinalities)
  row.names(eff.cardinalities) <- NULL
  
  return(eff.cardinalities)
}


estimate.r.and.k <- function(frequencies, metadata.group.column, markers, k, r) {
  krhats <- list()
  
  for (group in names(frequencies)) {
    simulated_Ys <-
      simulate_Ys(frequencies[[group]], markers[["distances"]], k, r)
    
    krhat <-
      estimate_r_and_k(frequencies[[group]], markers[["distances"]], simulated_Ys)
    
    khat <- krhat['khat']
    rhat <- krhat['rhat']

    krhat <- data.frame(group, k, r, khat, rhat)
    names(krhat)[1] <- metadata.group.column
    row.names(krhat) <- NULL

    krhats[[group]] <- krhat
  }
  
  krhats <- do.call(rbind, krhats)
  row.names(krhats) <- NULL
  
  return(krhats)
}


compute.r.and.k.CIs <-
  function(frequencies, metadata.group.column, markers, krhats, confidence = 95) {
    kCIs <- vector("list", nrow(krhats))
    rCIs <- vector("list", nrow(krhats))
    
    for (i in 1:nrow(krhats)) {
      print(paste("Calculating CIs for row", i))
      
      group <- krhats[i, metadata.group.column]
      khat <- krhats[i, "khat"]
      rhat <- krhats[i, "rhat"]
      
      CIs <-
        compute_r_and_k_CIs(
          frequencies[[group]],
          markers[["distances"]],
          khat,
          rhat,
          confidence = confidence
        )
      
      kCI <-
        setNames(CIs['khat', ], paste("khat", names(CIs['khat', ]), sep = "_"))
      rCI <-
        setNames(CIs['rhat', ], paste("rhat", names(CIs['rhat', ]), sep = "_"))
      
      kCIs[[i]] <- t(kCI)
      rCIs[[i]] <- t(rCI)
    }
    
    kCIs <- do.call(rbind, kCIs)
    rCIs <- do.call(rbind, rCIs)
    krhats.CIs <- cbind(krhats, kCIs, rCIs)
    
    return(krhats.CIs)
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


# FIXME: change according sample and group column in metadata
metadata.sample.column <- "ID"
metadata.group.column <- "Year"

dlong <-
  merge(
    dlong,
    metadata[, c(metadata.sample.column, metadata.group.column)],
    by.x = "sample_id",
    by.y = metadata.sample.column
  )


markers <- create.markers(dlong)
frequencies <- create.frequencies(dlong, markers, metadata.group.column)

## check how many markers are used for all groups
frequencies <- synchronise.used.markers(frequencies)
nrow(frequencies[[1]])


# FIXME: change to path of diversities
diversities.file <- "location/to/mhap_paneljudge_diversities.tsv"

diversities <- compute.diversities(frequencies, metadata.group.column)
write.table(
  diversities,
  diversities.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of effective cardinalities
eff.cardinalities.file <- "location/to/mhap_paneljudge_eff_cardinalities.tsv"

eff.cardinalities <- compute.eff.cardinalities(frequencies, metadata.group.column)
write.table(
  eff.cardinalities,
  eff.cardinalities.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of estimates of k and r
krhats.file <- "location/to/mhap_paneljudge_k_r_estimates.tsv"

# FIXME: number of pairs to simulate
n <- 100

# FIXME: switch rates (k) and relatedness values (r) parameter to simulate
ks <- c(5)
rs <- c(0.01, seq(0.25, 0.75, 0.25), 0.99)

krhats <- vector("list", n * length(ks) * length(rs))
j <- 1

for (i in 1:n) {
  for (k in ks) {
    for (r in rs) {
      print(paste("i =", i, "k =", k, "r =", r))
      
      krhat <-
        estimate.r.and.k(frequencies, metadata.group.column, markers, k, r)
      
      krhats[[j]] <- cbind(Pair = i, krhat)
      j <- j + 1
    }
  }
}

krhats <- do.call(rbind, krhats)

write.table(
  krhats,
  krhats.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


## perform the lines below if confidence intervals are needed
## this might take a while to perform

# FIXME: change to path of estimates of k and r with confidence intervals
krhats.CIs.file <- "location/to/mhap_paneljudge_k_r_estimates_CIs.tsv"

krhats.CIs <- compute.r.and.k.CIs(frequencies, metadata.group.column, markers, krhats)

write.table(
  krhats.CIs,
  krhats.CIs.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
