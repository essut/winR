#!/usr/bin/env Rscript
library(dcifer)

## Load all functions by running everything from line 5-164
# use case: alleles in dsmp are not present in the population allele 
pad.afreq <- function(afreq, dsmp) {
  # assumes the sample allele frequencies are more diverse
  smpafreq <- dsmp[[1]]
  
  # normalise number of locus in population
  locus.diff <- setdiff(names(smpafreq), names(afreq))
  afreq <- c(afreq, setNames(vector("list", length(locus.diff)), locus.diff))
  afreq <- afreq[order(names(afreq))]
  
  # normalise number of allele in population
  for (i in seq_along(smpafreq)) {
    allele.diff <- setdiff(names(smpafreq[[i]]), names(afreq[[i]]))
    
    # set the frequencies for the missing alleles
    small <- min(afreq[[i]]) / length(allele.diff)
    smalls <- setNames(rep(small, length(allele.diff)), allele.diff)
    
    afreq[[i]] <- c(afreq[[i]], smalls)
    afreq[[i]] <- afreq[[i]][order(names(afreq[[i]]))]
    
    # normalise frequencies to 1
    afreq[[i]] <- afreq[[i]] / sum(afreq[[i]])
  }
  
  return(afreq)
}


# old function to get relatedness estimate from ibdDat
get.m1.estimate <- function(dres0) {
  dmat <- dres0[, , "estimate"]
  m1.estimate <- as.data.frame(as.table(dmat))
  m1.estimate <- m1.estimate[!is.na(m1.estimate[["Freq"]]), ]
  names(m1.estimate) <- c("sample_id2", "sample_id1", "M1")
  
  return(m1.estimate)
}


# old function to analyse relatedness of "significantly related" pairs
analyse.significantly.related.samples <-
  function(dsmp, coi, afreq, dres0, alpha = 0.05) {
    isig <- which(dres0[, , "p_value"] <= alpha, arr.ind = TRUE)[, 2:1]
    
    # if no significant relatedness found
    if (nrow(isig) == 0) {
      return(NULL)
    }
    
    sig2 <- vector("list", nrow(isig))
    for (i in 1:nrow(isig)) {
      sig2[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, equalr = TRUE)
    }
    M2 <- sapply(sig2, length)
    rtotal2 <- sapply(sig2, sum)
    
    samples <- names(dsmp)
    sig <- data.frame(
      sample_id1 = samples[isig[, 1]],
      sample_id2 = samples[isig[, 2]],
      M = M2,
      rtotal = rtotal2
    )
    
    return(sig)
  }


# old function to retrieve relatedness of all pairs
calculate.overall.relatedness.estimate <- function(m1.estimate, sig, coi) {
  if (is.null(sig)) {
    mall.estimate <- cbind(m1.estimate, M = NA, rtotal = NA)
  } else {
    mall.estimate <- merge(m1.estimate, sig, all = TRUE)
  }
  
  coi.df <-
    data.frame(
      sample_id = sapply(strsplit(names(coi), ".", fixed = TRUE), "[[", 1),
      coi = coi,
      row.names = NULL
    )
  
  mall.estimate <-
    merge(mall.estimate, coi.df, by.x = "sample_id2", by.y = "sample_id")
  mall.estimate <-
    merge(
      mall.estimate,
      coi.df,
      by.x = "sample_id1",
      by.y = "sample_id",
      suffixes = c("2", "1")
    )
  
  mall.estimate[["scaled_r"]] <-
    mall.estimate[["rtotal"]] /
    (pmin(mall.estimate[["coi1"]], mall.estimate[["coi2"]]))
  
  mall.estimate[["relatedness"]] <- mall.estimate[["scaled_r"]]
  mall.estimate[is.na(mall.estimate[["relatedness"]]), "relatedness"] <-
    mall.estimate[is.na(mall.estimate[["relatedness"]]), "M1"]
  
  return(mall.estimate)
}


# nr is used to control the precision of the relatedness estimate
# where precision = 1 / nr, be warned that increasing this value
# will also increase the time and memory required to finish
analyse.all.pairs.relatedness <-
  function(dsmp, coi, afreq, spec, alpha = 0.05, nr = 1000) {
    # setup to run Dcifer in parallel
    cl <- parallel::makeCluster(spec)
    on.exit(parallel::stopCluster(cl = cl))
    parallel::clusterExport(
      cl = cl, c("dsmp", "coi", "afreq", "alpha", "nr"), envir = environment()
    )
    
    indices <- combn(1:length(dsmp), 2, simplify = FALSE)
    
    # do not limit number of related pairs (Mmax)
    res <-
      parallel::parLapply(
        cl = cl,
        indices,
        function(x) {
          dcifer::ibdEstM(
            dsmp[x],
            coi[x],
            afreq,
            Mmax = max(coi),
            confreg = TRUE,
            alpha = alpha,
            equalr = TRUE,
            nrs = nr
          )
        }
      )
    
    pairs <- vapply(indices, function(x) names(dsmp)[x], character(2))
    coi1 <- vapply(indices, function(x) coi[x[1]], numeric(1))
    coi2 <- vapply(indices, function(x) coi[x[2]], numeric(1))
    rhats <- lapply(res, function(x) x[["rhat"]])
    Ms <- lengths(rhats)
    rtotals <- vapply(rhats, sum, numeric(1))
    confints <- vapply(res, function(x) range(x[["confreg"]]), numeric(2))
    
    data.frame(
      sample_id1 = pairs[1, ],
      sample_id2 = pairs[2, ],
      coi1 = coi1,
      coi2 = coi2,
      M = Ms,
      rtotal = rtotals,
      lower_CI = confints[1, ],
      upper_CI = confints[2, ],
      scaled_r = rtotals / pmin(coi1, coi2)
    )
  }


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to folder to store outputs
output.dir <- "location/to/Dcifer"
dir.create(output.dir, recursive = TRUE)

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/data/mhap_filtered.tsv"

dlong <- read.delim(sfile)
dsmp <- formatDat(dlong, svar = "sample_id", lvar = "locus", avar = "allele")
coi   <- getCOI(dsmp)
afreq <- calcAfreq(dsmp, coi)

# FIXME: adjust the number of threads to use if necessary
# if not set, will try to use as many threads as possible
spec <- parallel::detectCores() - 2L

mall.estimate <- analyse.all.pairs.relatedness(dsmp, coi, afreq, spec)

mall.estimate.file <-
  paste0(output.dir, "/", "mhap_between_relatedness_estimate.tsv")

write.table(
  mall.estimate,
  mall.estimate.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
