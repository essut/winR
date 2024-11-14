#!/usr/bin/env Rscript
library(dcifer)

## Change FIXME lines to the appropriate parameters

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/mhap_filtered.tsv"

dlong <- read.delim(sfile)
dsmp <- formatDat(dlong, svar = "sample_id", lvar = "locus", avar = "allele")
coi   <- getCOI(dsmp)
afreq <- calcAfreq(dsmp, coi)


## only one pair of strains between two infections can be related
dres0 <- ibdDat(dsmp, coi, afreq)

dmat <- dres0[, , "estimate"]
m1.estimate <- as.data.frame(as.table(dmat))
m1.estimate <- m1.estimate[!is.na(m1.estimate[["Freq"]]), ]
names(m1.estimate) <- c("sample_id1", "sample_id2", "M1")


## allow multiple pairs of strains to be related between two infections
alpha <- 0.05
isig <- which(dres0[, , "p_value"] <= alpha, arr.ind = TRUE)[, 2:1] 

revals <- mapply(generateReval, 1:5, nr = c(1e3, 1e2, 32, 16, 12))

sig2 <- vector("list", nrow(isig))
for (i in 1:nrow(isig)) {
  sig2[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, equalr = TRUE)
}
M2 <- sapply(sig2, length)
rtotal2 <- sapply(sig2, sum)

samples <- names(dsmp)
sig <- data.frame(
  sample_id1 = samples[isig[, 2]],
  sample_id2 = samples[isig[, 1]],
  M = M2,
  rtotal = rtotal2
)


## combine results from restrained and unrestrained M
mall.estimate <- merge(m1.estimate, sig, all = TRUE)
coi.df <-
  data.frame(
    sample_id = sapply(strsplit(names(coi), ".", fixed = TRUE), "[[", 1),
    coi = coi,
    row.names = NULL
  )
mall.estimate <-
  merge(mall.estimate, coi.df, by.x = "sample_id1", by.y = "sample_id")
mall.estimate <-
  merge(mall.estimate, coi.df, by.x = "sample_id2", by.y = "sample_id")

mall.estimate[["scaled_r"]] <-
  mall.estimate[["rtotal"]] /
  (pmin(mall.estimate[["coi.x"]], mall.estimate[["coi.y"]]))

mall.estimate[["relatedness"]] <- mall.estimate[["scaled_r"]]
mall.estimate[is.na(mall.estimate[["relatedness"]]), "relatedness"] <-
  mall.estimate[is.na(mall.estimate[["relatedness"]]), "M1"]


# FIXME: change to path of between-infection relatedness estimate
mall.estimate.file <- "location/to/mhap_relatedness_estimate.tsv"

write.table(
  mall.estimate,
  mall.estimate.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
