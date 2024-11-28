#!/usr/bin/env Rscript

## Change FIXME lines to the appropriate parameters

## If you already have the RDS from previous MOIRE MCMC,
## skip until line 22

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/mhap_filtered.tsv"

df <- read.delim(sfile)
data <- moire::load_long_form_data(df)
mcmc_results <- moire::run_mcmc(data, is_missing = data[["is_missing"]])


# FIXME: change to path of MOIRE MCMC 
file <- "location/to/mhap_MOIRE.rds"

saveRDS(mcmc_results, file = file)


# FIXME: change to path of COI summary
coi_summary_file <- "location/to/mhap_COI_summary.tsv"

coi_summary <- moire::summarize_coi(mcmc_results)
write.table(
  coi_summary,
  coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of effective COI summary
effective_coi_summary_file <- "location/to/mhap_effective_COI_summary.tsv"

effective_coi_summary <- moire::summarize_effective_coi(mcmc_results)
write.table(
  effective_coi_summary,
  effective_coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of locus heterozygosity summary
he_summary_file <- "location/to/mhap_he_summary.tsv"

he_summary <- moire::summarize_he(mcmc_results)
write.table(
  he_summary,
  he_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of allele frequencies summary
allele_freqs_summary_file <- "location/to/mhap_allele_freqs_summary.tsv"

allele_freqs_summary <- moire::summarize_allele_freqs(mcmc_results)
write.table(
  allele_freqs_summary,
  allele_freqs_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of within-infection relatedness estimate
relatedness_summary_file <- "location/to/mhap_within_relatedness_summary.tsv"

relatedness_summary <- moire::summarize_relatedness(mcmc_results)
write.table(
  relatedness_summary,
  relatedness_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
