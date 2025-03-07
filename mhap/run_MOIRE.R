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


# FIXME: change to path of polyclonal status
polyclonal_status_file <- "location/to/mhap_polyclonal_status.tsv"

polyclonal_status <-
  data.frame(
    sample_id = coi_summary[["sample_id"]],
    is_polyclonal = coi_summary[["prob_polyclonal"]] >= 0.5
  )

write.table(
  polyclonal_status,
  polyclonal_status_file,
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
