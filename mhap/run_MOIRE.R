#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-46
plot.posterior.burnin.trace <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    plot(
      chain[["posterior_burnin"]],
      type = "l",
      main = paste0("Chain ", i),
      xlab = "Iteration",
      ylab = "Posterior (burn-in)"
    )
  }
}

plot.posterior.sample.trace <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    plot(
      chain[["posterior_sample"]],
      type = "l",
      main = paste0("Chain ", i),
      xlab = "Iteration",
      ylab = "Posterior (sample)"
    )
  }
}

plot.posterior.sample.density <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    plot(
      density(chain[["posterior_sample"]], bw = "SJ"),
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) density"
    )
  }
}

plot.posterior.sample.acf <- function(mcmc_results) {
  for (chain in mcmc_results[["chains"]]) {
    acf(
      chain[["posterior_sample"]],
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) ACF"
    )
  }
}


## Change FIXME lines to the appropriate parameters

## If you already have the mcmc_results saved as RDS,
## you don't need to re-do run_mcmc (use readRDS() instead)

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/mhap_filtered.tsv"

df <- read.delim(sfile)
data <- moire::load_long_form_data(df)
mcmc_results <- moire::run_mcmc(data, is_missing = data[["is_missing"]])

# FIXME: change to path of MOIRE MCMC 
file <- "location/to/mhap_MOIRE.rds"

saveRDS(mcmc_results, file = file)


## produce MCMC diagnostics for assessment, see:
## https://sbfnk.github.io/mfiidd/mcmc_diagnostics.html

nchain <- length(mcmc_results[["chains"]])
m <- ceiling(sqrt(nchain))

# FIXME: change to path of posterior (burn-in) trace plot
posterior.burnin.trace.plot.file <- "location/to/mhap_MOIRE_burnin_trace.pdf"
pdf(file = posterior.burnin.trace.plot.file, width = 7, height = 7)
par(mfrow = c(m, m))
plot.posterior.burnin.trace(mcmc_results)
dev.off()

# FIXME: change to path of posterior (sample) trace plot
posterior.sample.trace.plot.file <- "location/to/mhap_MOIRE_sample_trace.pdf"
pdf(file = posterior.sample.trace.plot.file, width = 7, height = 7)
par(mfrow = c(m, m))
plot.posterior.sample.trace(mcmc_results)
dev.off()

# FIXME: change to path of posterior (sample) density plot
posterior.sample.density.plot.file <- "location/to/mhap_MOIRE_sample_density.pdf"
pdf(file = posterior.sample.density.plot.file, width = 7, height = 7)
par(mfrow = c(m, m))
plot.posterior.sample.density(mcmc_results)
dev.off()

# FIXME: change to path of posterior (sample) acf plot
posterior.sample.acf.plot.file <- "location/to/mhap_MOIRE_sample_acf.pdf"
pdf(file = posterior.sample.acf.plot.file, width = 7, height = 7)
par(mfrow = c(m, m))
plot.posterior.sample.acf(mcmc_results)
dev.off()

## STOP and determine whether the MCMC mixes well


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
