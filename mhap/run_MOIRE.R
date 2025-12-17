#!/usr/bin/env Rscript

## Load all functions by running everything from line 4-61
plot.posterior.burnin.trace <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]

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
    chain <- mcmc_results[["chains"]][[i]]

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
    chain <- mcmc_results[["chains"]][[i]]

    plot(
      density(chain[["posterior_sample"]], bw = "SJ"),
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) density"
    )
  }
}

plot.posterior.sample.acf <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]

    acf <-
      acf(
        chain[["posterior_sample"]],
        lag.max = length(chain[["posterior_sample"]]) / 5,
        plot = FALSE
      )

    plot(
      acf,
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) ACF"
    )
  }
}


## Change FIXME lines to the appropriate parameters

## If you already have the mcmc_results saved as RDS,
## you don't need to re-do run_mcmc (use readRDS() instead)

# FIXME: change to folder to store outputs
output_dir <- "location/to/MOIRE"
dir.create(output_dir, recursive = TRUE)

# FIXME: change to path of filtered microhaplotype data
sfile <- "location/to/data/mhap_filtered.tsv"

df <- read.delim(sfile)
data <- moire::load_long_form_data(df)


# these default parameters can be adjusted accordingly:
# burn-in = 1,000 [burnin], samples per chain = 1,000 [samples_per_chain]
# parallel tempering chains = 40 chains between 0.5 to 1 [pt_chains]
# adjust the number of threads [pt_num_threads] accordingly
burnin <- 1000
pt_num_threads <- 40
pt_chains <- seq(0.5, 1, length.out = 40)
samples_per_chain <- 1000

mcmc_results <-
  moire::run_mcmc(
    data = data,
    is_missing = data[["is_missing"]],
    burnin = burnin,
    samples_per_chain = samples_per_chain,
    r_alpha = 3,
    r_beta = 3,
    pt_chains = pt_chains,
    pt_num_threads = pt_num_threads
  )

file <- paste0(output_dir, "/", "mhap_MOIRE.rds")

saveRDS(mcmc_results, file = file)


# Ideally, all swap_rate should be above 0.25
chain_swaps_plot_file <-
  paste0(output_dir, "/", "mhap_MOIRE_chain_swaps.pdf")
pdf(file = chain_swaps_plot_file, width = 3, height = 2)

print(moire::plot_chain_swaps(mcmc_results))
dev.off()


## produce MCMC diagnostics for assessment, see:
## https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics
nchain <- length(mcmc_results[["chains"]])
m <- ceiling(sqrt(nchain))

posterior_burnin_trace_plot_file <-
  paste0(output_dir, "/", "mhap_MOIRE_posterior_burnin_trace.pdf")
pdf(file = posterior_burnin_trace_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.burnin.trace(mcmc_results)
dev.off()

posterior_sample_trace_plot_file <-
  paste0(output_dir, "/", "mhap_MOIRE_posterior_sample_trace.pdf")
pdf(file = posterior_sample_trace_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.sample.trace(mcmc_results)
dev.off()

posterior_sample_acf_plot_file <-
  paste0(output_dir, "/", "mhap_MOIRE_posterior_sample_ACF.pdf")
pdf(file = posterior_sample_acf_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.sample.acf(mcmc_results)
dev.off()

## STOP and determine if the MCMC parameters need to be re-adjusted
## https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics
## section: Trace plots & Autocorrelation function (ACF) plots

coi_summary_file <- paste0(output_dir, "/", "mhap_COI_summary.tsv")

coi_summary <- moire::summarize_coi(mcmc_results)
write.table(
  coi_summary,
  coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


polyclonal_status_file <- paste0(output_dir, "/", "mhap_polyclonal_status.tsv")

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


effective_coi_summary_file <- paste0(
  output_dir,
  "/",
  "mhap_effective_COI_summary.tsv"
)

effective_coi_summary <- moire::summarize_effective_coi(mcmc_results)
write.table(
  effective_coi_summary,
  effective_coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
