Before running the scripts, install 'moire', 'dcifer', 'readxl', 'ggnetwork', 'ape':

pkgs <- character()
for (package in c("moire", "dcifer", "readxl", "ggnetwork", "ape")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    pkgs <- c(pkgs, package)
  }
}
install.packages(pkgs, repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))

This might take around 15 minutes to install everything :)

Steps:
1. prepare_mhap_data.R
2. run_MOIRE.R
3. run_Dcifer.R
4. plot_MOIRE.R
5. plot_Dcifer.R
6. plot_PCoA_NJ.R

After each step, please restart R.


To evaluate the performance of a genetic marker panel in inferring relatedness, install 'paneljudge':

install.packages("devtools")
devtools::install_github("aimeertaylor/paneljudge", build_vignettes = TRUE)

After prepare_mhap_data.R, do run_paneljudge.R then plot_paneljudge.R.
