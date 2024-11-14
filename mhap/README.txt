Before running the scripts, install 'dcifer', 'RColorBrewer', 'ape', 'network', 'sna', 'readxl', 'moire', 'GGally':

pkgs <- character()
for (package in c("dcifer", "RColorBrewer", "ape", "network", "sna", "readxl", "moire", "GGally")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    pkgs <- c(pkgs, package)
  }
}
install.packages(pkgs, repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))

This might take around ten minutes to install everything :)

Steps:
1. prepare_mhap_data.R
2. run_MOIRE.R
3. run_Dcifer.R
4. plot_MOIRE.R
5. plot_Dcifer.R
6. plot_PCoA_NJ.R

After each step, please restart R.
