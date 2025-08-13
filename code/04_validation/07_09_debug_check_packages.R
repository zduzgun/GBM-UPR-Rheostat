# 07_09_debug_check_packages.R
# Amaç: Gerekli paketlerin yüklü olup olmadýðýný saniyeler içinde kontrol etmek.

cat("--- PACKAGE CHECK STARTING ---\n")

# clusterProfiler paketini yüklemeyi dene
cat("Checking for 'clusterProfiler'...\n")
load_cp <- try(suppressPackageStartupMessages(library(clusterProfiler)), silent = TRUE)
if (inherits(load_cp, "try-error")) {
  stop("FATAL: 'clusterProfiler' package is NOT installed correctly.")
} else {
  cat("SUCCESS: 'clusterProfiler' loaded.\n")
}

# org.Hs.eg.db paketini yüklemeyi dene
cat("Checking for 'org.Hs.eg.db'...\n")
load_orgdb <- try(suppressPackageStartupMessages(library(org.Hs.eg.db)), silent = TRUE)
if (inherits(load_orgdb, "try-error")) {
  stop("FATAL: 'org.Hs.eg.db' annotation package is NOT installed correctly.")
} else {
  cat("SUCCESS: 'org.Hs.eg.db' loaded.\n")
}

cat("\n--- PACKAGE CHECK COMPLETED SUCCESSFULLY ---\n")