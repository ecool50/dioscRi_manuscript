# Fetch the manuscript data archive from Zenodo into $BIOHEART_ROOT/data/.
#
# Single 7.2 GB zip (Zenodo record 15694581, "Analysis files for the dioscRi
# manuscript"). MD5 is verified after download.
#
# Usage:
#   Rscript revision/scripts/helpers/fetch_data.R
#   Rscript revision/scripts/helpers/fetch_data.R --force      # re-download
#   BIOHEART_ROOT=/custom/path Rscript revision/scripts/helpers/fetch_data.R

ZENODO_FILE_URL <- "https://zenodo.org/api/records/15694581/files/dioscRi_analysis_data.zip/content"
ZIP_NAME        <- "dioscRi_analysis_data.zip"
ZIP_MD5         <- "621fbb0d58a88e47b6019add1ce33412"
ZIP_SIZE_GB     <- 7.2
SENTINEL_SUBDIR <- "dioscRi_analysis_data"  # produced by the unzip

args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

bioheart_root <- Sys.getenv(
  "BIOHEART_ROOT",
  unset = file.path(Sys.getenv("HOME"),
                    "Documents/Academic/PhD/bioheart_analysis")
)
data_dir <- file.path(bioheart_root, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

cat("BIOHEART_ROOT: ", bioheart_root, "\n", sep = "")
cat("Data dir:      ", data_dir, "\n\n", sep = "")

sentinel <- file.path(data_dir, SENTINEL_SUBDIR)
if (dir.exists(sentinel) && !force) {
  cat("Data already present at ", sentinel, ".\n",
      "Pass --force to re-download.\n", sep = "")
  quit(status = 0)
}

zip_path <- file.path(data_dir, ZIP_NAME)

md5 <- function(path) {
  if (!requireNamespace("tools", quietly = TRUE)) return(NA_character_)
  unname(tools::md5sum(path))
}

needs_download <- TRUE
if (file.exists(zip_path) && !force) {
  cat("Found existing ", ZIP_NAME, ", verifying md5...\n", sep = "")
  if (identical(md5(zip_path), ZIP_MD5)) {
    cat("  md5 matches; skipping download.\n")
    needs_download <- FALSE
  } else {
    cat("  md5 mismatch; re-downloading.\n")
  }
}

if (needs_download) {
  cat("Downloading ~", ZIP_SIZE_GB, " GB from Zenodo...\n", sep = "")
  cat("  ", ZENODO_FILE_URL, " -> ", zip_path, "\n", sep = "")
  ok <- tryCatch({
    utils::download.file(ZENODO_FILE_URL, zip_path,
                         mode = "wb", quiet = FALSE)
    0L
  }, error = function(e) {
    cat("  R download.file failed: ", conditionMessage(e), "\n", sep = "")
    1L
  })
  if (ok != 0L || !file.exists(zip_path)) {
    if (nzchar(Sys.which("curl"))) {
      cat("  Falling back to curl...\n")
      ok <- system2("curl",
                    c("-L", "--fail", "-o", shQuote(zip_path),
                      shQuote(ZENODO_FILE_URL)))
    }
    if (ok != 0L) stop("Download failed. Network or Zenodo issue.")
  }

  cat("Verifying md5...\n")
  got <- md5(zip_path)
  if (!identical(got, ZIP_MD5)) {
    stop("md5 mismatch.\n  expected: ", ZIP_MD5,
         "\n  got:      ", got,
         "\nThe download may be corrupted; delete and re-run.")
  }
  cat("  md5 OK (", ZIP_MD5, ").\n", sep = "")
}

cat("Unzipping into ", data_dir, "...\n", sep = "")
utils::unzip(zip_path, exdir = data_dir)

if (!dir.exists(sentinel)) {
  stop("Unzip completed but expected directory ", sentinel,
       " not found. Archive layout may have changed.")
}

cat("\nDone. Data available under ", sentinel, ".\n", sep = "")
cat("To free disk space you can now: rm ", zip_path, "\n", sep = "")
