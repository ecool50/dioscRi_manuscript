# Setup Environment for dioscRi
# Run this script once to install the correct versions of keras3 and TensorFlow.
# Requires: R >= 4.5.0, Python 3.10 (via pyenv or system install)
#
# Usage: Rscript scripts/setup_environment.R

cat("=== dioscRi Environment Setup ===\n\n")

# --- 1. Install correct R packages ---
cat("Step 1: Installing R packages...\n")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

# `here` is used by every revision script to resolve revision_root (the dioscRi_manuscript repo root).
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

# Bioconductor packages (treekoR, ggtree, SummarizedExperiment, SingleCellExperiment)
# are pulled in by the dioscRi package. install.packages() alone cannot resolve them,
# so install BiocManager up front and pre-install the Bioc dependencies.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_deps <- c("SummarizedExperiment", "SingleCellExperiment",
               "treekoR", "ggtree", "BiocStyle")
missing_bioc <- bioc_deps[!vapply(bioc_deps, requireNamespace,
                                   logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  cat("  Installing Bioconductor deps:",
      paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# keras3 1.2.0 (must match manuscript)
if (packageVersion("keras3") != "1.2.0") {
  cat("  Installing keras3 1.2.0...\n")
  remotes::install_version("keras3", version = "1.2.0",
                           repos = "https://cloud.r-project.org",
                           upgrade = "never")
} else {
  cat("  keras3 1.2.0 already installed.\n")
}

# tensorflow R package
if (!requireNamespace("tensorflow", quietly = TRUE)) {
  install.packages("tensorflow", repos = "https://cloud.r-project.org")
}

# reticulate
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org")
}

# --- 2. Create Python virtualenv with correct TF version ---
cat("\nStep 2: Setting up Python virtualenv...\n")

library(reticulate)

venv_name <- "r-reticulate"
venv_path <- file.path("~/.virtualenvs", venv_name)

# Find Python 3.10 (preferred) or fall back to available Python 3
python_path <- NULL
candidates <- c(
  Sys.which("python3.10"),
  "/usr/local/bin/python3.10",
  file.path(Sys.getenv("HOME"), ".pyenv/versions/3.10.15/bin/python"),
  file.path(Sys.getenv("HOME"), ".pyenv/shims/python3.10"),
  Sys.which("python3")
)
for (p in candidates) {
  if (nzchar(p) && file.exists(p)) {
    python_path <- p
    break
  }
}

if (is.null(python_path)) {
  stop("Could not find Python 3. Please install Python 3.10 via pyenv or your system package manager.")
}
cat("  Using Python:", python_path, "\n")

# Create virtualenv if it doesn't exist
if (!dir.exists(normalizePath(venv_path, mustWork = FALSE))) {
  cat("  Creating virtualenv:", venv_path, "\n")
  virtualenv_create(venv_name, python = python_path)
} else {
  cat("  Virtualenv already exists:", venv_path, "\n")
}

# Install Python packages with pinned versions
cat("  Installing Python packages...\n")
virtualenv_install(venv_name, packages = c(
  "tensorflow==2.16.2",
  "numpy==1.26.4",
  "pandas",
  "scipy"
), pip_options = "--quiet")

# --- 3. Verify installation ---
cat("\nStep 3: Verifying installation...\n")

Sys.setenv(RETICULATE_PYTHON = file.path(venv_path, "bin/python"))
reticulate::use_virtualenv(venv_name, required = TRUE)

tf <- reticulate::import("tensorflow")
np <- reticulate::import("numpy")

cat("  keras3 R package:", as.character(packageVersion("keras3")), "\n")
cat("  tensorflow R package:", as.character(packageVersion("tensorflow")), "\n")
cat("  TensorFlow Python:", tf$version$VERSION, "\n")
cat("  NumPy:", np$`__version__`, "\n")
cat("  Python:", reticulate::py_config()$python, "\n")

# Validate versions
errors <- c()
if (packageVersion("keras3") != "1.2.0") errors <- c(errors, "keras3 should be 1.2.0")
if (tf$version$VERSION != "2.16.2") errors <- c(errors, "TensorFlow should be 2.16.2")

if (length(errors) > 0) {
  cat("\nWARNING: Version mismatches detected:\n")
  for (e in errors) cat("  -", e, "\n")
  cat("Results may not be reproducible.\n")
} else {
  cat("\nAll versions match manuscript requirements.\n")
}

# --- 4. Quick smoke test ---
cat("\nStep 4: Smoke test...\n")

library(keras3)
input <- layer_input(shape = 10L)
output <- input %>%
  layer_dense(5L, activation = "relu") %>%
  layer_dense(10L, activation = "sigmoid")
model <- keras_model(input, output)
model %>% compile(optimizer = "rmsprop", loss = "binary_crossentropy")

x <- matrix(runif(100), nrow = 10)
model %>% fit(x, x, epochs = 2L, verbose = 0L)
pred <- model %>% predict(x, verbose = 0L)

cat("  Model training and prediction OK.\n")

cat("\n=== Setup complete. ===\n")
cat("To use this environment in scripts, add this at the top:\n")
cat("  Sys.setenv(RETICULATE_PYTHON = '~/.virtualenvs/r-reticulate/bin/python')\n")
cat("  library(reticulate)\n\n")

cat("Next step: download the manuscript data (~7.2 GB) from Zenodo:\n")
cat("  export BIOHEART_ROOT=/path/to/bioheart_analysis\n")
cat("  Rscript revision/scripts/helpers/fetch_data.R\n\n")
