# conda_create("r-reticulate")

devtools::install_github("rstudio/reticulate")
library(reticulate)
conda_version()
conda_list()
# conda_python(envname = "r-reticulate")    
conda_install(envname = "r-reticulate",packages="r-reticulate")
conda_install(envname = "r-reticulate",packages="r-tensorflow")
conda_install(envname = "r-reticulate",packages="r-keras")

library(reticulate)
library(tensorflow)
library(keras)
conda_python(envname = "r-reticulate")

# reticulate::virtualenv_create("keras_env", pip_version = 3, python = "3.10")
library(reticulate)
library(tensorflow)
library(keras)
conda_install(envname = "keras_env",packages="r-reticulate")
conda_install(envname = "keras_env",packages="r-tensorflow")
conda_install(envname = "keras_env",packages="r-keras")
# conda_python(envname = "keras_env")
install_tensorflow(method = "conda" ,envname = "keras_env")
install_keras(method = "conda" ,envname = "keras_env", version = "2.16")
k <- backend()


