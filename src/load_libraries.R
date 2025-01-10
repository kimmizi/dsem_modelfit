# Function to centralize library loading
load_libraries <- function() {
  # List of required libraries
  required_packages <- c(
    "doParallel", "foreach", "doSNOW", "future", "profvis", "furrr",
    "mvtnorm", "matrixcalc", "lavaan", "blavaan", "rstan"
  )
  
  # Install and load each library
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg, repos = "http://cran.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
  
  # Specific configurations for certain packages
  if ("rstan" %in% installed.packages()) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores()) # Enable parallel processing for rstan
  }
  
  message("All required libraries loaded and configured.")
}

# Call the library-loading function
load_libraries()
