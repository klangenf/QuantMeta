#!/usr/bin/env Rscript
# Convert R RDS regression models to Python pickle format

library(tidyverse)

# Function to convert R lm object to Python-compatible format
convert_lm_to_python <- function(model, bin_name) {
  # Extract coefficients
  coeffs <- coef(model)

  # Determine response type based on bin
  response_type <- case_when(
    bin_name == "gte1000" ~ "linear",
    bin_name == "100to1000" ~ "log",
    bin_name == "10to100" ~ "log",
    bin_name == "0to10" ~ "log",
    TRUE ~ "linear"
  )

  # Create Python-compatible structure
  python_model <- list(
    coefficients = coeffs,
    response_type = response_type,
    bin = bin_name,
    # Note: X_test, y_original, residuals, rmse will be computed in Python
    X_test = NULL,
    y_original = NULL,
    residuals = NULL,
    rmse = NULL
  )

  return(python_model)
}

# File paths
base_path <- "/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/Regressions/read_depth_variability"
output_path <- "/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/REVAMP/Regressions/read_depth_variability"

# Create output directory if it doesn't exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Convert each model
models_to_convert <- c(
  "Langenfeld_2024_0to10readsperbp" = "bin1",
  "Langenfeld_2024_10to100readsperbp" = "bin2",
  "Langenfeld_2024_100to1000readsperbp" = "bin3",
  "Langenfeld_2024_gte1000readsperbp" = "bin4"
)

for (rds_file in names(models_to_convert)) {
  bin_name <- models_to_convert[rds_file]

  # Load RDS file
  model_path <- file.path(base_path, rds_file)
  if (file.exists(model_path)) {
    cat("Loading", rds_file, "\n")
    model <- readRDS(model_path)

    # Convert to Python format
    python_model <- convert_lm_to_python(model, bin_name)

    # Save as RDS first (temporary), then we'll convert to pickle in Python
    temp_rds_path <- file.path(output_path, paste0(bin_name, "_temp.rds"))
    saveRDS(python_model, temp_rds_path)

    cat("Converted", rds_file, "to", temp_rds_path, "\n")
  } else {
    cat("File not found:", model_path, "\n")
  }
}

cat("R conversion complete. Now run the Python script to create pickle files.\n")