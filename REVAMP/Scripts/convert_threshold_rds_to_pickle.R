#!/usr/bin/env Rscript
# Convert R RDS RMSE threshold functions to Python pickle format

# Function to convert R lm object to Python-compatible format
convert_lm_to_python <- function(model, bin_name) {
  # Extract coefficients
  coeffs <- coef(model)

  # Create Python-compatible structure for RMSE limit functions
  python_model <- list(
    coefficients = coeffs,
    type = "log_linear",
    bin = bin_name
  )

  return(python_model)
}

# Function to convert constant value to Python format
convert_constant_to_python <- function(value, bin_name) {
  python_model <- list(
    value = value,
    type = "constant",
    bin = bin_name
  )

  return(python_model)
}

# File paths
base_path <- "/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/Regressions/threshold_read_depth_variability"
output_path <- "/Users/kathrynlangenfeld/Documents/UMich_Postdoc/QuantMeta/REVAMP/Regressions/threshold_read_depth_variability"

# Create output directory if it doesn't exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Convert each threshold function
models_to_convert <- c(
  "Langenfeld_2024_0to10readsperbp" = "func1",
  "Langenfeld_2024_10to100readsperbp" = "func2",
  "Langenfeld_2024_100to1000readsperbp" = "func3",
  "Langenfeld_2024_gte1000readsperbp" = "func4"
)

for (rds_file in names(models_to_convert)) {
  func_name <- models_to_convert[rds_file]

  # Load RDS file
  model_path <- file.path(base_path, rds_file)
  if (file.exists(model_path)) {
    cat("Loading", rds_file, "\n")
    model <- readRDS(model_path)

    # Print model summary for debugging
    cat("Model summary for", func_name, ":\n")
    if (is.numeric(model)) {
      # It's a constant value
      cat("Constant value:", model, "\n")
      python_model <- convert_constant_to_python(model, func_name)
    } else {
      # It's a linear model
      print(summary(model))
      python_model <- convert_lm_to_python(model, func_name)
    }

    # Save as RDS first (temporary), then we'll convert to pickle in Python
    temp_rds_path <- file.path(output_path, paste0(func_name, "_temp.rds"))
    saveRDS(python_model, temp_rds_path)

    cat("Converted", rds_file, "to", temp_rds_path, "\n")
  } else {
    cat("File not found:", model_path, "\n")
  }
}

cat("R conversion complete. Now run the Python script to create pickle files.\n")