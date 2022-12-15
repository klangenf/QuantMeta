##### Title: "Quant_Unknowns"
##### by: Kathryn Langenfeld
##### date: 10/28/2021

### NOTICE ####
# If the number of Regressions/threshold_read_depth_variability regressions is altered to have different read depth bins
# or have more regressions, the quant_correction_v5,quant_correction_contigs, and quant_correction functions will need 
# to be updated to account for the differences (these differences will be based on decisions made while performing the 
# quant_correct_regression_builder.Rmd)

# Inputs: 
# From your working directory:
#   + database_lengths: file with the lengths of every sequence mapped to in the dataset/fasta file (see Map_Indexes/README), for db: Map_Indexes/{sample}/{target}_length.txt
#   + mapping_results: file from Scripts/organize_mapping_db.R or organize_mapping_contigs.R, "Mapping/{sample}/{target}_mapping.txt" 
#   + path_to_quant_regression: regression output file from quantmeta.R, "Regressions/{sample}/quant_lr
#   + sliding_window_file: file from Scripts/GC_read_depth_var_db.R or GC_read_depth_var_contigs.R, "Mapping/{sample}/{target}_window49_GC_readdepth.txt"
#   + bin_assignment: ONLY FOR CONTIGS (if analyzing results of mapping to a database, insert NA), file listing "contig_ID"	"ID" where ID is the bin name for the respective contig
# sample_name: Sample Name
# descript: Description of mapping targets for filenames
# target_type: insert "database" or "contigs"
# DNA_input: Mass of DNA used in library prep (ng DNA), located in Spike-ins/sample_info.txt (update example for your samples)
# DNA_conc: concentration of DNA in the DNA extract (ng DNA/µL DNA extract), located in Spike-ins/sample_info.txt (update example for your samples)
# R: recovery of microbes in metagenome (if not determined, use value of 1), located in Sample_Characteristics/sample_extraction_info.txt (update example for your samples)
# CF: concentration factor from sample collection to DNA extraction (if not determined, use value of 1), located in Sample_Characteristics/sample_extraction_info.txt (update example for your samples)
# 
# order of inputs: sample_name, descript, database_lengths, path_to_quant_regression, mapping_results, sliding_window_file, bin_assignment, target_type,  DNA_input, DNA_conc, R, CF

##### Required R Packages #############################################
#install.packages("dplyr")
#install.packages("MASS")
#install.packages("scales")
library(dplyr)
library(MASS)
library(scales)

##### Determine coverage, average read depth (i.e. gene copies), and detection threshold parameters of the standards ##############################
detection_threshold <- function(mapping_results, length) {
  E_detect <- readRDS(snakemake@params[[4]])
  
  # make a list of the standards ID in the mapping file
  mapping_results <- as.data.frame(read.table(mapping_results, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  targets = data.frame(unique(mapping_results$ID))
  colnames(targets) <- "ID"
  targets <- merge.data.frame(targets, length, by = "ID")
  targets$B_G <- 0
  targets$gene_copies <- 0
  targets$I_G <- 0
  targets$E_rel <- 0
  
  for(i in 1:nrow(targets)) {
    # select genome information on specific genome from mapping table
    target = subset(mapping_results, ID == targets$ID[i])
    target$I_G_x <- 0
    
    # Calculate total number of bases mapping to genome and the average read depth (gene copies)
    targets$B_G[i] = sum(target$read_depth)
    targets$gene_copies[i] = mean(target$read_depth)
    
    # Calculate I_G and number of positions covered by at least 1 read
    target$I_G_x <- (target$read_depth/targets$B_G[i])*log(target$read_depth/targets$B_G[i])
    targets$I_G[i] <- -sum(na.omit(target$I_G_x))
    
    # Calculate E_rel
    targets$E_rel[i] = targets$I_G[i]/log(targets$length[i])
    
    # Calculate E_detect
    targets$E_detect <- predict(E_detect, targets)
  }
  
  targets$detection_status <- "not_detected"
  targets$detection_status[targets$E_rel >= targets$E_detect] <- "detected"
  
  mapping_targets_analysis = cbind.data.frame("ID" = targets$ID, "E_rel" = targets$E_rel, "gene_copies" = targets$gene_copies, "E_detect" = targets$E_detect, "detection_status" = targets$detection_status)
  
  write.table(mapping_targets_analysis, file = snakemake@output[[3]], append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(mapping_targets_analysis)
} 

##### Non-specific mapping detection and correction #############################################################
# If the number of Regressions/threshold_read_depth_variability regressions is altered to have different read depth bins
# or have more regressions, this function will need to be updated to account for the differences (these differences will 
# be based on decisions made while performing the quant_correct_regression_builder.Rmd)
quant_correction_v5 <- function(input_info, sliding_window, quad_reg1, quad_reg2, quad_reg3, quad_reg4, 
                                cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4) {
  ### DETERMINE TARGETS WITH NON-SPECIFIC MAPPING
  # Determine the maximum RMSE for a particular average read depth across a whole target
  input_info$RMSE_limit <- 0
  ### Targets with average read depth greater than 1,000 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 1000] <- 
    predict(cutoff_function4, newdata = subset(input_info, total_avg_depth >= 1000))
  ### Targets with average read depth between 100 and 1,000 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 100 & input_info$total_avg_depth < 1000] <-
    predict(cutoff_function3, newdata = subset(input_info, total_avg_depth >= 100 & total_avg_depth < 1000))
  ### Targets with average read depth between 10 and 100 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 10 & input_info$total_avg_depth < 100] <-
    predict(cutoff_function2, newdata = subset(input_info, total_avg_depth >= 10 & total_avg_depth < 100))
  ### Targets with average read depth less than 10 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth < 10] <-
    predict(cutoff_function1, newdata = subset(input_info, total_avg_depth < 10))
  # Predict returns log(RMSE_limit) and the RMSE_limit has to be calculated
  input_info$RMSE_limit <- exp(input_info$RMSE_limit)
  # Number of targets detected with non-specific mapping that require correction
  print(nrow(subset(input_info, RMSE > RMSE_limit)))
  # Save the initial RMSE to compare to the corrected RMSE
  input_info$prev_RMSE <- input_info$RMSE
  # Establish new columns to assess the fraction corrected
  input_info$frac_corrected <- 0
  input_info$correction_status <- "correctable"
  
  # List the targets requiring correction
  ID_list <- data.frame(unique(subset(input_info, RMSE > RMSE_limit)))
  # Only undergo correction if there are targets with non-specific mapping
  if (nrow(ID_list) > 0) {
    # Cycle through each target requiring correction 
    for (i in 1:nrow(ID_list)) {
      print("loop number:")
      print(i)
      # Pull out just the relevant sliding window information for the focus target
      wind_clean <- subset(sliding_window, ID == ID_list$ID[i])
      # save the initial read depths along the sliding windows
      wind_clean$initial_avg_depth <- wind_clean$avg_depth
      c = 1
      
      # Correct until the RMSE is less than the allowable RMSE or 20 cycles are completed 
      # (also leave while loop if the total_avg_depth becomes zero)
      while ((input_info$RMSE[input_info$ID == ID_list$ID[i]] > input_info$RMSE_limit[input_info$ID == ID_list$ID[i]]) & (c <= 20) & 
             (input_info$total_avg_depth[input_info$ID == ID_list$ID[i]] > 0)) {
        
        ### IDENTIFY REGIONS OF NON-SPECIFIC MAPPING ALONG TARGET SEQUENCE
        # For each range of average read depths, calculate the predicted read depths and std dev of predictions with 
        # respective regression. Then develop upper and lower thresholds (above thresholds indicate outlier points where
        # non-specific mapping occurs). Currently setting upper and lower thresholds as 1.5 std deviations of observed
        # read depths around predicted read depths --> not setting it around mean observed read depths because there should
        # be some fluctuation of the allowable regions and not using std deviations of predictions because those are not
        # specific to this particular target (may be unnecessarily large or small depending on situation)
        if (wind_clean$total_avg_depth[1] < 10) {
          predict.list <- predict(quad_reg1, newdata = wind_clean, se.fit = TRUE)
          wind_clean <- cbind.data.frame(wind_clean, "pred_depth" = predict.list[[1]], "sd_pred_depth" = predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]]))
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
          predict.list <- predict(quad_reg2, newdata = wind_clean, se.fit = TRUE)
          wind_clean <- cbind.data.frame(wind_clean, "pred_depth" = predict.list[[1]], "sd_pred_depth" = predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]]))
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
          predict.list <- predict(quad_reg3, newdata = wind_clean, se.fit = TRUE)
          wind_clean <- cbind.data.frame(wind_clean, "pred_depth" = predict.list[[1]], "sd_pred_depth" = predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]]))
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) 
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 1000) {
          predict.list <- predict(quad_reg4, newdata = wind_clean, se.fit = TRUE)
          wind_clean <- cbind.data.frame(wind_clean, "pred_depth" = predict.list[[1]], "sd_pred_depth" = predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]]))
          wind_clean$pred_depth <- exp(wind_clean$pred_depth)
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth)
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        # Ensure there is not a "less than zero" mapping option as that is impossible
        wind_clean$lower[wind_clean$lower < 0] <- 0
        wind_clean$upper[wind_clean$upper < 0] <- 0
        # Create variable to easily separate non-specific mapping regions from regions not requiring corrections
        wind_clean$nonspec_region <- "NO"
        ### non-specific mapping regions are those with a read depth greater or less than the previously established cut-offs
        ### do not "invent" reads, therefore, never change areas with a read depth of 0
        wind_clean$nonspec_region[(wind_clean$avg_depth > wind_clean$upper | wind_clean$avg_depth < wind_clean$lower) & wind_clean$avg_depth != 0] <- "YES"
        # Calculate the average read depth across a target without the non-specific mapping regions
        wind_clean$total_avg_depth <- mean(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
        # Occasionally, there will be such low abundance that it barely met detection thresholds and the predicted depth and upper
        # and lower thresholds are all determined to be 0 (originally less than 0), in these instances the total_avg_depth is
        # reset to 0 and the target cannot be quantified
        if (nrow(subset(wind_clean, nonspec_region == "NO")) == 0) {
          wind_clean$total_avg_depth <- 0
        }
        
        # DETERMINE CORRECTED READ DEPTHS FOR NON-SPECIFIC MAPPING REGIONS
        # With the updated average read depth that does not include non-specific mapping regions, recalculate upper and lower
        # limits as was done previously
        if (wind_clean$total_avg_depth[1] < 10 & wind_clean$total_avg_depth[1] > 0) {
          predict.list <- predict(quad_reg1, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
        }
        if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
          predict.list <- predict(quad_reg2, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
        }
        if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
          predict.list <- predict(quad_reg3, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) 
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
        }
        if (wind_clean$total_avg_depth[1] >= 1000) {
          predict.list <- predict(quad_reg4, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth)
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth)
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
        }
        
        # Again, reset negative values to 0 
        wind_clean$lower[wind_clean$lower < 0] <- 0
        wind_clean$upper[wind_clean$upper < 0] <- 0
        # If this isn't corrected, the read mapping end up with "cat ears"
        wind_clean$nonspec_region[wind_clean$avg_depth > wind_clean$upper | wind_clean$avg_depth < wind_clean$lower] <- "YES"
        # Set the avg_depth of non-specific mapping regions to the upper or lower limits if it is too high or too low, respectively
        wind_clean$avg_depth[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth > wind_clean$upper] <- 
          wind_clean$upper[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth > wind_clean$upper]
        wind_clean$avg_depth[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth < wind_clean$lower] <- 
          wind_clean$lower[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth < wind_clean$lower]
        # Remove unnecessary columns from wind_clean
        wind_clean <- cbind.data.frame("ID" = wind_clean$ID, "avg_GC" = wind_clean$avg_GC, "avg_depth" = wind_clean$avg_depth, "E_rel" = wind_clean$E_rel, "total_avg_depth" = wind_clean$total_avg_depth, "E_detect" = wind_clean$E_detect, "detection_status" = wind_clean$detection_status,
                                       "length" = wind_clean$length, "initial_avg_depth" = wind_clean$initial_avg_depth)
        # Recalculate average read depth based on corrected values and save it to input_info data frame
        wind_clean$total_avg_depth <- mean(wind_clean$avg_depth)
        input_info$total_avg_depth[input_info$ID == ID_list[i,1]] <- wind_clean$total_avg_depth[1]
        
        # Update the RMSE and RMSE_limit based on the new total_avg_depth
        if (wind_clean$total_avg_depth[1] < 10) {
          input_info$RMSE[input_info$ID == ID_list[i,1]] <- 
            sqrt(sum((predict(quad_reg1, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <-
            predict(cutoff_function1, newdata = subset(input_info, ID == ID_list[i,1]))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list[i,1]])
        }
        if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
          input_info$RMSE[input_info$ID == ID_list[i,1]] <- 
            sqrt(sum((predict(quad_reg2, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <-
            predict(cutoff_function2, newdata = subset(input_info, ID == ID_list[i,1]))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list[i,1]])
        }
        if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
          input_info$RMSE[input_info$ID == ID_list[i,1]] <- 
            sqrt(sum((predict(quad_reg3, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <-
            predict(cutoff_function3, newdata = subset(input_info, ID == ID_list[i,1]))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list[i,1]])
        }
        if (wind_clean$total_avg_depth[1] >= 1000) {
          input_info$RMSE[input_info$ID == ID_list[i,1]] <- 
            sqrt(sum((predict(quad_reg4, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <- 
            predict(cutoff_function4, newdata = subset(input_info, ID == ID_list[i,1]))
          input_info$RMSE_limit[input_info$ID == ID_list[i,1]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list[i,1]])
        }
        
        # Update loop counter
        print(c)
        c = c + 1
      }
      # Save the initial and corrected mapping results for later comparisons
      if (i == 1) {
        corrected_mapping <- wind_clean
      } else {
        corrected_mapping <- rbind.data.frame(corrected_mapping, wind_clean)
      }
      
      # calculate the fraction of sliding windows were corrected
      input_info$frac_corrected[input_info$ID == ID_list[i,1]] <- nrow(subset(wind_clean, avg_depth != initial_avg_depth))/nrow(wind_clean)
      input_info$correction_status[input_info$ID == ID_list[i,1] & input_info$frac_corrected > 0.2] <- "error"
    }
    
    
    # Write the data frame of initial and corrected mapping to a .txt file
    write.table(corrected_mapping, snakemake@output[[4]], sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    corrected_mapping <- "No targets required correction!"
    write.table(corrected_mapping, snakemake@output[[4]], row.names = FALSE, col.names = FALSE)
  }
  
  return(input_info)
}

quant_correction_contigs <- function(input_info, sliding_window, bin_assignment, quad_reg1, quad_reg2, quad_reg3, quad_reg4, 
                                     cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4) {
  ### DETERMINE TARGETS WITH NON-SPECIFIC MAPPING
  # Determine the maximum RMSE for a particular average read depth across a whole target
  input_info$RMSE_limit <- 0
  ### Targets with average read depth greater than 1,000 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 1000] <- 
    predict(cutoff_function4, newdata = subset(input_info, total_avg_depth >= 1000))
  ### Targets with average read depth between 100 and 1,000 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 100 & input_info$total_avg_depth < 1000] <-
    predict(cutoff_function3, newdata = subset(input_info, total_avg_depth >= 100 & total_avg_depth < 1000))
  ### Targets with average read depth between 10 and 100 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth >= 10 & input_info$total_avg_depth < 100] <-
    predict(cutoff_function2, newdata = subset(input_info, total_avg_depth >= 10 & total_avg_depth < 100))
  ### Targets with average read depth less than 10 reads/bp
  input_info$RMSE_limit[input_info$total_avg_depth < 10] <-
    predict(cutoff_function1, newdata = subset(input_info, total_avg_depth < 10))
  # Predict returns log(RMSE_limit) and the RMSE_limit has to be calculated
  input_info$RMSE_limit <- exp(input_info$RMSE_limit)
  # Save the initial RMSE to compare to the corrected RMSE
  input_info$prev_RMSE <- input_info$RMSE
  # Establish new columns to assess the fraction corrected
  input_info$frac_corrected <- 0
  input_info$correction_status <- NA
  input_info$correction_status[input_info$RMSE > input_info$RMSE_limit] <- "correctable"
  
  
  # List the targets requiring correction
  ID_list <- subset(input_info, RMSE > RMSE_limit)
  ID_list <- data.frame(unique(ID_list$ID))
  colnames(ID_list) <- c("ID")
  # Number of targets detected with non-specific mapping that require correction
  print(nrow(ID_list))
  # Only undergo correction if there are targets with non-specific mapping
  if (nrow(ID_list) > 0) {
    # Cycle through each target requiring correction 
    for (i in 1:nrow(ID_list)) {
      print("loop number:")
      print(i)
      # Pull out just the relevant sliding window information for the focus target
      wind_clean <- subset(sliding_window, ID == ID_list$ID[i])
      # save the initial read depths along the sliding windows
      wind_clean$initial_avg_depth <- wind_clean$avg_depth
      c = 1
      
      # Correct until the RMSE is less than the allowable RMSE or 20 cycles are completed 
      # (also leave while loop if the total_avg_depth becomes zero)
      while ((mean(input_info$RMSE[input_info$ID == ID_list$ID[i]]) > mean(input_info$RMSE_limit[input_info$ID == ID_list$ID[i]])) & (c <= 20) & 
             (mean(input_info$total_avg_depth[input_info$ID == ID_list$ID[i]]) > 0)) {
        
        ### IDENTIFY REGIONS OF NON-SPECIFIC MAPPING ALONG TARGET SEQUENCE
        # For each range of average read depths, calculate the predicted read depths and std dev of predictions with 
        # respective regression. Then develop upper and lower thresholds (above thresholds indicate outlier points where
        # non-specific mapping occurs). Currently setting upper and lower thresholds as 1.5 std deviations of observed
        # read depths around predicted read depths --> not setting it around mean observed read depths because there should
        # be some fluctuation of the allowable regions and not using std deviations of predictions because those are not
        # specific to this particular target (may be unnecessarily large or small depending on situation)
        if (wind_clean$total_avg_depth[1] < 10) {
          predict.list <- predict(quad_reg1, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
          predict.list <- predict(quad_reg2, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
          predict.list <- predict(quad_reg3, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) 
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        if (wind_clean$total_avg_depth[1] >= 1000) {
          predict.list <- predict(quad_reg4, newdata = wind_clean, se.fit = TRUE)
          wind_clean$pred_depth <- predict.list[[1]]
          wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
          wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
          wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth)
          wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth)
          wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth)
        }
        # Ensure there is not a "less than zero" mapping option as that is impossible
        wind_clean$lower[wind_clean$lower < 0] <- 0
        wind_clean$upper[wind_clean$upper < 0] <- 0
        # Create variable to easily separate non-specific mapping regions from regions not requiring corrections
        wind_clean$nonspec_region <- "NO"
        ### non-specific mapping regions are those with a read depth greater or less than the previously established cut-offs
        ### do not "invent" reads, therefore, never change areas with a read depth of 0
        wind_clean$nonspec_region[(wind_clean$avg_depth > wind_clean$upper | wind_clean$avg_depth < wind_clean$lower) & wind_clean$avg_depth != 0] <- "YES"
        # Occasionally, there will be such low abundance that it barely met detection thresholds and the predicted depth and upper
        # and lower thresholds are all determined to be 0 (originally less than 0), in these instances the total_avg_depth is
        # reset to 0 and the target cannot be quantified
        if (nrow(subset(wind_clean, nonspec_region == "NO")) == 0 | nrow(subset(wind_clean, nonspec_region == "NO")) == nrow(wind_clean)) {
          c = 21
          print("No convergence to a solution, evaluate contig(s) quality")
          input_info$correction_status[input_info$ID == ID_list[i,1]] <- "error_no convergence"
        } else {
          # Calculate the average read depth across a target without the non-specific mapping regions
          wind_clean$total_avg_depth <- mean(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          
          # DETERMINE CORRECTED READ DEPTHS FOR NON-SPECIFIC MAPPING REGIONS
          # With the updated average read depth that does not include non-specific mapping regions, recalculate upper and lower
          # limits as was done previously
          if (wind_clean$total_avg_depth[1] < 10 & wind_clean$total_avg_depth[1] > 0) {
            predict.list <- predict(quad_reg1, newdata = wind_clean, se.fit = TRUE)
            wind_clean$pred_depth <- predict.list[[1]]
            wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
            wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
            wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
            wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
            wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          }
          if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
            predict.list <- predict(quad_reg2, newdata = wind_clean, se.fit = TRUE)
            wind_clean$pred_depth <- predict.list[[1]]
            wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
            wind_clean$pred_depth <- exp(wind_clean$pred_depth) - 1
            wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) - 1
            wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
            wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          }
          if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
            predict.list <- predict(quad_reg3, newdata = wind_clean, se.fit = TRUE)
            wind_clean$pred_depth <- predict.list[[1]]
            wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
            wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
            wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth) 
            wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
            wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          }
          if (wind_clean$total_avg_depth[1] >= 1000) {
            predict.list <- predict(quad_reg4, newdata = wind_clean, se.fit = TRUE)
            wind_clean$pred_depth <- predict.list[[1]]
            wind_clean$sd_pred_depth <- predict.list[[2]]*sqrt(predict.list[[3]]+predict.list[[4]])
            wind_clean$pred_depth <- exp(wind_clean$pred_depth) 
            wind_clean$sd_pred_depth <- exp(wind_clean$sd_pred_depth)
            wind_clean$upper <- wind_clean$pred_depth + 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
            wind_clean$lower <- wind_clean$pred_depth - 1.5*sd(wind_clean$avg_depth[wind_clean$nonspec_region == "NO"])
          }
          
          # Again, reset negative values to 0 
          wind_clean$lower[wind_clean$lower < 0] <- 0
          wind_clean$upper[wind_clean$upper < 0] <- 0
          # If this isn't corrected, the read mapping end up with "cat ears"
          wind_clean$nonspec_region[wind_clean$avg_depth > wind_clean$upper | wind_clean$avg_depth < wind_clean$lower] <- "YES"
          # Set the avg_depth of non-specific mapping regions to the upper or lower limits if it is too high or too low, respectively
          wind_clean$avg_depth[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth > wind_clean$upper] <- 
            wind_clean$upper[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth > wind_clean$upper]
          wind_clean$avg_depth[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth < wind_clean$lower] <- 
            wind_clean$lower[wind_clean$nonspec_region == "YES" & wind_clean$avg_depth < wind_clean$lower]
          
          # Recalculate average read depth based on corrected values and save it to input_info data frame
          wind_clean$total_avg_depth <- mean(wind_clean$avg_depth)
          input_info$total_avg_depth[input_info$ID == ID_list$ID[i]] <- wind_clean$total_avg_depth[1]
          
          # Update the RMSE and RMSE_limit based on the new total_avg_depth
          if (wind_clean$total_avg_depth[1] < 10) {
            input_info$RMSE[input_info$ID == ID_list$ID[i]] <- sqrt(sum((predict(quad_reg1, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <-
              predict(cutoff_function1, newdata = subset(input_info, ID == ID_list$ID[i]))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list$ID[i]])
          }
          if (wind_clean$total_avg_depth[1] >= 10 & wind_clean$total_avg_depth[1] < 100) {
            input_info$RMSE[input_info$ID == ID_list$ID[i]] <- 
              sqrt(sum((predict(quad_reg2, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <-
              predict(cutoff_function2, newdata = subset(input_info, ID == ID_list$ID[i]))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list$ID[i]])
          }
          if (wind_clean$total_avg_depth[1] >= 100 & wind_clean$total_avg_depth[1] < 1000) {
            input_info$RMSE[input_info$ID == ID_list$ID[i]] <- 
              sqrt(sum((predict(quad_reg3, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <-
              predict(cutoff_function3, newdata = subset(input_info, ID == ID_list$ID[i]))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list$ID[i]])
          }
          if (wind_clean$total_avg_depth[1] >= 1000) {
            input_info$RMSE[input_info$ID == ID_list$ID[i]] <- 
              sqrt(sum((predict(quad_reg4, newdata = wind_clean) - wind_clean$avg_depth)^2/length(wind_clean$avg_depth)))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <- 
              predict(cutoff_function4, newdata = subset(input_info, ID == ID_list$ID[i]))
            input_info$RMSE_limit[input_info$ID == ID_list$ID[i]] <- exp(input_info$RMSE_limit[input_info$ID == ID_list$ID[i]])
          }
          # Update loop counter
          print(c)
          c = c + 1
        }
        
      }
      # Remove unnecessary columns from wind_clean
      if (is.na(bin_assignment)) {
        wind_clean <- cbind.data.frame("ID" = wind_clean$ID, "avg_GC" = wind_clean$avg_GC, "avg_depth" = wind_clean$avg_depth, "E_rel" = wind_clean$E_rel, "total_avg_depth" = wind_clean$total_avg_depth,
                                     "E_detect" = wind_clean$E_detect, "detection_status" = wind_clean$detection_status, "length" = wind_clean$length, "initial_avg_depth" = wind_clean$initial_avg_depth)
      } else {
        wind_clean <- cbind.data.frame("contig_ID" = wind_clean$contig_ID, "ID" = wind_clean$ID, "avg_GC" = wind_clean$avg_GC, "avg_depth" = wind_clean$avg_depth, "E_rel" = wind_clean$E_rel, "total_avg_depth" = wind_clean$total_avg_depth,
                                       "E_detect" = wind_clean$E_detect, "detection_status" = wind_clean$detection_status, "length" = wind_clean$length, "initial_avg_depth" = wind_clean$initial_avg_depth)
      }
      # Save the initial and corrected mapping results for later comparisons
      if (i == 1) {
        corrected_mapping <- wind_clean
      } else {
        corrected_mapping <- rbind.data.frame(corrected_mapping, wind_clean)
      }
      
      # calculate the fraction of sliding windows were corrected
      input_info$frac_corrected[input_info$ID == ID_list[i,1]] <- nrow(subset(wind_clean, avg_depth != initial_avg_depth))/nrow(wind_clean)
      input_info$correction_status[input_info$ID == ID_list[i,1] & input_info$frac_corrected > 0.2] <- "error_high fraction corrected"
    }
    
    # Write the data frame of initial and corrected mapping to a .txt file
    write.table(corrected_mapping, snakemake@output[[4]], sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    corrected_mapping <- "No targets required correction!"
    write.table(corrected_mapping, snakemake@output[[4]], row.names = FALSE, col.names = FALSE)
  }
  
  return(input_info)
}

## results = mapping_targets_analysis from the detection threshold function
quant_correction <- function(sample_name, descript, sliding_window_file, lengths, results, bin_assignment, target_type) {
  sliding_window <- as.data.frame(read.table(sliding_window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  
  quad_reg1 <- readRDS(snakemake@params[[5]])
  quad_reg2 <- readRDS(snakemake@params[[6]])
  quad_reg3 <- readRDS(snakemake@params[[7]])
  quad_reg4 <- readRDS(snakemake@params[[8]])
  
  cutoff_function1 <- readRDS(snakemake@params[[9]])
  cutoff_function2 <- readRDS(snakemake@params[[10]])
  cutoff_function3 <- readRDS(snakemake@params[[11]])
  cutoff_function4 <- readRDS(snakemake@params[[12]])
  
  targets_list <- subset(lengths, ID %in% sliding_window$ID)
  
  if (!is.na(bin_assignment)) {
    bins <- data.frame(read.table(bin_assignment, header = TRUE, sep = "\t"))
    colnames(bins) <- c("contig_ID", "ID")
    
    colnames(results) <- c("contig_ID", "E_rel", "total_avg_depth", "E_detect", "detection_status")
    results <- merge(bins, results, by = "contig_ID")
    
    colnames(targets_list) <- c("contig_ID", "length")
    results <- merge(results, targets_list, by  = "contig_ID")
  } else {
    results <- merge(results, targets_list, by  = "ID")
  }
  
  if (!(is.na(bin_assignment))) {
    colnames(sliding_window) <- c("contig_ID", "avg_GC", "avg_depth")
    sliding_window <- merge(sliding_window, results, by = "contig_ID")
    colnames(sliding_window) <- c("contig_ID", "avg_GC", "avg_depth", "ID", "E_rel", "total_avg_depth", 
                                  "E_detect", "detection_status", "length")
    colnames(results) <- c("contig_ID", "ID", "E_rel", "total_avg_depth", "E_detect", "detection_status", "length")
    temp <- aggregate(sliding_window[c("avg_depth")], sliding_window[c("ID")], FUN = function(x) mean(x))
    colnames(temp) <- c("ID", "total_avg_depth")
    sliding_window <- merge(sliding_window, temp, by = "ID")
    sliding_window <- cbind.data.frame(sliding_window[,1:5], "total_avg_depth" = sliding_window$total_avg_depth.y, sliding_window[,7:9])
    results <- merge(results, temp, by = "ID")
    results <- cbind.data.frame(results[,1:3], "total_avg_depth" = results$total_avg_depth.y, results[,5:7])
  } else {
    sliding_window <- merge(sliding_window, results, by = "ID")
    colnames(sliding_window) <- c("ID", "avg_GC", "avg_depth", "E_rel", "total_avg_depth", 
                                "E_detect", "detection_status", "length")
    colnames(results) <- c("ID", "E_rel", "total_avg_depth", "E_detect", "detection_status", "length")
  }
  
  ID_list <- data.frame(unique(results$ID))
  for (j in 1:nrow(ID_list)) {
    test_set <- subset(sliding_window, ID == ID_list[j,1])
    if(test_set$total_avg_depth[1] < 10) {
      pred <- predict.lm(quad_reg1, newdata = test_set)
      pred <- exp(pred) - 1
    } 
    if(test_set$total_avg_depth[1] >= 10 & test_set$total_avg_depth[1] < 100) {
      pred <- predict.lm(quad_reg2, newdata = test_set)
      pred <- exp(pred) - 1
    } 
    if(test_set$total_avg_depth[1] >= 100 & test_set$total_avg_depth[1] < 1000) {
      pred <- predict.lm(quad_reg3, newdata = test_set)
      pred <- exp(pred)
    }
    if(test_set$total_avg_depth[1] >= 1000) {
      pred <- predict.lm(quad_reg4, newdata = test_set)
    }   
    results$RMSE[results$ID == ID_list[j,1]] <- sqrt(sum((pred - test_set$avg_depth)^2/length(test_set$avg_depth)))
  }
  
  if(target_type == "database") {
    results_v2 <- quant_correction_v5(results, sliding_window, quad_reg1, quad_reg2, quad_reg3, quad_reg4, 
                                      cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4)
  } else {
    results_v2 <- quant_correction_contigs(results, sliding_window, bin_assignment, quad_reg1, quad_reg2, quad_reg3, quad_reg4, 
                                      cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4)
  }
  
  
  results_v2$reliability <- "low conf/high error"
  results_v2$reliability[results_v2$RMSE <= results_v2$RMSE_limit] <- "high conf/low error"
  if (!(is.na(bin_assignment))) {
    colnames(results_v2) <- c("ID", "contig_ID", "E_rel", "gene_copies", "E_detect", 
                              "detection_status", "length", "RMSE", "RMSE_limit", "prev_RMSE", 
                              "frac_corrected", "correction_status", "reliability")
  } else {
    colnames(results_v2) <- c("ID", "E_rel", "gene_copies", "E_detect", 
                              "detection_status", "length", "RMSE", "RMSE_limit", "prev_RMSE", 
                              "frac_corrected", "correction_status", "reliability")
  }
  
  write.table(results_v2, snakemake@output[[2]], sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(results_v2)
}

##### Predicted Concentration (gene copies/µL DNA extract) Calculator ########################################3
prediction <- function (mapping, DNA_input, DNA_conc, target_type) {
  pred_conc_calc <- function(gc, mass, conc) {
    predicted <- gc*conc/mass
    # (gene copy standard_x)/(ng DNA library insert) = gene copies standard_x/ng total DNA*(ng DNA/µL DNA extract)
    return(predicted)
  }
  
  if (target_type == "database") {
    result = data.frame(mapping$ID)
    result$predicted_conc <- pred_conc_calc(mapping$gene_copies, DNA_input, DNA_conc)
    colnames(result) <- c("ID", "predicted_conc")
  } else {
    result = cbind.data.frame("ID" = mapping$ID, "contig_ID" = mapping$contig_ID)
    result$predicted_conc <- pred_conc_calc(mapping$gene_copies, DNA_input, DNA_conc)
    colnames(result) <- c("ID", "contig_ID", "predicted_conc")
  }
  
  return(result)
}

##### Calculate concentrations in original sample ################################
in_original <- function(concentration, recovery, conc_factor) {
  concentration$`concentration (gc/µL)` <- concentration$`concentration (gc/µL)`/(recovery*conc_factor)
  concentration$`std deviation (gc/µL)` <- concentration$`std deviation (gc/µL)`/(recovery*conc_factor)
  return(concentration)
}

##### Quantify unknown targets ###################################################
quant_unknown <- function(database_lengths, path_to_quant_regression, mapping_results, sliding_window_file, bin_assignment, target_type, DNA_input, DNA_conc, R, CF) {
  lengths <- data.frame(read.table(database_lengths, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
  colnames(lengths) <- c("ID", "length")
  
  quant_rel <- readRDS(path_to_quant_regression)
  
  # Assess which targets are above the detection threshold
  mapping <- detection_threshold(mapping_results, lengths)
  mapping[,2:4] <- lapply(mapping[,2:4], as.numeric)
  mapping <- subset(mapping, E_rel >= E_detect)
  
  # Detect and correct non-specific mapping
  if (nrow(mapping) != 0) {
    mapping <- quant_correction(sample_name, descript, sliding_window_file, lengths, mapping, bin_assignment, target_type)
    
    # Remove targets that are not quantifiable (correction_status == "error_no convergence" or "error_high fraction corrected)
    mapping <- subset(mapping, !(correction_status %in% c("error_no convergence", "error_high fraction corrected")))
  } else {
    temp <- "No targets above detection!"
    write.table(temp, snakemake@output[[4]], row.names = FALSE, col.names = FALSE)
    write.table(temp, snakemake@output[[2]], row.names = FALSE, col.names = FALSE)
  } 
  
  
  # Convert targets' relative abundances to units of (gene copies/ng DNA)
  results <- prediction(mapping, DNA_input, DNA_conc, target_type)
  
  # Convert relative abundance to absolute abundance
  pred <- predict(quant_rel, newdata = results, se.fit = TRUE)
  if (target_type == "database") {
    results <- cbind.data.frame(results$ID, 10^(pred[[1]]), 10^(pred[[2]]*sqrt(pred[[3]]+1)))
    colnames(results) <- c("ID", "concentration (gc/µL)",  "std deviation (gc/µL)")
  } else {
    results <- cbind.data.frame(results$ID, results$contig_ID, 10^(pred[[1]]), 10^(pred[[2]]*sqrt(pred[[3]]+1)))
    colnames(results) <- c("ID", "contig_ID", "concentration (gc/µL)",  "std deviation (gc/µL)")
  }
  
  # Convert to concentrations in wastewater (gc/µL)
  results <- in_original(results, R, CF)
  
  write.table(results, snakemake@output[[1]], sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(results)
}

### Setup for snakemake
sample_info <- data.frame(read.table(snakemake@input[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE))
char <- data.frame(read.table(snakemake@input[[2]], header = TRUE, sep = "\t"))

if (snakemake@params[[3]] == "database") {
  results <- quant_unknown(snakemake@input[[3]], 
                           snakemake@input[[6]], snakemake@input[[4]], snakemake@input[[5]], NA, snakemake@params[[3]], 
                           sample_info$lib_mass[sample_info$Sample == snakemake@params[[1]]], sample_info$DNA_conc[sample_info$Sample == snakemake@params[[1]]], 
                           char$recovery[char$sample == snakemake@params[[1]]], char$CF[char$sample == snakemake@params[[1]]])
} else {
  results <- quant_unknown(snakemake@input[[3]],
                           snakemake@input[[6]], snakemake@input[[4]], snakemake@input[[5]], snakemake@input[[7]], snakemake@params[[3]],
                           sample_info$lib_mass[sample_info$Sample == snakemake@params[[1]]], sample_info$DNA_conc[sample_info$Sample == snakemake@params[[1]]],
                           char$recovery[char$sample == snakemake@params[[1]]], char$CF[char$sample == snakemake@params[[1]]])
}

