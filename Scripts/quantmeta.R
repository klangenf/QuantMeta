##### Title: "QuantMeta"
##### by: Dr. Kathryn Langenfeld
##### date: 02/15/2022

# Inputs: 
  # From your working directory:
  #   + STD_MIXES: location of "STD_MIXES.txt"
  #   + mapping_output: location of output ".txt" files from mapping and subsequent samtools/mapping steps 
  # sample_name: Sample Name
  # DNA_input: Mass of DNA used in library prep (ng DNA)
  # DNA_conc: concentration of DNA in the DNA extract (ng DNA/µL DNA extract)
  # MIX: Mix of sequins spiked into sample (Hardwick et al. 2019) (refer to as "MIX_A", "MIX_B", "MIX_C")
  # sequins_spike: Sequins spike-in concentration (ng sequins/ng total DNA)
  # PCR_quant_file: Measured spike-in concentrations with qPCR or ddPCR (gc sequins/µL DNA extract) file, if ddPCR/qPCR was not performed, insert NA
  # ssDNA_stds: Were ssDNA standards used? (YES/NO)
  # ssDNA_spike: ssDNA standards spike-in concentration (ng ssDNA standards mix/ng total DNA), if "NO" for ssDNA_stds insert 0

# Regression: 
  # Detection threshold regression should be available in the Regressions/detection directory
  # Regressions/detection/Langenfeld_2022_E_detect can be used unless more or less stringent detection
  # thresholds are warranted. New detection thresholds can be created with confident_detection_regression_builder.R

# Outputs will be saved to Results/{sample_name}
  # make each Results/{sample_name} output folder before running the script

##### Required R Packages #############################################
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("MASS")
#install.packages("scales")
#install.packages("ggpmisc")
library(dplyr)
library(ggplot2)
library(MASS)
library(scales)
library(ggpmisc)

##### Determine coverage, average read depth (i.e. gene copies), and detection threshold parameters of the standards ##############################
detection_threshold <- function(sample_name, mapping_results, target_length, detect_thresh) {
  E_detect <- readRDS(detect_thresh)
  
  # make a list of the standards ID in the mapping file
  mapping_results <- as.data.frame(read.table(mapping_results, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  targets = data.frame(unique(mapping_results$ID))
  colnames(targets) <- "ID"
  targets <- merge.data.frame(targets, target_length, by = "ID")
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
    targets$E_detect[i] <- predict(E_detect, newdata = targets[i,])
  }
  
  targets$detection_status <- "not_detected"
  targets$detection_status[targets$E_rel >= targets$E_detect] <- "detected"
  
  mapping_targets_analysis = cbind.data.frame(targets$ID, targets$E_rel, targets$gene_copies, targets$E_detect, targets$detection_status)
  colnames(mapping_targets_analysis) <- c("ID", "E_rel", "gene_copies", "E_detect", "detection_status")
  
  output = "Mapping/sample/standards_mapping_analysis.txt"
  output = gsub("sample", sample_name, output)
  write.table(mapping_targets_analysis, file = output, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(mapping_targets_analysis)
} 

##### Known Spike-in Concentration (gene copies/µL DNA extract) Calculator ##################################################
expected_conc <- function (MIX, spike, DNA_conc) {
  expected_conc_calc <- function(std_frac, all_std_mass_spike, DNA_conc, mass_std) {
    exp_conc <- std_frac*all_std_mass_spike*DNA_conc/mass_std
    # (fraction by gene copy standard_x/all_standards)*(ng all standards/ng total DNA)*(ng total DNA/µL DNA extract)/(ng standard_x/gene copy standard_x) 
    # = gene copy standard_x/ng total DNA
    return(exp_conc)
  }
  
  result = data.frame(MIX$ID)
  result$known_conc <- expected_conc_calc(MIX$MIX, spike, DNA_conc, MIX$Mass)
  colnames(result) <- c("ID", "known_conc")
  
  return(result)
}

##### Predicted Concentration (gene copies/µL DNA extract) Calculator ########################################3
prediction <- function (mapping, DNA_input, DNA_conc) {
  pred_conc_calc <- function(gc, mass, conc) {
    predicted <- gc*conc/mass
    # (gene copy standard_x)/(µL DNA extract) = gene copies standard_x/ng total DNA*(ng DNA/µL DNA extract)
    return(predicted)
  }
  
  result = data.frame(mapping$ID)
  result$predicted_conc <- pred_conc_calc(mapping$gene_copies, DNA_input, DNA_conc)
  colnames(result) <- c("ID", "predicted_conc")
  
  return(result)
}

##### Visualize the known concentration vs. predicted concentration results ##########################################
visualize <- function(sample_name, results) {
  pretty_plot <- theme_classic() + theme(
  text = element_text(family = "Lucinda Sans", color = "black"),
  plot.margin = margin(2,2,2,2, "cm"),
  axis.line.x.bottom = element_line(color = "black"),
  axis.line.y.left = element_line(color = "black"),
  panel.border = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 15, face = "bold"),
  axis.title = element_text(size = 12, face = "bold"), 
  axis.text.y = element_text(size = 12, color = "#000000"),
  axis.text.x = element_text(size = 12, color = "#000000"))
  
  my.formula <- y ~ x
  ggplot(results, aes(x = predicted_conc, y = known_conc)) + geom_point(shape = 1) +
    pretty_plot + geom_smooth(method = "lm", se = FALSE, color = "black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~"), size = 12),
                 parse = TRUE) +
    labs(x = "Metagenome Gene Copies (gc/µL)", y = "Concentration (gc/µL)", color = "Standard\nType") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
  
  filename = "Results/sample/standards_rel_to_abs.png"
  filename = gsub("sample", sample_name, filename)
  ggsave(filename, plot = last_plot(), scale = 2, dpi = 320, limitsize = TRUE)
}

###### Take output from mapping to determine the predicted concentration for each standard to compare to the known spike-in concentration of each standard ######
quantmeta <- function(sample_name, STD_MIXES, mapping_results, DNA_input, DNA_conc, MIX, sequins_spike, PCR_quant_file, ssDNA_stds, ssDNA_spike, detect_thresh) {
  # Read in STD_MIXES.txt table
  STD_MIXES = as.data.frame(read.table(file = STD_MIXES, sep ="\t", header = TRUE, stringsAsFactors = FALSE))

  # Reduce table to only the relevant sequins spiked into the samples
  keep <- c("ID", "length", "Mass", MIX)
  STD_MIX <- STD_MIXES[keep]
  colnames(STD_MIX) <- c("ID", "length", "Mass", "MIX")
  STD_MIX <- na.omit(STD_MIX)
  
  if (ssDNA_stds == "YES") {
    # make MIX table for ssDNA standards
    keep <- c("ID","length", "Mass", "ssDNA")
    ssDNA_STD_MIX <- STD_MIXES[keep]
    colnames(ssDNA_STD_MIX) <- c("ID", "length", "Mass", "MIX")
    ssDNA_STD_MIX <- na.omit(ssDNA_STD_MIX)
    STD_MIX <- rbind.data.frame(STD_MIX, ssDNA_STD_MIX)
  } else {
    ssDNA_STD_MIX <- data.frame(matrix(ncol=4, nrow=0))
    colnames(ssDNA_STD_MIX) <- c("ID", "length", "Mass", "MIX")
  }
  
  std_length <- cbind.data.frame(STD_MIX$ID, STD_MIX$length)
  colnames(std_length) <- c("ID", "length")
  
  # assess which standards are above the detection threshold
  mapping <- detection_threshold(sample_name, mapping_results, std_length, detect_thresh)
  mapping <- transform(mapping, E_rel = as.numeric(as.character(E_rel)), gene_copies = as.numeric(as.character(gene_copies)), E_detect = as.numeric(as.character(E_detect)))
  mapping <- subset(mapping, E_rel >= E_detect)
  
  # Read in mapping_output, add in 0 for sequins not in the metagenome
  dsDNA_mapping = subset(mapping, !(ID %in% ssDNA_STD_MIX$ID))
  
  # Calculate the known spike-in concentration for each standard (gene copies/ng DNA)
  if (file.size(PCR_quant_file)<500) {
    results <- expected_conc(subset(STD_MIX, !(ID %in% ssDNA_STD_MIX$ID)), sequins_spike, DNA_conc)
  } else {
    # column 1: list of IDs (ID), column 2: list of standard concentrations in the DNA extract (known_conc)
    results <- as.data.frame(read.table(PCR_quant_file, header = TRUE, stringsAsFactors = FALSE))
    results <- cbind.data.frame(results$ID, results[c(sample_name)])
    colnames(results) <- c("ID", "known_conc")
    results <- subset(results, !(ID %in% ssDNA_STD_MIX))
  }
  
  # Calculate the predicted spike-in concentration based on the metagenome for each standard (gene copies/ng DNA)
  predict <- prediction(dsDNA_mapping, DNA_input, DNA_conc)
  results <- merge(results, predict, by = "ID")
  
  # Analyze ssDNA standards if the standards were added into the sample
  if (ssDNA_stds == "YES") {
    # read in ssDNA mapping results
    ssDNA_mapping = subset(mapping, ID %in% ssDNA_STD_MIX$ID)
    
    # Calculate the known spike-in concentration for each standard (gene copies/ng DNA)
    if (is.na(PCR_quant_file)) {
      ssDNA_results <- expected_conc(ssDNA_STD_MIX, ssDNA_spike, DNA_conc)
    } else {
      # column 1: list of IDs (ID), column 2: list of standard concentrations in the DNA extract (known_conc)
      ssDNA_results <- as.data.frame(read.table(PCR_quant_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
      ssDNA_results <- cbind.data.frame(ssDNA_results$ID, ssDNA_results[c(sample_name)])
      colnames(ssDNA_results) <- c("ID", "known_conc")
      ssDNA_results <- subset(ssDNA_results, ID %in% ssDNA_STD_MIX)
    }
    
    # Calculate the predicted spike-in concentration based on the metagenome for each standard (gene copies/ng DNA)
    predict <- prediction(ssDNA_mapping, DNA_input, DNA_conc)
    ssDNA_results <- merge(ssDNA_results, predict, by = "ID")
    
    results <- data.frame(rbind(results, ssDNA_results))
  }
  
  ### Create linear regression relating predicted and known concentrations of standards
  quant_rel <- lm(log10(known_conc) ~ log10(predicted_conc), data = results)
  
  reg_name <- "Regressions/quantification/sample_rel_to_abs"
  reg_name <- gsub("sample", sample_name, reg_name)
  saveRDS(quant_rel, file = reg_name)
  
  visualize(sample_name, results)
  
  return(results)
}

##### Run quantmeta for a sample ############################
sample_info <- as.data.frame(read.table(file = snakemake@input[[2]], sep ='\t', header = TRUE, stringsAsFactors = FALSE))
sample_info <- subset(sample_info, Sample == snakemake@params[[1]])

results <- quantmeta(sample_info$Sample, snakemake@input[[3]], snakemake@input[[1]], sample_info$lib_mass, sample_info$DNA_conc, sample_info$mix, sample_info$dsDNA_spike, snakemake@input[[4]], sample_info$ssDNA, sample_info$ssDNA_spike, snakemake@input[[5]])

write.table(results, snakemake@output[[1]], sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)
