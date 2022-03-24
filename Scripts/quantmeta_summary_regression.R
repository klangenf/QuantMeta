##### Title: "QuantMeta: Abs-Rel Regression for All Samples"
##### by: Dr. Kathryn Langenfeld
##### date: 02/15/2022

# Purpose: Create a single linear regression relating absolute and relative 
# abundances based on the results of spiking in standards to all samples

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

##### Visualize the known concentration vs. predicted concentration results ##########################################
visualize <- function(results) {
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
  
  filename = "Results/standards_rel_to_abs.png"
  ggsave(filename, plot = last_plot(), scale = 2, dpi = 320, limitsize = TRUE)
}

num_samples <- nrow(as.data.frame(read.table(snakemake@params[[1]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)))

for (i in 1:num_samples) {
  temp <- as.data.frame(read.table(snakemake@input[[i]], header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  
  if (i == 1) {
    qm_results <- temp
  } else {
    qm_results <- rbind.data.frame(qm_results, temp)
  }
}

quant_rel <- lm(log10(known_conc) ~ log10(predicted_conc), data = qm_results)

saveRDS(quant_rel, snakemake@output[[1]])

visualize(qm_results)
