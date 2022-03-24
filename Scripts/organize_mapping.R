# Organize mapping depth text files (Mapping/{target}_mapping.txt) for mapping to sequencing from database or standards
library(Biostrings)
library(stringr)
library(fuzzyjoin)

organize_mapping <- function(file, nucleic_acids_file) {
  mapping = read.table(file, header = TRUE, sep = "\t")
  colnames(mapping) <- c("ID", "position", "read_depth")
  mapping$nucleic_acid <- "N"
  
  fasta2dataframe=function(nucleic_acids_file, mapping){
    s = readDNAStringSet(nucleic_acids_file)
    seq_names = names(s)
    sequence = paste(s)
    seq <- cbind.data.frame(seq_names, sequence)
    colnames(seq) <- c("init_ID", "nucleic_acid")
    IDs <- data.frame(unique(mapping$ID))
    colnames(IDs) <- "ID"
    seq <- seq %>% fuzzy_inner_join(IDs, by = c("init_ID" = "ID"), match_fun = str_detect)
    for (i in 1:nrow(seq)){
      seq_sep = data.frame(strsplit(seq$nucleic_acid[i], ""))
      seq_sep <- cbind.data.frame(seq$ID[i], 1:nrow(seq_sep), seq_sep)
      colnames(seq_sep) <- c("ID", "position", "nucleic_acid")
      if (nrow(subset(mapping, ID == seq_sep$ID[1])) != nrow(seq_sep)){
        seq_sep <- seq_sep[1:nrow(subset(mapping, ID == seq_sep$ID[1])),]
        mapping$nucleic_acid[mapping$ID == seq_sep$ID[1]] <- 
          seq_sep$nucleic_acid
      } else {
        mapping$nucleic_acid[mapping$ID == seq_sep$ID[1]] <- 
          seq_sep$nucleic_acid
      }
    }
    
    return(mapping)
  }
  
  mapping <- fasta2dataframe(nucleic_acids_file, mapping)
  
  mapping$nucleic_acid[mapping$nucleic_acid == "a"] <- "A"
  mapping$nucleic_acid[mapping$nucleic_acid == "t"] <- "T"
  mapping$nucleic_acid[mapping$nucleic_acid == "g"] <- "G"
  mapping$nucleic_acid[mapping$nucleic_acid == "c"] <- "C"
  
  return(mapping)
}

#run function
mapping <- organize_mapping(snakemake@input[[1]], snakemake@params[[1]])

write.table(mapping, snakemake@output[[1]], append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

