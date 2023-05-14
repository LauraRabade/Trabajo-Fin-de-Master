setwd("/home/nodotea/LAURA")
getwd()
library(stringr)
library(tidyverse)

# --------------------------------------------------------------------------
# FUNCTION FOR DMRCATE OUTPUT CLEANING
# --------------------------------------------------------------------------
# Cleans the output of DMRCate analysis

# Arguments: DMR = Data frame containing DMR information


dmrcate_cleaning <- function(DMR){
    
  # DATA CLEANING--------------------------------------------------------------
  # Sort DMRs by meandiff in descending order and select relevant columns
  DMR_genes <- DMR %>%
    arrange(desc(meandiff)) %>%
    select(overlapping.genes, meandiff)
  
  # Filter DMRs with positive differential methylation
  SIGN_POS <- DMR_genes %>%
    filter(meandiff > 0) 
  
  # Filter DMRs with negative differential methylation
  SIGN_NEG <- DMR_genes %>%
    filter(meandiff < 0) 
  
  
  # SIGN_POS--------------------------------------------------------------------
  # Extract overlapping genes without commas
  pos_genes <- SIGN_POS$overlapping.genes
  pos_genes <- pos_genes[!is.na(pos_genes)] # Remove rows containing NA
  temp_columns <- grep(",", pos_genes, invert = TRUE) # Remove rows containing ","
  lista1 <- pos_genes[temp_columns] # Add genes without "," to list
  
  # Replace commas with blank spaces in genes and split into separate genes
  temp_genes <- gsub(",", " ", pos_genes) 
  temp_columns <- grep(" ", temp_genes) 
  temp_genes2 <- temp_genes[temp_columns] 
  temp_genes3 <- strsplit(temp_genes2, " +") 
  lista2 <- unlist(temp_genes3, recursive = FALSE) # New list with all splitted genes
  
  # Join the two lists of genes
  GENES_SIGN_POS <- c(lista1, lista2)
  
  
  # SIGN_NEG--------------------------------------------------------------------
  # Extract overlapping genes without commas
  neg_genes <- SIGN_NEG$overlapping.genes
  neg_genes <- neg_genes[!is.na(neg_genes)] # Remove rows containing NA
  temp_columns <- grep(",", neg_genes, invert = TRUE) # Remove rows containing ","
  lista3 <- neg_genes[temp_columns] # Add genes without "," to list
  
  # Replace commas with blank spaces in genes and split into separate genes
  temp_genes <- gsub(",", " ", neg_genes)
  temp_columns <- grep(" ", temp_genes)
  temp_genes2 <- temp_genes[temp_columns]
  temp_genes3 <- strsplit(temp_genes2, " +")
  lista4 <- unlist(temp_genes3, recursive = FALSE) # New list with all splitted genes
  
  # Join the two lists of genes
  GENES_SIGN_NEG <- c(lista3, lista4)
  
  
  # DATA FILE-------------------------------------------------------------------
  # Write positive genes to file
  pruebapos <- t(GENES_SIGN_POS)
  write.table(pruebapos, file = "genes_pos_dmrcate.txt", sep = " ", quote = FALSE, row.names = "SIGN_POS", col.names = FALSE)
  
  # Write negative genes to file
  pruebaneg <- t(GENES_SIGN_NEG)
  write.table(pruebaneg, file = "genes_neg_dmrcate.txt", sep = " ", quote = FALSE, row.names = "SIGN_NEG", col.names = FALSE)
  
}

######### Obtention of the 2 lists using the function ######### 
dmrcate <- "./ChAMP_Def_analisis/RESULTS_T0_HC_vs_T0_HPRE/DMRs/myDMR_DMRCate.csv"
DMR_dmrcate <- read.table(dmrcate, header = TRUE, sep = ",", dec = ".")
dmrcate_cleaning(DMR_dmrcate)





