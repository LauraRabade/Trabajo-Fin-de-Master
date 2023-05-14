setwd("/home/nodotea/LAURA")

###############################################################################
##################### ESTIMATION OF GESTATIONAL DNAm AGE ######################
###############################################################################
# Step repeated for each comparison: HC-0 vs HPRE-N, HC-0 vs HPRE-0, HPRE-0 and HPRE-N

# --- LOAD METHYLATION DATA --- 
library("ChAMP")
myLoad <- champ.load("/home/nodotea/LAURA/ChAMP_Def_analisis/DATOS_21_22_DEF", # location of idat files of each comparison
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     beadCutoff=0.05,
                     ProbeCutoff=0, 
                     filterBeads=TRUE,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     arraytype="EPIC")

# --- QUALITY CONTROL ---
myQC <- champ.QC(beta = myLoad$beta,
                 pheno=myLoad$pd$Sample_Group,
                 Feature.sel="None",
                 resultsDir="None")

# --- NORMALIZATION ---
myNorm <- champ.norm(beta=myLoad$beta,
                     method="BMIQ", # BMIQ Method
                     arraytype="EPIC", 
                     cores=16)

# --- EPIGENETIC CLOCK ---
# Load packages
library(methylclockData)
library(methylclock)
library(tibble)
library(tidyverse)
library(Biobase)
library(impute)
library(minfi)
library(ggpmisc)
library(GEOquery)

# Load data (coef, intercept) of gestational DNA methylation clocks
load_DNAmGA_Clocks_data()

# Check wheter input data contains the required CpGs for the implemented clocks for Gestational Age
cpgs.missing.GA <- checkClocksGA(myNorm)

# Gestational DNAm age estimation (in weeks) using different DNA methylation clocks
age <- child_weeks$Dias.totales # load chronological age
estimatedGA_age <- DNAmGA(myNorm, 
                         toBetas = FALSE,
                         fastImp = FALSE, # imputation of missing CpGs for the clock
                         normalize = FALSE, # normalization of data: already done
                         age,
                         cell.count = FALSE,
                         min.perc = 0.8,) # min 80% of CpGs included in the clock 

# Create csv with gestational age estimation
write.table(x = estimatedGA_age, 
            file = "estGA_file.csv",
            sep = ",", 
            row.names = FALSE, 
            col.names = TRUE)

#------------------------------------------------------------------------------
########## STATISTICS FOR GESTATIONAL AGE ESTIMATED WITH EPIC CLOCK ########## 
#------------------------------------------------------------------------------
library("psych")

# --- STATS FOR T0_HC_TN_HPRE ---
# Read file with data
T0_HC_TN_HPRE <- read.csv("T0_HC_TN_HPRE/estGA_TN_HPRE.csv")

# Divide by sample group
HC <- T0_HC_TN_HPRE[81:125, ]
HPRE <- T0_HC_TN_HPRE[1:80, ]

# Calculate statistics for HC
resultadosEPIC <- data.frame(describe(HC$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(HC$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(HC$age)) # chronological age

# Calculate statistics for HPRE
resultadosEPIC <- data.frame(describe(HPRE$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(HPRE$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(HPRE$age)) # chronological age

# T-student for T0_HC_TN_HPRE comparison
t.test(EPIC ~ Sample_Group, data = T0_HC_TN_HPRE) # for gestational age
t.test(ageAcc.EPIC ~ Sample_Group, data = T0_HC_TN_HPRE) # for acceleration


# --- STATS FOR T0_HC_T0_HPRE ---
# Read file with data
T0_HC_T0_HPRE <- read.csv("T0_HC_T0_HPRE/estGA_T0_HC_T0_HPRE.csv")

# Divide by sample group
HC <- T0_HC_T0_HPRE[61:105, ]
HPRE <- T0_HC_T0_HPRE[1:60, ]

# Calculate statistics for HC
resultadosEPIC <- data.frame(describe(HC$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(HC$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(HC$age)) # chronological age

# Calculate statistics for HPRE
resultadosEPIC <- data.frame(describe(HPRE$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(HPRE$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(HPRE$age)) # chronological age

# T-student for T0_HC_T0_HPRE comparison
t.test(EPIC ~ Sample_Group, data = T0_HC_T0_HPRE) # for gestational age
t.test(ageAcc.EPIC ~ Sample_Group, data = T0_HC_T0_HPRE) # for acceleration


# --- STATS FOR T0_HPRE ---
# Read file with data
T0_HPRE <- read.csv("T0_HPRE/estGA_T0_HPRE.csv")

# Divide by sample group
P1 <- T0_HPRE[1:34, ]
P2 <- T0_HPRE[35:60, ]

# Calculate statistics for P1
resultadosEPIC <- data.frame(describe(P1$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(P1$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(P1$age)) # chronological age

# Calculate statistics for P2
resultadosEPIC <- data.frame(describe(P2$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(P2$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(P2$age)) # chronological age

# T-student for T0_HPRE comparison
t.test(EPIC ~ Sample_Group, data = T0_HPRE) # for gestational age
t.test(ageAcc.EPIC ~ Sample_Group, data = T0_HPRE) # for acceleration


# --- STATS FOR TN_HPRE ---
# Read file with data
TN_HPRE <- read.csv("TN_HPRE/estGA_TN_HPRE.csv")

# Divide by sample group
P1 <- TN_HPRE[1:51, ]
P2 <- TN_HPRE[52:80, ]

# Calculate statistics for P1
resultadosEPIC <- data.frame(describe(P1$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(P1$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(P1$age)) # chronological age

# Calculate statistics for P2
resultadosEPIC <- data.frame(describe(P2$EPIC)) # EPIC clock age
resultadosAcc <- data.frame(describe(P2$ageAcc.EPIC)) # EPIC clock acceleration
resultadosAge <- data.frame(describe(P2$age)) # chronological age

# T-student for TN_HPRE comparison
t.test(EPIC ~ Sample_Group, data = TN_HPRE) # for gestational age
t.test(ageAcc.EPIC ~ Sample_Group, data = TN_HPRE) # for acceleration 


