
## SCRIPT SUMMARY STATISTICS

# Load summary statistics to have variants of interest with rsID, N and P

## SCZ - PGC3
dataDiscovery <- read.table('PGC3_leave_out_cibersam', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8)  
dataDiscovery <- temp
temp <- dataDiscovery[,c("SNP", "P", "Neff")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "SCZ.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## BIP - PGC 2021
dataDiscovery <- read.table('BIP_2021_PGC_rsID', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8)  
dataDiscovery <- temp
temp <- dataDiscovery[,c("SNP", "P", "Neff")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "BIP.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## MDD - UKB-PGC 2019
dataDiscovery <- read.table('MDD_2021', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8)  
dataDiscovery <- temp
temp <- dataDiscovery[,c("SNP", "P", "Neff")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "MDD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## TEA - PGC-iPSYCH 2019
dataDiscovery <- read.table('iPSYCH-PGC_ASD_Nov2017', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8)  
dataDiscovery <- temp
dataDiscovery$N = 46350
temp <- dataDiscovery[,c("SNP", "P", "N")]
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "ASD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## TDAH - PGC 2019
dataDiscovery <- read.table('daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8) 
dataDiscovery <- temp
temp <- dataDiscovery[,c("SNP", "P", "Neff")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "ADHD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## TOURETTE - 2018
dataDiscovery <- read.table('TS_Oct2018', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8) 
dataDiscovery <- temp
dataDiscovery$N = 14397
temp <- dataDiscovery[,c("SNP", "P", "N")]
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "TS.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## OCD - 2017 PGC
dataDiscovery <- read.table('ocd_aug2017', header=TRUE, stringsAsFactors = FALSE) 
temp <-subset(dataDiscovery, INFO >= 0.8)
dataDiscovery <- temp
dataDiscovery$N = 9725
temp <- dataDiscovery[,c("SNP", "P", "N")]
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "OCD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## CROSS DISORDER 8 TRAITS - PGC 2019
dataDiscovery <- read.table('pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt', header=TRUE, stringsAsFactors = FALSE) 
temp <- dataDiscovery[,c("ID", "PVAL", "NEFFDIV2")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "CROSSD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## ALCOHOL USE DISORDER - AUDIT
dataDiscovery <- read.table('AUDIT_UKB_2018_AJP.txt', header=TRUE, stringsAsFactors = FALSE) 
library(dplyr)
temp = dataDiscovery %>% filter(grepl('rs', rsid ))
dataDiscovery <- temp
temp <-subset(dataDiscovery, info >= 0.8) 
dataDiscovery <- temp
temp <- dataDiscovery[,c("rsid", "p_P", "N")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "ALCOH.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## CANNABIS USE DISORDER 2020 
dataDiscovery <- read.table('CUD_EUR_full_public_11.14.2020', header=TRUE, stringsAsFactors = FALSE) 
temp <- dataDiscovery[,c("SNP", "P", "N")]
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "CUD.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## NEUROTICISM - Luciano 2019  
dataDiscovery <- read.table('Neuroticism_Full.txt', header=TRUE, stringsAsFactors = FALSE) 
dataDiscovery$N = 329821
temp <- dataDiscovery[,c("MarkerName", "Pval", "N")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "NEUR.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)


###########################################################################################


## GENERAL COGNITITVE FUNTION - Davies 2018 Natcomms
dataDiscovery <- read.table('Davies2018_OPEN_DATASET_summary_results.txt', header=TRUE, stringsAsFactors = FALSE) 
dataDiscovery$N = 300486
temp <- dataDiscovery[,c("Markername", "P.value", "N")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "COGN.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## EDUCATIONAL ATTAINMENT - Obkay 2022
dataDiscovery <- read.table('EA4_additive_excl_23andMe', header=TRUE, stringsAsFactors = FALSE) 
dataDiscovery$N = 3037499
temp <- dataDiscovery[,c("rsID", "P", "N")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "EA.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)



## FETAL EARLY PRETERM BIRTH  - EGG 2019
dataDiscovery <- read.table('Fetal_early_preterm_birth_NComms2019.txt', header=TRUE, stringsAsFactors = FALSE) 
dataDiscovery$N = dataDiscovery$N_case + dataDiscovery$N_control
temp <- dataDiscovery[,c("Rsid", "P", "N")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "PRETERM_BIRTH.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)


## FETAL BIRTH WEIGHT - EGG 2019
dataDiscovery <- read.table('Fetal_BW_European_meta.NG2019.txt', header=TRUE, stringsAsFactors = FALSE) 
temp <- dataDiscovery[,c("rsid", "p", "n")]
colnames(temp)<-c("SNP", "P", "N")
dataDiscovery <- temp
head(dataDiscovery)
write.table(dataDiscovery, "BIRTH_WEIGHT.magma_rsID.txt", sep="\t", quote=FALSE, row.names=FALSE)


