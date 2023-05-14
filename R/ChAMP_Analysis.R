setwd("/home/nodotea/LAURA/ChAMP_Def_analisis")
getwd()

###############################################################################
################         ChAMP analysis pipeline         ######################        
###############################################################################
# Pipeline done for each of the 4 comparisons: 
# HC-0 vs HPRE-N, HC-0 vs HPRE-0, HPRE-0 and HPRE-N.

#------------------------------------------------------------------------------
# LOAD DATA
#------------------------------------------------------------------------------
library("ChAMP")
myLoad <- champ.load("DATOS_21_22_DEF", #load idat and pd files of each comparison
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

#------------------------------------------------------------------------------
# QUALITY CONTROL
#------------------------------------------------------------------------------
myQC <- champ.QC(beta = myLoad$beta,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./RESULTS_T0_vs_TN_HPRE/CHAMP_QCimages/")
QC.GUI(beta=myLoad$beta,
       arraytype="EPIC")

#------------------------------------------------------------------------------
# NORMALIZATION
#------------------------------------------------------------------------------
# Using BMIQ Method
myNorm <- champ.norm(beta=myLoad$beta,
                     method="BMIQ",
                     plotBMIQ=TRUE, 
                     arraytype="EPIC", 
                     cores=16,
                     resultsDir="./RESULTS_T0_vs_TN_HPRE/CHAMP_Normalization/")
QC.GUI(beta=myNorm,
       arraytype="EPIC")

#------------------------------------------------------------------------------
# SVD (Singular Value Decomposition analysis)
#------------------------------------------------------------------------------
mySVD <- champ.SVD(beta=myNorm %>% as.data.frame(), 
                   pd=myLoad$pd, 
                   PDFplot=TRUE,
                   Rplot=TRUE,
                   resultsDir="./RESULTS_T0_vs_TN_HPRE/CHAMP_SVDimages/") 
write.table(x = mySVD, file = "./RESULTS_T0_vs_TN_HPRE/CHAMP_SVDimages/mySVD.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE)

#------------------------------------------------------------------------------
# BATCH EFFECTS CORRECTION
#------------------------------------------------------------------------------
myCombat <- champ.runCombat(beta=myNorm,
                            pd=myLoad$pd,
                            variablename="Sample_Group",
                            batchname=c("Epithelial_cells", "Array"),
                            logitTrans=TRUE)

# SVD plot to check batch effects correction has been done:
mySVD_corrected <- champ.SVD(beta = myCombat %>% as.data.frame(),
                             pd=myLoad$pd, 
                             PDFplot=TRUE,
                             Rplot=TRUE,
                             resultsDir="./RESULTS_T0_vs_TN_HPRE/CHAMP_SVDimages_combat/")

write.table(x = mySVD_corrected, file = "./RESULTS_T0_vs_TN_HPRE/CHAMP_SVDimages_combat/mySVD_corrected.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE)


#------------------------------------------------------------------------------
# CALCULATION OF DMPs
#------------------------------------------------------------------------------
myDMP <- champ.DMP(beta = myCombat,
          pheno = myLoad$pd$Sample_Group,
          compare.group = NULL,
          adjPVal = 0.05,
          adjust.method = "fdr", # default method is BH
          arraytype = "EPIC")
head(myDMP[[1]])
DMP.GUI(DMP=myDMP[[1]],
        beta=myCombat,
        pheno=myLoad$pd$Sample_Group)

write.table(x = myDMP, file = "./RESULTS_T0_vs_TN_HPRE/DMPs/myDMP.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE)


#------------------------------------------------------------------------------
# CALCULATION OF DMRs
#------------------------------------------------------------------------------
# Before doing this step, load the function found in champDMR_funct.R script.

myDMR_DMRCate <- champ.DMR(beta=myCombat, 
                           pheno=myLoad$pd$Sample_Group,
                           compare.group=NULL,
                           method="DMRcate", # Method to find DMRs
                           arraytype="EPIC",
                           minProbes=3, 
                           adjPvalDmr=0.05,
                           cores=16,
                           ## following parameters are specifically for DMRcate method:
                           rmSNPCH=T, # Filters a matrix of M-values (or beta values) by distance to SNP
                           fdr=0.05, # FDR cutoff 
                           dist=3, # Maximum distance (from CpG to SNP) of probes to be filtered out
                           mafcut=0.05,
                           lambda=1000,
                           C=2)
head(myDMR_DMRCate$DMRcateDMR)
DMR.GUI(DMR=myDMR_DMRCate,
        beta=myCombat,
        pheno=myLoad$pd$Sample_Group,
        runDMP=FALSE,
        compare.group=NULL,
        arraytype="EPIC")

write.table(x = myDMR_DMRCate, file = "./RESULTS_T0_vs_TN_HPRE/DMRs/myDMR_DMRCate.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE)

