setwd("/home/nodotea/LAURA")
getwd()

##############################################################################
##################### EPITHELIAL CELLS PROPORTION ############################
##############################################################################
### Step repeated for each comparison: HC-0 vs HPRE-N, HC-0 vs HPRE-0, HPRE-0 and HPRE-N

#------------------------------------------------------------------------------
# Methylation Matrix of beta values
#------------------------------------------------------------------------------
library("ChAMP")
myLoad <- champ.load("DATOS_21_22", ProbeCutoff=0, arraytype="EPIC") # load data of each comparison
beta_child <- myLoad$beta # methyl matrix
class(beta_child)
pd_child <- myLoad$pd # data of pd file
class(pd_child)


#------------------------------------------------------------------------------
# Estimation of Cell Population
#------------------------------------------------------------------------------
# estimateLC() function use the reference panel of saliva EPIC referenced below:
#' @references{Middleton LYM, et al. Saliva cell type DNA methylation reference panel for epidemiology studies in children. 2020 Sep;}

# The use of the package can be found here:
# https://github.com/bakulskilab/Saliva_Reference/blob/master/13.2%20Child%20Saliva%20Estimation.R

devtools::install_github("hhhh5/ewastools@master")
library(ewastools)
library(data.table)

saliva_est <- estimateLC(meth = beta_child, ref = "salivaEPIC", constrained=TRUE)

# make it a dataframe
saliva_est_df <- as.data.frame(saliva_est)

# join to original pd file
pd_child_cell_est <- cbind(pd_child, saliva_est_df$Epithelial.cells)
colnames(pd_child_cell_est)[9] <- "Epithelial_cells"
pd_child_cell_est
write.table(x = pd_child_cell_est, file = "child_cell_est.csv", sep = ",", 
            row.names = FALSE, col.names = TRUE)
