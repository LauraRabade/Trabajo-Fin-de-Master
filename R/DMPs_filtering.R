setwd("C:/Users/lalaa/Desktop/UAM/TFM/Analisis/Resultados/Champ")

##############################################################################
#####################        DMPs FILTERING      ############################
##############################################################################
library("dplyr")
library("ggplot2")

# Load DMP data for each comparison:
myDMP_0_HPRE <- "RESULTS_T0_HPRE/DMPs/myDMP.csv"
myDMP_N_HPRE <- "RESULTS_TN_HPRE/DMPs/myDMP.csv"
myDMP_0_HCvs0_HPRE <- "RESULTS_T0_HC_vs_T0_HPRE/DMPs/myDMP.csv"
myDMP_0_HCvsN_HPRE <- "RESULTS_T0_HC_vs_TN_HPRE/DMPs/myDMP.csv"


######## T0_HPRE comparison #########
# to obtain filtered DMPS
DMP_0_HPRE <- read.csv(myDMP_0_HPRE, header = TRUE, sep = ",", dec = ".")
dmp_abs <- DMP_0_HPRE %>% # filter
  mutate(abs_beta = abs(P1_to_P2.deltaBeta)) %>% # create new column with absolute values of column 'P1_to_P2.deltaBeta'
  arrange(desc(abs_beta)) %>%
  filter(abs_beta > 0.2) %>%
  filter(P1_to_P2.adj.P.Val < 0.05) %>%
  select(cg_id, P1_to_P2.deltaBeta, P1_to_P2.adj.P.Val, P1_to_P2.gene, P1_to_P2.CHR, P1_to_P2.MAPINFO, P1_to_P2.feat.cgi)
write.csv(dmp_abs, file = "RESULTS_T0_HPRE/DMPs/DMP_abs_filter.csv", quote = FALSE, row.names = FALSE)

# to obtain how many DMP are hyper and hypomethylated, before and after filtering:
sum(DMP_0_HPRE$P1_to_P2.deltaBeta > 0) #hypo
sum(DMP_0_HPRE$P1_to_P2.deltaBeta < 0) #hyper

sum(dmp_abs$P1_to_P2.deltaBeta > 0) #hypo
sum(dmp_abs$P1_to_P2.deltaBeta < 0) #hyper


######## TN_HPRE comparison #########
DMP_N_HPRE <- read.csv(myDMP_N_HPRE, header = TRUE, sep = ",", dec = ".")
dmp_abs <- DMP_N_HPRE %>% 
  mutate(abs_beta = abs(P1_to_P2.deltaBeta)) %>% 
  arrange(desc(abs_beta)) %>%
  filter(abs_beta > 0.2) %>%
  filter(P1_to_P2.adj.P.Val < 0.05) %>%
  select(cg_id, P1_to_P2.deltaBeta, P1_to_P2.adj.P.Val, P1_to_P2.gene, P1_to_P2.CHR, P1_to_P2.MAPINFO, P1_to_P2.feat.cgi)
write.csv(dmp_abs, file = "RESULTS_TN_HPRE/DMPs/DMP_abs_filter.csv", quote = FALSE, row.names = FALSE)

# to obtain how many DMP are hyper and hypomethylated, before and after filtering:
sum(DMP_N_HPRE$P1_to_P2.deltaBeta > 0) #hypo
sum(DMP_N_HPRE$P1_to_P2.deltaBeta < 0) #hyper

sum(dmp_abs$P1_to_P2.deltaBeta > 0) #hypo
sum(dmp_abs$P1_to_P2.deltaBeta < 0) #hyper


######## T0_HC_vs_T0_HPRE comparison #########
DMP_0_HCvs0_HPRE <- read.csv(myDMP_0_HCvs0_HPRE, header = TRUE, sep = ",", dec = ".")
dmp_abs <- DMP_0_HCvs0_HPRE %>% 
  mutate(abs_beta = abs(HPRE_to_HC.deltaBeta)) %>% 
  arrange(desc(abs_beta)) %>%
  filter(abs_beta > 0.2) %>%
  filter(HPRE_to_HC.adj.P.Val < 0.05) %>%
  select(cg_id, HPRE_to_HC.deltaBeta, HPRE_to_HC.adj.P.Val, HPRE_to_HC.gene, HPRE_to_HC.CHR, HPRE_to_HC.MAPINFO, HPRE_to_HC.feat.cgi)
write.csv(dmp_abs, file = "RESULTS_T0_HC_vs_T0_HPRE/DMPs/DMP_abs_filter.csv", quote = FALSE, row.names = FALSE)

# to obtain how many DMP are hyper and hypomethylated, before and after filtering:
sum(DMP_0_HCvs0_HPRE$HPRE_to_HC.deltaBeta > 0) #hypo
sum(DMP_0_HCvs0_HPRE$HPRE_to_HC.deltaBeta < 0) #hyper
sum(dmp_abs$HPRE_to_HC.deltaBeta > 0) #hypo
sum(dmp_abs$HPRE_to_HC.deltaBeta < 0) #hyper


######## T0_HC_vs_TN_HPRE comparison #########
DMP_0_HCvsN_HPRE <- read.csv(myDMP_0_HCvsN_HPRE, header = TRUE, sep = ",", dec = ".")
dmp_abs <- DMP_0_HCvsN_HPRE %>% 
  mutate(abs_beta = abs(HPRE_to_HC.deltaBeta)) %>% 
  arrange(desc(abs_beta)) %>%
  filter(abs_beta > 0.2) %>%
  filter(HPRE_to_HC.adj.P.Val < 0.05) %>%
  select(cg_id, HPRE_to_HC.deltaBeta, HPRE_to_HC.adj.P.Val, HPRE_to_HC.gene, HPRE_to_HC.CHR, HPRE_to_HC.MAPINFO, HPRE_to_HC.feat.cgi)
write.csv(dmp_abs, file = "RESULTS_T0_HC_vs_TN_HPRE/DMPs/DMP_abs_filter.csv", quote = FALSE, row.names = FALSE)

# to obtain how many DMP are hyper and hypomethylated, before and after filtering:
sum(DMP_0_HCvsN_HPRE$HPRE_to_HC.deltaBeta > 0) #hypo
sum(DMP_0_HCvsN_HPRE$HPRE_to_HC.deltaBeta < 0) #hyper
sum(dmp_abs$HPRE_to_HC.deltaBeta > 0) #hypo
sum(dmp_abs$HPRE_to_HC.deltaBeta < 0) #hyper



##############################################################################
#####################   Plots for probe characteristics    ###################
##############################################################################

######## T0_HPRE comparison #########
DMP_0_HPRE <- DMP_0_HPRE %>% 
  mutate(methylation = if_else(P1_to_P2.deltaBeta < 0, "HyperProbe", "HypoProbe"))
p <- ggplot(DMP_0_HPRE, aes(x = P1_to_P2.feat.cgi, fill = methylation)) +
  geom_bar(position = "dodge") +
  scale_x_discrete("probe") +
  scale_y_continuous("counts") +
  scale_fill_brewer(palette = "Dark2") +
  theme(rect = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        legend.title = element_blank())
ggplotly(p)

######## TN_HPRE comparison #########
DMP_N_HPRE <- DMP_N_HPRE %>% 
  mutate(methylation = if_else(P1_to_P2.deltaBeta < 0, "HyperProbe", "HypoProbe"))
p <- ggplot(DMP_N_HPRE, aes(x = P1_to_P2.feat.cgi, fill = methylation)) +
  geom_bar(position = "dodge") +
  scale_x_discrete("probe") +
  scale_y_continuous("counts") +
  scale_fill_brewer(palette = "Dark2") +
  theme(rect = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        legend.title = element_blank())
ggplotly(p)

######## HC-0 vs HPRE-0 comparison #########
DMP_0_HCvs0_HPRE <- DMP_0_HCvs0_HPRE %>% 
  mutate(methylation = if_else(HPRE_to_HC.deltaBeta < 0, "HyperProbe", "HypoProbe"))
p <- ggplot(dmps3, aes(x = HPRE_to_HC.feat.cgi, fill = methylation)) +
  geom_bar(position = "dodge") +
  scale_x_discrete("probe") +
  scale_y_continuous("counts") +
  scale_fill_brewer(palette = "Dark2") +
  theme(rect = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        legend.title = element_blank())
ggplotly(p)

######## HC-0 vs HPRE-N comparison #########
DMP_0_HCvsN_HPRE <- DMP_0_HCvsN_HPRE %>% 
  mutate(methylation = if_else(HPRE_to_HC.deltaBeta < 0, "HyperProbe", "HypoProbe"))
p <- ggplot(DMP_0_HCvsN_HPRE, aes(x = HPRE_to_HC.feat.cgi, fill = methylation)) +
  geom_bar(position = "dodge") +
  scale_x_discrete("probe") +
  scale_y_continuous("counts") +
  scale_fill_brewer(palette = "Dark2") +
  theme(rect = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        legend.title = element_blank())
ggplotly(p)

