###########################################################################
##### SCRIPT FOR GENE-SET ENRICHEMNT ANALYSIS WITH MAGMA #######
###########################################################################

###########################
##### STEP A: ANNOTATION ###
###########################

# Annotate SNPS to GENES:
# SNPs annotation of gene regions with window of 35,10 kb.
# SNP list with 1000 european population genomes
# NCBI 37.3 gene database, with brain expresion (GTEx v8 TPM > 0.5)

./magma --annotate window=35,10 --snp-loc ./g1000_eur_rsID/g1000_eur.bim --gene-loc ./Gene_Location/NCBI/NCBI37.3_GENE_NAME.gene.brainGTExV8_05_TPM.loc --out ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID



###############################
##### STEP B: GENE ANALYSIS ###
###############################

## Gene p-values and other gene-level metrics are computed:

./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/SCZ.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/SCZ.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/BIP.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/BIP.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/MDD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/MDD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/ASD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/ASD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/ADHD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/ADHD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/TS.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/TS.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/OCD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/OCD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/ALCOH.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/ALCOH.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/CROSSD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/CROSSD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/CUD.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/CUD.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/NEUR.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/NEUR.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/COGN.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/COGN.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/EA.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/EA.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/PRETERM_BIRTH.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/PRETERM_BIRTH.ws35_10_gene_analysis
./magma --bfile ./g1000_eur_rsID/g1000_eur synonyms=./dbsnp151.synonyms/dbsnp151.synonyms --pval ./magma_summaries/rsID/BIRTH_WEIGHT.magma_rsID.txt use=SNP,P ncol=N --gene-annot ./annotation/NCBI37_brainGTExV8_05_TPM_hg19_eur_ws35_10_rsID.genes.annot --out ./GENE_RESULTS/rsID/BIRTH_WEIGHT.ws35_10_gene_analysis



###################################
##### STEP C: GENE SET ANALYSIS ###
###################################

## Competitive analysis to see if genes with associated variants are enriched significantly in two gene sets, the hipermetilated (negative gene set) and the hipometilated (positive gene set).
## Step repeated with gene sets of the 4 comparisons:

##### T0HC_vs_TNHPRE comparison:
declare -a StringArray_Selection=("SCZ" "BIP" "MDD" "ASD" "ADHD" "TS" "OCD" "ALCOH" "CROSSD" "CUD" "NEUR" "COGN" "EA" "PRETERM_BIRTH" "BIRTH_WEIGHT")
for TRAIT in ${StringArray_Selection[@]}
	do
	./magma \
	--gene-results ./GENE_RESULTS/rsID/$TRAIT.ws35_10_gene_analysis.genes.raw \
	--set-annot ./GENE_SETS/PREM_TEA_DEF/GeneSets_T0HC_vs_TNHPRE \
	--out ./GENE_SET_RESULTS_DEF/T0HC_vs_TNHPRE/$TRAIT.ws35
done

##### T0HC_vs_T0HPRE comparison:
declare -a StringArray_Selection=("SCZ" "BIP" "MDD" "ASD" "ADHD" "TS" "OCD" "ALCOH" "CROSSD" "CUD" "NEUR" "COGN" "EA" "PRETERM_BIRTH" "BIRTH_WEIGHT")
for TRAIT in ${StringArray_Selection[@]}
	do
	./magma \
	--gene-results ./GENE_RESULTS/rsID/$TRAIT.ws35_10_gene_analysis.genes.raw \
	--set-annot ./GENE_SETS/PREM_TEA_DEF/GeneSets_T0HC_vs_T0HPRE \
	--out ./GENE_SET_RESULTS_DEF/T0HC_vs_T0HPRE/$TRAIT.ws35
done

##### T0_HPRE comparison:
declare -a StringArray_Selection=("SCZ" "BIP" "MDD" "ASD" "ADHD" "TS" "OCD" "ALCOH" "CROSSD" "CUD" "NEUR" "COGN" "EA" "PRETERM_BIRTH" "BIRTH_WEIGHT")
for TRAIT in ${StringArray_Selection[@]}
	do
	./magma \
	--gene-results ./GENE_RESULTS/rsID/$TRAIT.ws35_10_gene_analysis.genes.raw \
	--set-annot ./GENE_SETS/PREM_TEA_DEF/GeneSets_T0_HPRE \
	--out ./GENE_SET_RESULTS_DEF/T0_HPRE/$TRAIT.ws35
done

##### TN_HPRE comparison:
declare -a StringArray_Selection=("SCZ" "BIP" "MDD" "ASD" "ADHD" "TS" "OCD" "ALCOH" "CROSSD" "CUD" "NEUR" "COGN" "EA" "PRETERM_BIRTH" "BIRTH_WEIGHT")
for TRAIT in ${StringArray_Selection[@]}
	do
	./magma \
	--gene-results ./GENE_RESULTS/rsID/$TRAIT.ws35_10_gene_analysis.genes.raw \
	--set-annot ./GENE_SETS/PREM_TEA_DEF/GeneSets_TN_HPRE \
	--out ./GENE_SET_RESULTS_DEF/TN_HPRE/$TRAIT.ws35_10 
done

