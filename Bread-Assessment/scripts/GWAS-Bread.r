#########################################################
# Author: Rafat Eissa 
# Date: 20th Sept 2024
# Project: Bread Wheat Assesment 
#########################################################

##################################################################################################################
# 1. Calculating Broad-Sense Heritability 
##################################################################################################################
library(lme4)

# Example dataset
data <- read.csv("wheat_yield_data.csv")

# Fit the linear mixed model
model <- lmer(Yield ~ (1|Genotype) + (1|Environment) + (1|Genotype:Environment), data = data)

# Get the variance components
var_comp <- as.data.frame(VarCorr(model))
genetic_var <- var_comp$vcov[1]
env_var <- var_comp$vcov[2]
gxe_var <- var_comp$vcov[3]
residual_var <- attr(VarCorr(model), "sc")^2

# Calculate total phenotypic variance
total_var <- genetic_var + env_var + gxe_var + residual_var

# Calculate broad-sense heritability
H2 <- genetic_var / total_var
H2

##################################################################################################################
# 1. GWAS using GAPIT-R
##################################################################################################################

#Intaliing GAPIT
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#Importing Genotyping data
myG<-read.table("bread_wheat_genotype.hmp.txt",head = FALSE)
#Import Phenotyping data
myY <- read.table("bread_wheat_phenotype.txt", head = TRUE)
#Running GWAS
myGAPIT<-GAPIT(
  Y=myY,
  G=myG,
  PCA.total=5,
  SNP.MAF =0.05,
  Multiple_analysis = TRUE,
  model=c("GLM","MLM"),
  Random.model=FALSE
)





