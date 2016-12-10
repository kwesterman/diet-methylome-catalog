library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args)==2)) {  # Check for exactly 2 arguments
  stop("2 arguments needed: exposure_var and cov_suffix", call.=FALSE)
}
envVar <- args[1]  # First argument is the dietary exposure of interest
cov_suffix <- args[2]  # Second argument is the suffix describing the covariates included
foodLib <- c(cocoa="FFD118", n3="NUT_OMEGA", folate="NUT_FOLEQ", soy="FFD49", beans="FFD60",   # For converting dietary factors into FFD codes
             french_fries="FFD99", coffee="FFD112", tea="FFD113", red_wine="FFD115", SSB="FFD145",
             nuts="FFD134", fat_pct="NUT_TFAT", vitk="NUT_VITK")
covCombos <- list(basic=c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),  # Associate covariate suffixes with specific sets of covariates
                  cov1=c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5","Aspirin","CVD_med","Menopause","cig_d","smoking_now","drink_d_pwk","BMI","PAS7","NUT_CALOR"))
if(!(cov_suffix %in% names(covCombos))) stop("Covariate list not recognized.", call.=FALSE)  # Check for existence of the requested covariate set

## Create methylation and additional variables files if they don't yet exist
if(!("methylation.txt" %in% list.files("../results/int"))) {  # Does the re-formatted methylation data exist at the moment?
  methCovData <- fread("../data/fram_methyl_cov.csv")  # data.table fread for faster read-in
  meth <- t(methCovData)
  colnames(meth) <- meth[1,]
  envData <- read.csv("../data/Off_Exam8_FFQ_diet.csv") %>%  # Environmental exposure data comes from FHS FFQ and derived vars
    dplyr::rename(id=shareid) %>%
    dplyr::select(-dbGaP_Subject_ID) %>%  # To avoid complications/confusion with the FHS ID
    t()
  colnames(envData) <- envData[1,]
  common_ids <- as.character(sort(base::intersect(meth[1,], envData[1,])))  # Get set of IDs common to methylation and dietary data
  write.table(meth[c(1,24:nrow(meth)),common_ids], "../results/int/methylation.txt", sep="\t", col.names=F, quote=F)  # Write only the actual methylation data
  all_vars <- rbind(meth[1:23,common_ids], envData[-1,common_ids])  # All_vars is a combination of the covariate data from the methylation file and the dietary data
  write.table(all_vars, "../results/int/all_vars.txt", sep="\t", col.names=F, quote=F)  # Write a file containing all possible exposure and covariate data
}

## Write experiment-specific data files for use by GEM
if(!(envVar %in% names(foodLib))) stop("dietary data not available")  # Check for existence of the requested dietary factor
varData <- read.delim("../results/int/all_vars.txt", header=F, row.names=1)  # all_vars contains all available covariates and dietary factors
env <- varData[c("id",foodLib[envVar]),]  # Subset to grab only food variable of interest
write.table(env, paste0("../results/int/",envVar,".txt"), sep="\t", col.names=F, quote=F)
covs <- varData[covCombos[[cov_suffix]],]  # Subset to grab only covariates of interest (based on above list)
write.table(covs, paste0("../results/int/cov_",cov_suffix,".txt"), sep="\t", col.names=F, quote=F)

## Run GEM EWAS analysis
library(GEM)
DATADIR <- getwd()
env_file_name <- paste(DATADIR, "../results/int/",paste0(envVar,".txt"), sep=.Platform$file.sep)
covariate_file_name <- paste(DATADIR, "../results/int", paste0("cov_",cov_suffix,".txt"), sep=.Platform$file.sep)
methylation_file_name <- paste(DATADIR, "../results/int/methylation.txt", sep=.Platform$file.sep)
Emodel_pv <- 1
Emodel_result_file_name <- paste0("../results/",envVar,"_",cov_suffix,"_result_Emodel.txt")
Emodel_qqplot_file_name <- paste0("../results/",envVar,"_",cov_suffix,"_qqplot_Emodel.jpg")
GEM_Emodel(env_file_name, covariate_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name, Emodel_qqplot_file_name, savePlot=T)
