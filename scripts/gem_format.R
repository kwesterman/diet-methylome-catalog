library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args)==2)) {
  stop("2 arguments needed: exposure_var and cov_suffix", call.=FALSE)
}
envVar <- args[1]  # First argument is the dietary exposure of interest
cov_suffix <- args[2]  # Second argument is the suffix describing the covariates included
foodLib <- c(cocoa="FFD118", n3="NUT_OMEGA", folate="NUT_FOLEQ", soy="FFD49", beans="FFD60",   # For converting dietary factors into FFD codes
             french_fries="FFD99", coffee="FFD112", tea="FFD113", red_wine="FFD115", SSB="FFD145", 
             nuts="FFD134", fat_pct="NUT_TFAT", vitk="NUT_VITK")

covCombos <- list(basic=c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),  # Associate covariate suffixes with specific sets of covariates
		  all=c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5","Aspirin","CVD_med","Menopause","cig_d","smoking_now","drink_d_pwk","BMI","PAS7"))
if(!(cov_suffix %in% names(covCombos))) stop("Covariate list not recognized.", call.=FALSE)  # Check for existence of the requested covariate set

if(!("methylation.txt" %in% list.files("../results/int"))) {
  methCovData <- fread("../data/fram_methyl_cov.csv")
  meth <- t(methCovData)
  colnames(meth) <- meth[1,]
  envData <- read.csv("../data/Off_Exam8_FFQ_diet.csv") %>%
    dplyr::rename(id=shareid) %>%
    t()
  colnames(envData) <- envData[1,]
  common_ids <- as.character(sort(base::intersect(meth[1,], envData[1,])))
  write.table(meth[c(1,24:nrow(meth)),common_ids], "../results/int/methylation.txt", sep="\t", col.names=F, quote=F)
  all_vars <- rbind(meth[1:23,common_ids], envData[,common_ids])
  write.table(all_vars, "../results/int/all_vars.txt", sep="\t", col.names=F, quote=F)
#   cov <- meth[c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),common_ids]
#   cov <- rbind(meth[c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),common_ids],
#                t(envData[common_ids,NULL]))  # Could put nutritional covariates in here
#   write.table(cov, "../results/int/cov.txt", sep="\t", col.names=F, quote=F)
}

if(!(envVar %in% names(foodLib))) stop("dietary data not available")

varData <- read.delim("../results/int/all_vars.txt", header=F, row.names=1)
env <- varData[c("id",foodLib[envVar]),]
write.table(env, paste0("../results/int/",envVar,".txt"), sep="\t", col.names=F, quote=F)
covs <- varData[covCombos[[cov_suffix]],]
write.table(cov, paste0("../results/int/cov_",cov_suffix,".txt"), sep="\t", col.names=F, quote=F)


#covData <- read.delim("../results/int/available_covariates.txt", header=F, row.names=1)
#colnames(covData) <- covData[1,]
#envData <- read.csv("../data/Off_Exam8_FFQ_diet.csv") %>%
#  rename(id=shareid) %>%
#  # mutate(fat_pct=NUT_TFAT/NUT_CALOR) %>%
#  # select(id, fat_pct)
#  select_("id", foodLib[envVar])
#common_ids <- as.character(sort(base::intersect(covData[1,], envData$id)))
#cov <- covData[covCombos[[cov_suffix]],common_ids]
#write.table(cov, paste0("../results/int/cov_",cov_suffix,".txt"), sep="\t", col.names=F, quote=F)
#a <- t(envData)
#colnames(a) <- a[1,]
#env <- a[c("id",foodLib[envVar]),common_ids]
#write.table(env, paste0("../results/int/",envVar,".txt"), sep="\t", col.names=F, quote=F)


# methCovData <- fread("toy2_fram_methyl_cov.csv")
# envData <- read.csv("Off_Exam8_FFQ_diet.csv") %>%
#   rename(id=shareid) %>%
#   select(id, FFD118)
# common_ids <- as.character(sort(intersect(methCovData$id, envData$id)))
# # a <- setNames(data.frame(t(methCovData[,-1,with=F])), methCovData[,1,with=F])
# a <- t(methCovData)
# colnames(a) <- a[1,]
# meth <- a[c(1,24:nrow(a)),common_ids]
# cov <- a[c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),common_ids]
# a <- t(envData)
# colnames(a) <- a[1,]
# env <- a[c("id","FFD118"),common_ids]
# 
# write.table(meth, "methylation.txt", sep="\t", col.names=F, quote=F)
# write.table(cov, "cov.txt", sep="\t", col.names=F, quote=F)
# write.table(env, "env.txt", sep="\t", col.names=F, quote=F)

# methCovData <- fread("toy2_fram_methyl_cov.csv") %>% 
#   data.frame()
# envData <- read.csv("Off_Exam8_FFQ_diet.csv") %>%
#   rename(id=shareid) %>%
#   select(id, FFD118)  # ID and chocolate consumption field from FFD
# common_ids <- data.frame(id=sort(intersect(methCovData$id, envData$id)))  # Sorted IDs common to methylation and environment
# methCovSort <- inner_join(common_ids, methData, by="id")
# meth <- methCovSort %>%
#   select(1, 24:ncol(methCovSort)) %>%
#   t()
# cov <- methCovSort %>%
#   rename(sex=SEX, age=AGE8) %>%
#   select(id, sex, age, pc1, pc2, pc3, pc4, pc5) %>%
#   t()
# env <- inner_join(common_ids, envData, by="id") %>%
#   t()
