library(data.table)
library(dplyr)

# setwd("~/ordovas/_0for_Kenny_2016Nov21/")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
envVar <- args[1]
foodLib <- c(cocoa="FFD118", n3="NUT_OMEGA", folate="NUT_FOLEQ", soy="FFD49", beans="FFD60", 
             french_fries="FFD99", coffee="FFD112", tea="FFD113", red_wine="FFD115", SSB="FFD145", 
             nuts="FFD134", fat_pct="NUT_TFAT", vitk="NUT_VITK")

if(!("methylation.txt" %in% list.files(".."))) {
  methCovData <- fread("../data/_0for_Kenny_2016Nov21/toy_fram_methyl_cov.csv")
  meth <- t(methCovData)
  colnames(meth) <- meth[1,]
  envData <- read.csv("../data/_0for_Kenny_2016Nov21/Off_Exam8_FFQ_diet.csv") %>%
    dplyr::rename(id=shareid)
  common_ids <- as.character(sort(base::intersect(meth[1,], envData$id)))
  write.table(meth[c(1,24:nrow(meth)),common_ids], "methylation.txt", sep="\t", col.names=F, quote=F)
  write.table(meth[1:23,common_ids], "available_covariates.txt", sep="\t", col.names=F, quote=F)
  cov <- meth[c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),common_ids]
#   cov <- rbind(meth[c("id","SEX","AGE8","pc1","pc2","pc3","pc4","pc5"),common_ids],
#                t(envData[common_ids,NULL]))  # Could put nutritional covariates in here
  write.table(cov, "cov.txt", sep="\t", col.names=F, quote=F)
}

if(!(envVar %in% names(foodLib))) stop("dietary data not available")
covData <- read.delim("../cov.txt", header=F)
envData <- read.csv("../data/_0for_Kenny_2016Nov21/Off_Exam8_FFQ_diet.csv") %>%
  dplyr::rename(id=shareid) %>%
#   mutate(fat_pct=NUT_TFAT/NUT_CALOR) %>%
#   select(id, fat_pct)
  select_("id", foodLib[envVar])
common_ids <- as.character(sort(base::intersect(covData[1,], envData$id)))
a <- t(envData)
colnames(a) <- a[1,]
env <- a[c("id",foodLib[envVar]),common_ids]
write.table(env, paste0("../",envVar,".txt"), sep="\t", col.names=F, quote=F)


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
