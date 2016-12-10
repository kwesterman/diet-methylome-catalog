# analysis of gem results

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(tibble)

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args)>=2)) {  # Check for exactly 2 arguments
  stop("Multiple arguments needed: all env variables and cov_suffix", call.=FALSE)
}

envVars <- args[1:(length(args)-1)]  # First argument is the dietary exposure of interest
cov_suffix <- args[length(args)]  # Second argument is the suffix describing the covariates included


anno <- data.frame(Other) %>% rownames_to_column()
islands <- data.frame(Islands.UCSC) %>% rownames_to_column()

add_annotation <- function(envVar, cov_suffix, anno, islands) {
  res <- read.delim(paste0("../results/", envVar, "_", cov_suffix, "_result_Emodel.txt"))
  #a <- read.delim("../results/cocoa_result_Emodel.txt")
  res_anno <- inner_join(res, anno, by=c("cpg"="rowname")) %>%
    inner_join(islands, by=c("cpg"="rowname")) %>%
    dplyr::select(cpg, beta, pvalue, FDR, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_Island)
  #res_anno <- inner_join(res, anno, by=c("cpg"="rowname"))
  res_annoTop <- filter(res_anno, FDR<0.05)
  write.csv(res_annoTop, paste0("../results/", envVar, "_", cov_suffix, "_result_Emodel_Annotated.csv"), row.names=F, quote=F)
}

lapply(envVars, function(x) add_annotation(x, cov_suffix, anno, islands))
#env <- read.delim("../env.txt", stringsAsFactors=F, header=F) %>%
#  t() %>%
#  data.frame(stringsAsFactors=F)
#colnames(env) <- env[1,]
#env <- env[-1,]
#meth <- read.delim("../methylation.txt", stringsAsFactors=F) %>%  # 1-import so row names come out right
#  dplyr::slice(25:30)
#  # filter(id %in% cocoaTop$cpg)
#meth <- setNames(data.frame(t(meth[,-1])), meth[,1]) %>%  # Convert so numeric don't become factor/character
#  rownames_to_column(var="id") %>%
#  mutate(id=gsub("X","",id))
#
#train <- inner_join(meth, env, by=c("id"="id")) %>%
#  dplyr::select(-id)
#train$FFD118 <- as.numeric(train$FFD118)
#cat(nrow(train))
#lm1 <- lm(FFD118~., data=train)


