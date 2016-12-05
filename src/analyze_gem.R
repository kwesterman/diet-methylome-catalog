# analysis of gem results

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


get_topCpGs <- function(topCpGs) {
  # Read in specific columns (CpGs) of interest from raw methylation data
  file <- "../data/_0for_Kenny_2016Nov21/fram_methyl_cov.csv"
  con <- file(description=file, open="r")  # Raw methylation data, CpGs across the top
  cols <- which(scan(file=con, nlines=1, what=character(), sep=",") %in% topCpGs)  # Get #s of relevant cols
  cmd <- paste0("cut -d ',' -f 1,", paste(cols, collapse=","),  # Prepare "cut" bash cmd to get specific cols
                " ../data/_0for_Kenny_2016Nov21/fram_methyl_cov.csv")
  read.csv(pipe(cmd))  # Execute bash cmd and import the resulting table
}


# file <- "../data/_0for_Kenny_2016Nov21/toy_fram_methyl_cov.csv"
# con <- file(description=file, open="r")
# com <- paste("wc -l ", file, " | awk '{ print $1 }'", sep="")
# n <- system(command=com, intern=TRUE)
# relevant_rows <- list(scan(file=con, nlines=1, what=character(), sep=","))
# for (i in 2:n) {
#   tmp <- scan(file=con, nlines=1, sep=",", quiet=T)
#   if (tmp[1] %in% cgs_of_interest) {relevant_rows <- c(relevant_rows, tmp)}
#   print(tmp[1])
# }


a <- read.delim("../results/vitk_result_Emodel.txt")
anno <- data.frame(Other) %>% rownames_to_column()
islands <- data.frame(Islands.UCSC) %>% rownames_to_column()
res <- inner_join(a, anno, by=c("cpg"="rowname")) %>%
  inner_join(islands, by=c("cpg"="rowname")) %>%
  dplyr::select(cpg, beta, pvalue, FDR, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_Island)
resTop <- filter(res, FDR<0.2)
env <- read.delim("../env.txt", stringsAsFactors=F, header=F) %>%
  t() %>%
  data.frame(stringsAsFactors=F)
colnames(env) <- env[1,]
env <- env[-1,]
system.time(meth11 <- get_topCpGs(resTop$cpg))
meth <- read.delim("../methylation.txt", stringsAsFactors=F) %>%  # 1-import so row names come out right
  dplyr::slice(25:30)
  # filter(id %in% cocoaTop$cpg)
meth <- setNames(data.frame(t(meth[,-1])), meth[,1]) %>%  # Convert so numeric don't become factor/character
  rownames_to_column(var="id") %>%
  mutate(id=gsub("X","",id))

train <- inner_join(meth11, env, by=c("id"="id")) %>%
  dplyr::select(-id)
cat(nrow(train))
lm1 <- lm(NUT_OMEGA~., data=train)


