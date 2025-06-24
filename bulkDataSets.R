library(DESeq2)
library(dplyr)
library(tidyr)
library(tidyverse)

# ------------------------------ starting with GSE162979 ----------------------------- #
# only interested in the acute injury
# macrophages during salmonella infection
# define the data dir ----
data_dir_985 <- '/Users/sm2949/Downloads/GSE224985_RAW/'
# define files and vars -----
file_info <- list(
  noninfected4hpi_rep1 = 'GSM7036521_6_NonINFectees_4hpi.geneCounts.txt',
  noninfected4hpi_rep2 = 'GSM7036522_10_NonINFectees_4hpi.geneCounts.txt',
  noninfected4hpi_rep3 = 'GSM7036523_14_NonINFectees_4hpi.geneCounts.txt',
  noninfected4hpi_rep4 = 'GSM7036524_18_NonINFectees_4hpi.geneCounts.txt',
  infected4hpi_rep1    = 'GSM7036525_5_INFectees_4hpi.geneCounts.txt',
  infected4hpi_rep2    = 'GSM7036526_9_INFectees_4hpi.geneCounts.txt',
  infected4hpi_rep3    = 'GSM7036527_13_INFectees_4hpi.geneCounts.txt',
  infected4hpi_rep4    = 'GSM7036528_17_INFectees_4hpi.geneCounts.txt'
)

# load and process each file ----
processed_list <- lapply(names(file_info), function(varname) {
  # read in
  df <- read.csv(file.path(data_dir_985, file_info[[varname]]), sep = "\t", header = FALSE)
  # clean gene names - only keep gene symbol
  df$V1 <- sub(".*\\|", "", df$V1)
  # remove duplicates 
  df <- df[!duplicated(df$V1), ]
  # rename columns
  colnames(df) <- c("gene", varname)
  return(df)
})

# merge together ----
salmonella_counts <- Reduce(function(x, y) merge(x, y, by = "gene", sort = FALSE), processed_list)
rownames(salmonella_counts) <- salmonella_counts[,1]
salmonella_counts[,1] <- NULL

# trying DESEQ2 ----
# create the coldata
colData <- data.frame(sample = colnames(salmonella_counts))
rownames(colData) <- colData[,1]
# create the cols
colData <- colData %>%
  separate(sample, into = c("condition", "replicate"), sep = "4hpi_")
# set cols as factors 
colData$condition <- factor(colData$condition)
colData$replicate <- factor(colData$replicate)

# create the obj
dds <- DESeqDataSetFromMatrix(countData = salmonella_counts,
                              colData = colData,
                              design = ~ replicate + condition) # var of interest at end of formula

# pre-filtering
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# set the reference level
dds$condition <- relevel(dds$condition, ref = "noninfected")

# variance stabalizing transformation
vsd <- vst(dds, blind=FALSE)

# QC
plotPCA(vsd, intgroup=c("condition", "replicate"))

# DE analysis
dds <- DESeq(dds)
res <- results(dds)

# apply log fold shrinkage
resLFC <- lfcShrink(dds, coef="condition_infected_vs_noninfected", type="apeglm")

# order
resOrdered <- resLFC[order(resLFC$pvalue),]
salmonella_results <- as.data.frame(resOrdered)

# save
write.csv(salmonella_results, file = "/Users/sm2949/Desktop/salmonella_results.csv", row.names = TRUE)

# ------------------------------ GSE140810 ----------------------------- #
# acute neural injury in zebrafish larvae
neural_results <- read.csv('/Users/sm2949/Downloads/GSE140810_extended_summary_gene.csv')
rownames(neural_results) <- neural_results[,1]
neural_results[,1] <- NULL
neural_results <- neural_results[,c(16:20)]
# save
write.csv(neural_results, file = "/Users/sm2949/Desktop/neural_results.csv", row.names = TRUE)

# ------------------------------ GSE100029 ------------------------------ #
# macrophages in heart injury - ventricular resection, acute injury
heart_data <- read.csv('/Users/sm2949/Downloads/GSE100029_count.table.txt', sep='\t')
heart_counts <- heart_data[,c(5:10)]

# create the coldata
colData <- data.frame(sample = colnames(heart_counts))
rownames(colData) <- colData[,1]
# create the cols
colData <- colData %>%
  mutate(
    last_part = str_extract(sample, "[^_]+$"),                 
    condition = str_extract(last_part, "[A-Za-z]+"),           
    replicate = str_extract(last_part, "[0-9]+")               
  ) %>%
  select(-last_part) 
colData[,1] <- NULL
# set cols as factors 
colData$condition <- factor(colData$condition)
colData$replicate <- factor(colData$replicate)

# create the obj
dds <- DESeqDataSetFromMatrix(countData = heart_counts,
                              colData = colData,
                              design = ~ replicate + condition) # var of interest at end of formula

# pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# set the reference level
dds$condition <- relevel(dds$condition, ref = "sham")

# variance stabalizing transformation
vsd <- vst(dds, blind=FALSE)

# QC
plotPCA(vsd, intgroup=c("condition", "replicate"))

# DE analysis
dds <- DESeq(dds)
res <- results(dds)

# apply log fold shrinkage
resLFC <- lfcShrink(dds, coef="condition_resec_vs_sham", type="apeglm")

# order
resOrdered <- resLFC[order(resLFC$pvalue),]
heart_results <- as.data.frame(resOrdered)

# save
write.csv(heart_results, file = "/Users/sm2949/Desktop/heart_results.csv", row.names = TRUE)

