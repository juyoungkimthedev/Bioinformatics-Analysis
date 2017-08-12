# hgu133a: treated and control samples

# read in the cel files
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/hgu133a")
library(affy)
hgu133a_affy <- ReadAffy()

# change margin for boxplots
op0 <- par()
op1 <- op0$mar
op1[1] <- 10 # change bottom line to 9.5 from 5.1
par(mar = op1)

# boxplot before normalization
# log2 intensities boxplot of hgu133a_affy is created
boxplot(hgu133a_affy,
        main = "HG-U133A before RMA normalization",
        las = 2,
        ylab = "Intensity",
        cex.axis = 0.8)

# robust multi-array normalization
hgu133a_normalized <- rma(hgu133a_affy)
# extract expression matrix
hgu133a <- exprs(hgu133a_normalized)
# probe names limited to 91952_at
hgu133a <- hgu133a[1:22215,]

# boxplot after normalization
boxplot(hgu133a,
        main = "HG-U133A after RMA normalization",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# save the hgu133a file
saveRDS(hgu133a, "./hgu133a.rds")
hgu133a # use this expression set when combining with hthgu133a set

# now you have normalized treated and control samples from hgu133a chip type
# use the expression matrix and combine it with the matrix from hthgu133a
# qunatile normalization will be performed
# after merging of matrix, extract the treated and control samples and check for batch effect etc.
# than differential expression analysis follows

# hthgu133a: treate and control samples
# read in the cel files
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/hthgu133a")
hthgu133a_affy <- ReadAffy()

# boxplot before normalization
boxplot(hthgu133a_affy,
        main = "HTHG-U133A before RMA normalization",
        las = 2,
        ylab = "Intensity",
        cex.axis = 0.8)

# robust multi-array normalization
hthgu133a_normalized <- rma(hthgu133a_affy)
# extract expression matrix
hthgu133a <- exprs(hthgu133a_normalized)
# probe names limited to 91952_at
hthgu133a <- hthgu133a[1:22215,]

# check dimensions for both type chips
dim(hgu133a); dim(hthgu133a)

# boxplot after normalization
boxplot(hthgu133a,
        main = "HTHG-U133A after RMA normalization",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# save the hthgu133a file
saveRDS(hthgu133a, "./hthgu133a.rds")
hthgu133a # use this expression when combining with hgu133a set

# merge hgu133a DataFrame and hthgu133a DataFrame
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/merged dataframe")
merged_DataFrame <- cbind(hgu133a, hthgu133a)
dim(merged_DataFrame)

# save the merged file
saveRDS(merged_DataFrame, "./merged_DataFrame before quantile normalization")

# boxplot after merging DataFrame
boxplot(merged_DataFrame,
        main = "Merged DataFrame before quantile normalization (hgu133a + hthgu133a)",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# quantile normalization
library(limma)
q_merged_DataFrame <- normalizeQuantiles(merged_DataFrame)

# boxplot after quantile normalization
boxplot(q_merged_DataFrame,
        main = "Merged DataFrame after quantile normalization (hgu133a + hthgu133a)",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# save the file
saveRDS(q_merged_DataFrame, "./q_merged_DataFrame")

# from here what you going to do is:
# Detect for batch effect using PCA batch id 
# perform ComBat analysis
# select dataframe for differential expression analysis
# perform and select the list of genes

# for the above purpose you can use the following:
DataFrame_GT <- q_merged_DataFrame

# I made mistake in C5000 naming with one control file thus must be fixed
# i.e wrong control must be deleted and re-saved "C500024035100021608460.A06.CEL"
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/merged dataframe")
DataFrame_GT <- readRDS("DataFrame_GT.rds")
DataFrame_GT <- as.data.frame(DataFrame_GT)
DataFrame_GT$C500024035100021608460.A06.CEL <- NULL
#saveRDS(DataFrame_GT, "./DataFrame_GT.rds")
dim(DataFrame_GT) # 22215 by 203

# next step is to average the necessary controls for each drug sample
# after averaging the controls we can have full data frame once again and check for variability
# after that check for batch effect etc

# extract control columns from the DataFrame_GT
control_df <- grep("^[C5]", colnames(DataFrame_GT), value = T)
# control_df <- DataFrame_GT[,control_df] # extract controls
# control_df <- control_df[,13:ncol(control_df)] # only select HTHGU133A
# use control_df to average the values later: hthgu133a type samples are the only ones needed to be averaged
control_df

# 1. first of all ComBat and batch effect check first + boxplot the whole values and second of all average the values + boxplot
DataFrame_GT # use this data for everything else

# batch id retrieval
library(dplyr)
library(stringr)
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/batch effect")
batchid <- read.delim("cmap_instances_01.txt", check.names = F)
batchid <- filter(batchid, cell2 == "MCF7")
batchid$perturbation_scan_id <- str_replace(batchid$perturbation_scan_id, pattern = "'", "")
batchid$perturbation_scan_id <- paste0(batchid$perturbation_scan_id, ".CEL")
colnames_GT <- colnames(DataFrame_GT)
treated_samples <- grep("^[T]", colnames_GT, value = T)
treated_samples <- str_replace(treated_samples, pattern = "T", replacement = "")
batchlevel <- which(batchid$perturbation_scan_id %in% treated_samples)
batchlevel <- batchid[batchlevel, ]
write.table(batchlevel, file = "./batchlevel.xls", quote = F, sep = "\t", col.names = NA)
# do some manual work in Excel to attach batchid to the controls for both hgu133a and hthgu133a
control <- as.data.frame(control_df)
write.table(control, file = "./control.xls", quote = F, sep = "\t", col.names = NA)

# label the samples with batch id
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/batch effect")

# identify batch level for each sample
batch_df <- as.data.frame(t(DataFrame_GT))
saveRDS(batch_df, "./batch_df.rds")

# retrieve batch level
t_batch <- read.delim("batchlevel.txt", check.names = F) # read in treated sample batch id
t_batch$perturbation_scan_id <- paste0("T", t_batch$perturbation_scan_id)
t_batch$cmap_name <- NULL
t_batch$cell2 <- NULL
t_batch$array3 <- NULL
t_batch$vehicle_scan_id4 <- NULL
t_batch <- t_batch[, c(2,3)]
c_batch <- read.delim("control batchid.txt", check.names = F) # read in control sample batch id
colnames(c_batch) <- c("batch_id", "perturbation_scan_id")
tc_batch <- rbind(t_batch, c_batch)
colnames(tc_batch)[2] <- "batch.r"
tc_batch # use this data frame

# merge the data frame
test <- match(rownames(batch_df), tc_batch$batch.r, nomatch = 0)
tc_batch_test <- tc_batch[test,]
colnames(tc_batch_test)[2] <- "common"
merged_DF <- cbind(tc_batch_test, batch_df)
rownames(merged_DF) <- merged_DF$common
merged_DF$common <- NULL
merged_DF # use this

# made a mistake by putting 4055 instead of 1055 must be replaced
merged_DF$batch_id <- as.character(merged_DF$batch_id)
ttt <- which(merged_DF$batch_id == "4055")
merged_DF$batch_id <- gsub("4055", "1055", merged_DF$batch_id)
saveRDS(merged_DF, "./merged_DF.rds")
merged_DF <- readRDS("merged_DF.rds")

# ComBat
library(sva)
bat_id <- merged_DF$batch_id
com <- ComBat(t(merged_DF[,-1]), bat_id, mod = NULL, par.prior = T, prior.plots = F) # t(merged_DF[,-1] columns: CEL names; rows: gene names
t.test(merged_DF[,-1],t(com)) # since p-value is bigger than 0.05 we cannot conclude there is significant difference that exists; we do not reject null hypothesis
# don't really need to use ComBat to remove batch effect at this stage

# boxplot
box_test <- as.data.frame(t(merged_DF[,-1]))
boxplot(box_test,
        main = "DataFrame before ComBat",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# boxplot after ComBat
box_test_ComBat <- com
boxplot(box_test_ComBat,
        main = "DataFrame after ComBat",
        las = 2, 
        ylab = "Intensity",
        cex.axis = 0.8)

# now use merged_DF to proceed with averaging the controls
merged_DF
# merged_DF has batchid column as well; if you'd like to remove it,
merged_DF[, -1]

# averaging the control values and test for consistency using boxplot
# if nothing much changes, use that to derive differential expression analyses
merged_DF01 <- merged_DF[, -1]

# CELS that must be averaged
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/averaging the controls")
ave_control <- read.delim("avereage columns.txt", check.names = F)
ave_control$Average1 <- as.character(ave_control$Average1)
ave_control$Average2 <- as.character(ave_control$Average2)
ave_control$Average3 <- as.character(ave_control$Average3)
ave_control$Average4 <- as.character(ave_control$Average4)
ave_control$Average5 <- as.character(ave_control$Average5)
ave_control_final <- c(ave_control$Average1,ave_control$Average2,ave_control$Average3,ave_control$Average4,ave_control$Average5)

# test
test_df <- merged_DF01
test_df01 <- as.data.frame(t(test_df))
test_df01$ave1 <- rowMeans(test_df01[, ave_control$Average1], na.rm = T)
test_df01$ave2 <- rowMeans(test_df01[, ave_control$Average2], na.rm = T)
test_df01$ave3 <- rowMeans(test_df01[, ave_control$Average3], na.rm = T)
test_df01$ave4 <- rowMeans(test_df01[, ave_control$Average4], na.rm = T)
test_df01$ave5 <- rowMeans(test_df01[, ave_control$Average5], na.rm = T)
# delete the other columns and boxplot
test_df02 <- test_df01
c1 <- colnames(test_df02)
c2 <- !(c1 %in% ave_control_final)
test_df02 <- test_df02[, c2]
# boxplot data with averaged control values to see if the data is different
# if data is not different continue like this
boxplot(test_df02,
        main = "Boxplot with averaged controls",
        las = 2, 
        ylab = "Intensity")

# test quantile normalization
library(limma)
ttttt <- normalizeQuantiles(test_df02)
boxplot(ttttt,
        main = "Boxplot with averaged controls - after Q normalization",
        las = 2, 
        ylab = "Intensity")

# check with test_df02 and ttttt which one selects differential expression better
# first of all just try with the test_df02 and go on
test_df02 # data without final quantile normalization (averaged values)
saveRDS(test_df02, "./test_df02.rds")
DEG_df_NQ <- readRDS("test_df02.rds")
ttttt # data after final quantile normalization (averaged values)
DEG_df_Q <- readRDS("ttttt.rds")
saveRDS(ttttt, "./ttttt.rds")

# boxplot after quantile normalization
boxplot(DEG_df_Q,
        main = "Boxplot with averaged controls - after Q normalization",
        las = 2, 
        ylab = "Intensity")
# boxplot with no quantile normalization
boxplot(DEG_df_NQ,
        main = "Boxplot with averaged controls - No Q normalization",
        las = 2, 
        ylab = "Intensity")

# now continue with DEG analysis by starting with design matrix
# DEG analyses must be done using data from DEG_df_NQ and DEG_df_Q
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Data preprocessing ver2/averaging the controls")

# model matrix 
# contrast
# perform differential expression analyses

# annotation of genes first
# work with DEG_df_NQ
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Gene name annotation")
myAnnot <- readRDS("myAnnot.rds")
myAnnot01 <- myAnnot[as.character(myAnnot[,"SYMBOL"]) != "NA",]
myAnnot02 <- myAnnot01[,c(1,2)]
myAnnot02$ACCNUM <- NULL
myAnnot03 <- myAnnot02
colnames(myAnnot03) <- "GENE"
DEG_df_NQ01 <- DEG_df_NQ[rownames(myAnnot03),]
dim(DEG_df_NQ01)

# average duplicated names
idlist <- tapply(row.names(myAnnot03), as.character(myAnnot03$GENE), c)
DEG_df_NQ02 <- t(sapply(names(idlist), function(x) colMeans(DEG_df_NQ01[idlist[[x]], ])))
dim(DEG_df_NQ02)

# glycosyltransferase gene names retrieval
gt <- read.table("GT Genes (210).txt", header = T, sep = "\t")
DEG_df_NQ03 <- DEG_df_NQ02[rownames(DEG_df_NQ02) %in% gt$Gene.Names,]
dim(DEG_df_NQ03)

DEG <- DEG_df_NQ03
saveRDS(DEG_df_NQ03, "./newDEG.rds")
DEG <- readRDS("newDEG.rds")

# boxplot after GT selection
#boxplot(DEG_df_NQ03,
#        main = "Boxplot after GT selection",
#        las = 2,
#        ylab = "Intensity",
#        cex.axis = 0.8)

#library(limma)
#test_q <- normalizeQuantiles(DEG_df_NQ03)
#boxplot(test_q,
#        main = "Boxplot after GT selection with quantile normalization",
#        las = 2,
#        ylab = "Intensity",
#        cex.axis = 0.8)

#View(DEG_df_NQ03)
#View(test_q)

# boxplot of DEG_df_NQ2 shows good result but after DEG_df_NQ3 it shuffles a bit, need to re-study
# from here use the DEG_df_NQ03 and quantile normalizae again and perform DEG analyses
#new_DEG <- normalizeQuantiles(DEG_df_NQ03)
#boxplot(new_DEG)
#dim(new_DEG)
#dim(DEG_df_NQ03)
#new_DEG

# for test purpose after coming back from the practical re run everything from the beginning by doing annotation first thing and averaging everything
# do testing later

# DEG analyses
# use DEG
DEG
dim(DEG)

# first differential expression analysis
# retrieve colnames for 15-prostaglandin J5
eset1_names <- c("T5500024024213121906562.C03.CEL",
                 "T5500024030403071907253.C03.CEL",
                 "T5500024035100021608460.C03.CEL",
                 "T5500024037497121008341.C03.CEL",
                 "T5500024037289120508304.C03.CEL",
                 "ave1",
                 "ave2",
                 "ave3",
                 "ave4",
                 "ave5")

# first drug expression set: treatment and controls
eset1 <- DEG[, eset1_names]
colnames(eset1)
# use
eset1

# differential expression analysis
library(limma)
# model.matrix
# without contrast
design <- model.matrix(~factor(rep(1:2, each = 5)))
colnames(design) <- c("treated", "control")
fit <- lmFit(eset1, design)
fit2 <- eBayes(fit)
t1 <- topTable(fit2, coef = 2, adjust = "fdr", number = 20)
t2 <- topTable(fit2, coef = 1, adjust = "fdr", number = 20)
t3 <- topTable(fit2, coef = 2, adjust = "none", number = 20)
t4 <- topTable(fit2, coef = 1, adjust = "none", number = 20)
t5 <- topTable(fit2, adjust = "fdr", number = 20)

# with the contrast
design2 <- model.matrix(~-1 + factor(rep(1:2, each = 5)))
colnames(design2) <- c("treated", "control")
contrast <- makeContrasts(treated - control, levels = design2)
fit <- lmFit(eset1, design2)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
t11 <- topTable(fit2)

# in the morning find out properly what the limma does and they mean?

j <- 1
for (i in 1:10) {
  j[i] = 10
}
j <- 1
j[2] <- 10
j[3] <- 10
j

# print ddao 10 times
for (i in 1:10) {
  print("This is ddao")
}

# heatmap after averaging - (selecting randomly the samples)
# Heatmap used to visualize differences between samples
# each row represents a gene and each column is an array 
# colored cell indicate expression levels of genes
# both samples and genes with similar expression profiles are clustered together
# Drawing heatmap in R uses a lot of memory and can take a long time
# typically data are filtered to include the genes which tell us the most
# about the biological variation

# OR alternatively, use a pre-selected list of genes of interets
# particular pathway, or gene-list from another study of interest

# heatmap after averaging
# using the DataFrame 'DEG'

# setwd
setwd("C:/Users/Ju Young/Desktop/University of Cape Town/MSc Comp Sci 2017/Research contents/Gene name annotation")
DEG <- readRDS("newDEG.rds")

# DataFrame information
dim(DEG)
rownames(DEG)
colnames(DEG)
str(DEG)
glimpse(DEG)

# Treatment names (selected based on the list_of_29_drugs.xls report)
# first rows selected

# 1. T5500024024213121906562.C03.CEL
# 2. T5500024024213121906562.H03.CEL
# 3. T5500024024213121906562.G04.CEL
# 4. T5500024024213121906562.F01.CEL
# 5. T5500024024213121906562.F05.CEL
# 6. T5500024024213121906562.C01.CEL
# 7. T5500024024213121906562.G01.CEL
# 8. T5500024024213121906562.B03.CEL
# 9. T5500024037497121008341.H04.CEL
# 10. T5500024024213121906562.D06.CEL
# 11. T5500024024213121906562.B01.CEL
# 12. T5500024035100021608460.A05.CEL
# 13. T5500024035100021608459.D06.CEL
# 14. TEC2003122315AA.CEL
# 15. TEC2004012203AA.CEL
# 16. T5500024024213121906562.E06.CEL
# 17. T5500024035100021608459.C04.CEL
# 18. T5500024024213121906562.F03.CEL
# 19. T5500024024213121906562.C05.CEL
# 20. T5500024037497121008341.E04.CEL
# 21. T5500024024213121906562.F06.CEL
# 22. T5500024024213121906562.C02.CEL
# 23. T5500024024213121906562.G05.CEL
# 24. T5500024037496121008326.F10.CEL
# 25. T5500024024213121906562.D01.CEL
# 26. T5500024024213121906562.C04.CEL
# 27. T5500024024213121906562.G03.CEL
# 28. T5500024024213121906562.E02.CEL
# 29. T5500024024213121906562.A04.CEL

# concatenate the treatment names
t_names <- c("T5500024024213121906562.C03.CEL",
             "T5500024024213121906562.H03.CEL",
             "T5500024024213121906562.G04.CEL",
             "T5500024024213121906562.F01.CEL",
             "T5500024024213121906562.F05.CEL",
             "T5500024024213121906562.C01.CEL",
             "T5500024024213121906562.G01.CEL",
             "T5500024024213121906562.B03.CEL",
             "T5500024037497121008341.H04.CEL",
             "T5500024024213121906562.D06.CEL",
             "T5500024024213121906562.B01.CEL",
             "T5500024035100021608460.A05.CEL",
             "T5500024035100021608459.D06.CEL",
             "TEC2003122315AA.CEL",
             "TEC2004012203AA.CEL",
             "T5500024024213121906562.E06.CEL",
             "T5500024035100021608459.C04.CEL",
             "T5500024024213121906562.F03.CEL",
             "T5500024024213121906562.C05.CEL",
             "T5500024037497121008341.E04.CEL",
             "T5500024024213121906562.F06.CEL",
             "T5500024024213121906562.C02.CEL",
             "T5500024024213121906562.G05.CEL",
             "T5500024037496121008326.F10.CEL",
             "T5500024024213121906562.D01.CEL",
             "T5500024024213121906562.C04.CEL",
             "T5500024024213121906562.G03.CEL",
             "T5500024024213121906562.E02.CEL",
             "T5500024024213121906562.A04.CEL")

# with the control
t_names <- c("T5500024024213121906562.C03.CEL",
             "T5500024024213121906562.H03.CEL",
             "T5500024024213121906562.G04.CEL",
             "T5500024024213121906562.F01.CEL",
             "T5500024024213121906562.F05.CEL",
             "T5500024024213121906562.C01.CEL",
             "T5500024024213121906562.G01.CEL",
             "T5500024024213121906562.B03.CEL",
             "T5500024037497121008341.H04.CEL",
             "T5500024024213121906562.D06.CEL",
             "T5500024024213121906562.B01.CEL",
             "T5500024035100021608460.A05.CEL",
             "T5500024035100021608459.D06.CEL",
             "TEC2003122315AA.CEL",
             "TEC2004012203AA.CEL",
             "T5500024024213121906562.E06.CEL",
             "T5500024035100021608459.C04.CEL",
             "T5500024024213121906562.F03.CEL",
             "T5500024024213121906562.C05.CEL",
             "T5500024037497121008341.E04.CEL",
             "T5500024024213121906562.F06.CEL",
             "T5500024024213121906562.C02.CEL",
             "T5500024024213121906562.G05.CEL",
             "T5500024037496121008326.F10.CEL",
             "T5500024024213121906562.D01.CEL",
             "T5500024024213121906562.C04.CEL",
             "T5500024024213121906562.G03.CEL",
             "T5500024024213121906562.E02.CEL",
             "T5500024024213121906562.A04.CEL",
             deg_col)

# extract t_names from DEG DataFrame
t_eset1 <- DEG[, t_names]
dim(t_eset1)
colnames(t_eset1)

# scale both rows and columns
x <- scale(t_eset1)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# FUT related GENES
# POFUT - polyexonic fucose
FUT <- rownames(DEG)
FUT1 <- grep("FUT", FUT, ignore.case = T)
FUT2 <- FUT[FUT1]
DEG_FUT <- DEG[FUT2,t_names]

# scale both rows and columns
x <- scale(DEG_FUT)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))

# -------------------------------------------------------------------------------

# SYIAL related genes i.e. STGAL family
STGAL <- rownames(DEG)
STGAL1 <- grep("ST", STGAL, ignore.case = T)
STGAL2 <- STGAL[STGAL1]
STGAL2 <- STGAL2[-c(1,16)]
DEG_STGAL <- DEG[STGAL2,t_names]

# scale both rows and columns
x <- scale(DEG_STGAL)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))

# -------------------------------------------------------------------------------

# BETA GALACTOSE FAMILY
BGAL <- rownames(DEG)
BGAL1 <- grep("B3|B4", BGAL, ignore.case = T)
BGAL2 <- BGAL[BGAL1]
BGAL2 <- BGAL2[-c(19,20)]
DEG_BGAL <- DEG[BGAL2,t_names]

# scale both rows and columns
x <- scale(DEG_BGAL)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))
# -------------------------------------------------------------------------------

# MGAT FAMILY
MGAT <- rownames(DEG)
MGAT1 <- grep("MGAT", MGAT, ignore.case = T)
MGAT2 <- MGAT[MGAT1]
DEG_MGAT <- DEG[MGAT2,t_names]

# scale both rows and columns
x <- scale(DEG_MGAT)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))
dev.off()
# -------------------------------------------------------------------------------

# ALL THE 'IMPORTANT' GENES
ALL <- c(BGAL2,MGAT2,FUT2,STGAL2)
DEG_ALL <- DEG[ALL,t_names]

# scale both rows and columns
x <- scale(DEG_ALL)
x <- t(scale(t(x)))

# draw heatmap
library(gplots)
heatmap.2(scale(x),
          dendrogram = "both", 
          trace = "none",
          col = greenred(100), 
          key = T,
          #main = "Heatmap",
          xlab = "Drugs",
          ylab = "GTs",
          cexRow = 0.3, cexCol = 0.2,
          margins = c(4,6),
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "ward.D2"))

# -------------------------------------------------------------------------------

