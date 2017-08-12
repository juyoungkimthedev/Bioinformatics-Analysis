##### preanalysis following the paper from Kihoon Cha et al 2014 #####
getwd()
wd <- "C:/Users/samsung/Desktop/University of Cape Town/2017 M.Sc/R SCRIPTS/preanalysis_1"
setwd(wd)
getwd()

# load cmap_instance.txt
cmap <- read.delim("cmap_instances_02.txt")
dim(cmap)
View(cmap)

# filter MCF7
library(dplyr)
cmap.mcf7 <- filter(cmap, cell2 == "MCF7")
dim(cmap.mcf7)
View(cmap.mcf7)

# filter HT_HG-U133A + HG-U133A
cmap.data <- filter(cmap.mcf7, array3 != "HT_HG-U133A_EA")
dim(cmap.data)
View(cmap.data)

# save in various forms
saveRDS(cmap.data, "./cmap.data.rds")
write.table(cmap.data, file = "./cmap.data.xls", quote = F, sep = "\t", col.names = NA)
dir()

# select only necessary columns
cmap.1 <- select(cmap.data, c(3,5,7,8,9,10))
dim(cmap.1)
View(cmap.1)

# replace unexpected ' in the perturbation scan id
class(cmap.1$vehicle_scan_id4) # factor
class(cmap.1$perturbation_scan_id) # factor
library(stringr)
cmap.1$perturbation_scan_id <- str_replace(cmap.1$perturbation_scan_id, "'", "")
View(cmap.1)

# repeated items?
# 1. test 
f <- data.frame(cmap_name = c("peter","peter","john","john","john","shiba"))
f$cmap_name <- as.character(f$cmap_name) # characterize
f$concentration <- c(2,2,3,3,4,5) # add another column for concentration
library(plyr)
ddply(f, .(cmap_name,concentration), nrow) # determines how many times each repeated
# 2. actual cmap
# first characterize
class(cmap.1$cmap_name)
class(cmap.1$concentration..M.)
cmap.1$cmap_name <- as.character(cmap.1$cmap_name)
colnames(cmap.1)[2] <- "concentration"
colnames(cmap.1)
cmap.2 <- select(cmap.1, c(1,2))
cmap.df <- as.data.frame(ddply(cmap.2, .(cmap_name, concentration), nrow))
View(cmap.df)
dim(cmap.df)

# excel report showing how many times the drug treated same concentration repeated
saveRDS(cmap.df, "./cmap.df.repeated.rds")
write.table(cmap.df, file = "./cmap.df.repeated.xls", quote = F, sep = "\t", col.names = NA)

# find out how many / which drugs repeated more than 5 times
colnames(cmap.df)[3] <- "v1"
cmap.df$v1 <- as.numeric(cmap.df$v1)
more.than.5 <- which(cmap.df$v1 >= 5)
cmap.df.repeat <- cmap.df[more.than.5,] # filtered drugs with biological replicates of more than 5
saveRDS(cmap.df.repeat, "./cmap.df.repeat.rds")
write.table(cmap.df.repeat, file = "./cmap.df.repeat.xls", quote = F, sep = "\t", col.names = NA)

### now filter out concentration (select highest concentration) 
### re-read the file into a data frame in R

# load the drug data and check against K.H. Cha et al 2014
# only difference is that KH Cha et al has metformin 
# and our study has monastrol
getwd()
setwd("C:/Users/samsung/Desktop/University of Cape Town/2017 M.Sc/R SCRIPTS/preanalysis_1/preanalysis_2")
cmap.drug.data <- read.delim("cmap.data.1.txt")
colnames(cmap.drug.data)[1] <- "index"
View(cmap.drug.data)

# [completed] 1. filter out the perturbation scan id for these drugs from cmap_instance file
# [completed] 2. randomly select 5 experiments from each drug to reduce bias by the number of experiments
# [completed] 3. gather information about perturbation vs vehicle id using functions
# 4. gather .CEL files for these items + normalize RMA + average the control values
# 5. construct matrix + normalize once again for different chip type 
# 6. annotate the genes
# 7. filter gt genes out / differential expression analysis

## cmap.data is hgu133a + hthgu133a + mcf7 files
## cmap.drug.data is 29 drug data information

##### pre filtering 2 for .CEL names #####
# filter the 29 cmap names
t1 <- cmap.data$cmap_name %in% cmap.drug.data$cmap_name
t2 <- cmap.data[t1,]
t3 <- t2[t2$concentration..M. %in% cmap.drug.data$concentration,]
length(unique(t3$cmap_name)) # returns 29

# filter the highest concentration in excel 
# randomly select 5 samples from each in excel
write.table(t3, file = "./t3.xls", quote = F, sep = "\t", col.names = NA)

# re-load final data with 29 drugs and 5 samples each
# perturbation.data.xls + perturbation.data.txt
perturbation <- read.delim("perturbation.data.txt")
dim(perturbation)
length(unique(perturbation$cmap_name))

# perturbation + vehicle scan id separation
colnames(perturbation)[1] <- "index"
# remove unexpected ' in the perturbation scan id
library(stringr)
perturbation$perturbation_scan_id <- str_replace(perturbation$perturbation_scan_id, "'", "")
View(perturbation)
dim(perturbation)
saveRDS(perturbation, "./perturbation.rds")
write.table(perturbation, file = "./perturbation.xls", quote = F, sep = "\t", col.names = NA)

# separate the perturbation + vehicle scan numbers
list.perturbation <- split(perturbation[, c("perturbation_scan_id","vehicle_scan_id4")],
                           as.character(perturbation$perturbation_scan_id))

# write functions that separates the name
# part of the function was adopted from the longevityTools project
sampleList1 <- function(x) {
  c <- unlist(strsplit(as.character(x[, 2]), "\\."))
  if (max(nchar(c) > 4)) {
    c <- gsub("'", "", as.character(x[, 2]))
  } else {
    c <- c[nchar(c) != 0]
    # the following code pastes the perturbation scan id with the vehicle scan tag i.e. H01 etc
    # gsub function replaces all matches of a string
    c <- paste0(gsub("\\..*", "", as.character(x[, 1])), ".", c)
    # the following code replaces unexpected "'" in the perturbation scan id
    c <- gsub("'", "", c)
  }
  t <- gsub("'", "", as.character(x[, 1]))
  return(list(list(t = t, c = c)))
}

# sapply
apply.perturbation <- sapply(list.perturbation, sampleList1)
class(apply.perturbation) # a list
apply.perturbation # shows the information of drug-treated vs control samples

# sapply application 
# show only control vehicle .CEL names
apply.1.perturbation <- sapply(apply.perturbation, function(x) x[2]) # greps the vehicle id 
vehicle.df <- as.data.frame(unlist(apply.1.perturbation))
saveRDS(vehicle.df, "./vehicle.df.rds")
write.table(vehicle.df, file = "./vehicle.df.xls", quote = F, sep = "\t", col.names = NA)
# the above excel sheet needs manual curation to separate .CEL names
# continue from step 4 above

# use (modified)List of 29 drug treated samples with 5 replicates.xls to work further