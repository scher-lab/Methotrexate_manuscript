#!/usr/bin/env Rscript
###############
# I have increased some parameters that guarantee some reproducibility:
#  - Number of boruta runs: the more runs, the less likelihood for a feature to remain tentative (tentative means that
#	there is not enough confidence to rejected or confirmed it).
#  - Number of trees for the random forests: the more trees, the more stability for the importance score associated with each feature.

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
p.vals.file <- as.character(parameters[3])
myfactors <- as.character(parameters[4])
mylevels <- as.character(parameters[5])
study <- as.character(parameters[6])
descr.column <- as.character(parameters[7])
mynorm <- as.logical(as.character(parameters[8]))
pval.cutoff <- as.numeric(as.character(parameters[9]))

source("myfunctions.downstream.R")

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat) <- dat[, 1]
rownames(dat) <- gsub("-", ".", rownames(dat))
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat)

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples

mylevels <- unlist(strsplit(mylevels, split=","))

mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]

ind <- which(is.element(GROUPS, mylevels))

dat <- dat[, ind]
GROUPS <- GROUPS[ind]
meta <- meta[ind, ]

if(mynorm==TRUE) dat <- prop.table(as.matrix(dat), margin=2)*100

mytb <- read.table(file=p.vals.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
mytb$feature <- gsub("-", ".", mytb$feature)
rownames(mytb) <- mytb$feature

# We remove features with desq2 adjusted p-value equal to NA
myfeatures <- intersect(mytb$feature[which(mytb$p.value<=pval.cutoff&!is.na(mytb$adj.p.value))], rownames(dat))

mytb <- mytb[myfeatures, ]
dat <- dat[myfeatures, ]

#######
## Boruta feature selection algorithm via random forest  

library(Boruta)

df <- as.data.frame(t(dat))
df$group <- as.factor(GROUPS)

# Run Boruta Algorithm
set.seed(777)
boruta <- Boruta(group ~ ., data = df, doTrace = 2, maxRuns = 500, ntree = 1000)

# In some circumstances (too short Boruta run, unfortunate mixing of shadow attributes, tricky dataset...), Boruta can leave some
# attributes Tentative. ‘TentativeRoughFix’ performs a simplified, weaker test for judging such attributes.
# This function claims as Confirmed those attributes that have median importance higher than the median importance of maximal shadow
# attribute, and the rest as Rejected. Depending of the user choice, medians for the test are count over last round, all rounds or N
# last importance source runs.
boruta <- TentativeRoughFix(boruta)

# Boruta importance statistics
res <- attStats(boruta)
rownames(res) <- gsub("[`]", "", rownames(res))
res$feature <- rownames(res)
res <- merge(res, mytb, by="feature")

res <- res[order(res$medianImp, decreasing=TRUE), ]
res <- rbind(res[which(res$decision=="Confirmed"), ], res[which(res$decision!="Confirmed"), ])

myfeatures <- res$feature[which(res$decision=="Confirmed")]
write(myfeatures, file=paste("RF.boruta", length(myfeatures), "important.variables", study, "txt", sep="."), sep="\n")

write.table(res, file=paste("results.boruta", study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)



