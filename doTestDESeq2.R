#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactors <- as.character(parameters[3])
mylevels <- as.character(parameters[4])
study <- as.character(parameters[5])
doFilt <- as.logical(as.character(parameters[6]))
descr.column <- as.character(parameters[7])
mypaired <- as.logical(as.character(parameters[8]))
mynorm <- as.logical(as.character(parameters[9]))
myfitType <- as.character(parameters[10])

source("myfunctions.downstream.R")
library(DESeq2)

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

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

# This leaves data unchanged if matrix was in absolute counts already, otherwise data turns into absolute counts
for(i in 1:ncol(dat)) dat[, i] <- round(dat[, i]/min(dat[which(dat[, i]!=0), i]))

doCCA <- FALSE

if(mynorm==TRUE) datnorm <- prop.table(as.matrix(dat), margin=2)*100 else datnorm <- dat

# Filter features by low signal
# A minimum of 3 samples in at least one group has to be above cutoff
if(doFilt==TRUE){
	mycut <- (100/nrow(dat))/5
	good <- doFiltering(datnorm, cutoff=mycut, group=GROUPS)
	dat <- dat[good, ]
	datnorm <- datnorm[good, ]
	if(mynorm==TRUE) mydata <- as.data.frame(datnorm) else mydata <- as.data.frame(dat)
	mydata$ID <- rownames(mydata)
	mydata <- mydata[, c("ID", colnames(dat))]
	write.table(mydata, file=paste("filtered.freqs", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}

# CCA and adonis test for more than 2 groups will be performed (this option was not used for this study)
if(length(mylevels)>2){
	if(doCCA==TRUE) correspondence.analysis(MATRIX=datnorm,
					GROUPS=GROUPS, FILE=paste("adonis.CCA", "all.groups", study, "pdf", sep="."))
}

mycombn <- combn(mylevels,2)

for(i in 1:ncol(mycombn)){

	g1 <- mycombn[1,i]
    	g2 <- mycombn[2,i]
	i1 <- which(GROUPS==g1)
	i2 <- which(GROUPS==g2)

	if(mypaired==TRUE){
		names(i1) <- meta$donor[i1]
		names(i2) <- meta$donor[i2]
		mydonors <- intersect(names(i1), names(i2))
		i1 <- i1[mydonors]; i2 <- i2[mydonors]
	}

	ind <- c(i1, i2)

	tag <- paste(g1, "vs", g2, sep=".")
	
	# CCA and adonis test for 2 groups (this option was not used for this study)
	if(doCCA==TRUE) correspondence.analysis(MATRIX=datnorm[, ind],
					GROUPS=GROUPS[ind], FILE=paste("adonis.CCA", tag, study, "pdf", sep="."))
	
	# deseq2
	mytb <- differential.tests(MATRIX=dat[, ind], GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="DESeq2.test",
						PAIRED=mypaired, CUTOFF=0.1, COLNAME="feature", FILE=paste("DESeq2.test", tag, study, sep="."),
						DESCRIPTION=DESCRIPTION, myfitType=myfitType, meta=meta[ind, ])

	myfeatures <- as.character(mytb[which(mytb$p.value<=0.05), "feature"])
	tem <- as.data.frame(dat[myfeatures, ind]); tem$ID <- rownames(tem); tem <- tem[, c("ID", setdiff(colnames(tem), "ID"))]
	write.table(tem, file=paste("significant.rel.freqs", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE) 
	
	if(mypaired==TRUE){
		pdf(file=paste("DESeq.test.paired.with.lines", tag, study, "pdf", sep="."))
		for(myfeature in myfeatures){
			plotPaired(dat=datnorm[, ind], feature=myfeature, group=GROUPS[ind],
					group.x=g1, group.y=g2, individual.factor=meta$donor[ind], DESCRIPTION=DESCRIPTION)
		}
		dev.off()
	}

	print(tag)
}




