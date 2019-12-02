#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.train.file <- as.character(parameters[2])
meta.vali.file <- as.character(parameters[3])
boruta.res.file <- as.character(parameters[4])
myfactors <- as.character(parameters[5])
mylevels <- as.character(parameters[6])
study <- as.character(parameters[7])
mynorm <- as.logical(as.character(parameters[8]))
threshold <- as.numeric(as.character(parameters[9]))
doROCfilt <- as.logical(as.character(parameters[10]))

assess.RF <- function(pred.test, prob.cut.positive=0.5, positive.level="YES", ALL=NULL, ALL.classes){
	true.test <- ALL.classes
	if(!is.null(ALL)) names(true.test) <- colnames(ALL) else names(true.test) <- names(ALL.classes)
	true.test <- true.test[names(pred.test)]

	pred.test.positive <- names(pred.test)[which(pred.test>prob.cut.positive)]
	pred.test.negative <- names(pred.test)[which(pred.test<=prob.cut.positive)]
	true.test.positive <- names(true.test)[which(true.test==positive.level)]
	true.test.negative <- names(true.test)[which(true.test!=positive.level)]
	VP <- length(intersect(pred.test.positive, true.test.positive))
	FP <- length(intersect(pred.test.positive, true.test.negative))
	VN <- length(intersect(pred.test.negative, true.test.negative))
	FN <- length(intersect(pred.test.negative, true.test.positive))
	res <- c(VP/(VP+FN), FP/(FP+VN), VP, FP, VN, FN)
	names(res) <- c("VPR", "FPR", "VP", "FP", "VN", "FN")
	return(res)
}

mydir <- system("ls -d prediction.extra.tables", intern=TRUE)
if(length(mydir)!=1) system("mkdir prediction.extra.tables")

# training data
dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat) <- dat[, 1]
rownames(dat) <- gsub("-", ".", rownames(dat))
meta <- read.table(file=meta.train.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
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
yes.group <- levels(as.factor(GROUPS))[1]
no.group <- levels(as.factor(GROUPS))[2]
meta <- meta[ind, ]

# validation data
dat.new <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat.new) <- dat.new[, 1]
rownames(dat.new) <- gsub("-", ".", rownames(dat.new))
meta.new <- read.table(file=meta.vali.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta.new) <- meta.new[, 1]
mysamples.new <- intersect(meta.new[, 1], colnames(dat.new))
dat.new <- as.matrix(dat.new[, mysamples.new])
meta.new <- meta.new[mysamples.new, ]

if(mynorm==TRUE){
	dat.all <- prop.table(as.matrix(dat), margin=2)*100
	dat.all.new <- prop.table(as.matrix(dat.new), margin=2)*100
} else {
	dat.all <- dat; dat.all.new <- dat.new
}

# take only confirmed boruta features
mytb <- read.table(file=boruta.res.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
mytb$feature <- gsub("-", ".", mytb$feature)
mytb <- mytb[which(mytb$decision=="Confirmed"), ]

numvars.v <- c()
AUC.v <- c()
pct.of.CA.for.YES <- c()
pct.of.CA.for.NO <- c()
abs.of.CA.for.YES <- c()
abs.of.CA.for.NO <- c()
tot.yes.v <- c()
tot.no.v <- c()

if(doROCfilt==TRUE){
	file.pdf <- paste("ROC.RF.for.threshold", threshold, study, "pdf", sep=".")
} else {
	file.pdf <- paste("ROC.RF", study, "pdf", sep=".")
}

pdf(file=file.pdf)

# We see if we can obtain better results by removing step by step from the model the less important feature until reach the top 3
for(numvars in nrow(mytb):3){

myvars <- intersect(mytb$feature[1:numvars], rownames(dat))
study.tem <- paste(study, length(myvars), sep=".")

dat <- dat.all[myvars, ]
dat.new <- dat.all.new[myvars, ]

library(randomForest)
library(pracma)

my.mtry <- round(sqrt(nrow(dat)-1))

training <- as.data.frame(t(dat))
training$y <- as.factor(GROUPS)
testing <- as.data.frame(t(dat.new))

set.seed(777)

# Fit the random forest model on the training samples
RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=1000, keep.forest=TRUE)

# Predict the status for the validation samples using fitted model
mypred <- predict(RF, type="prob", newdata=testing)

myclass <- apply(mypred, MARGIN=1, FUN=function(x){
		ind <- which(x==max(x))
		if(length(ind)>1) return(NA) else return(colnames(mypred)[which(x==max(x))])
	})

res <- as.data.frame(mypred)
res$class.predicted <- myclass
res$sample <- rownames(res)
res <- res[, c("sample", setdiff(colnames(res), "sample"))]

if(myfactor2!="") GROUPS.new <- paste(meta.new[, myfactor1], meta.new[, myfactor2], sep=".") else GROUPS.new <- as.character(meta.new[, myfactor1])
names(GROUPS.new) <- rownames(meta.new)

res$true.class <- GROUPS.new[rownames(res)]

write.table(res, file=file.path("./prediction.extra.tables", paste(study.tem, "prediction", "tsv", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)

ind <- which(!is.na(res$true.class))
res <- res[ind, ]; mypred <- mypred[ind, ]

# Filter samples based on threshold and compute pct.of.CA.for (filtered samples are not considered for the total sum of a group)
ind <- which(res[, 2]>=threshold|res[, 3]>=threshold)

if(length(ind)>0){
res.filt <- res[ind, ]; mypred.filt <- mypred[ind, ]
n.yes <- length(which(res.filt$class.predicted==yes.group&res.filt$true.class==yes.group))
n.no <- length(which(res.filt$class.predicted!=yes.group&res.filt$true.class!=yes.group))
tot.yes <- length(which(res.filt$true.class==yes.group))
tot.no <- length(which(res.filt$true.class!=yes.group))
pct.of.CA.for.YES <- c(pct.of.CA.for.YES, 100*n.yes/tot.yes)
pct.of.CA.for.NO <- c(pct.of.CA.for.NO, 100*n.no/tot.no)
abs.of.CA.for.YES <- c(abs.of.CA.for.YES, n.yes)
abs.of.CA.for.NO <- c(abs.of.CA.for.NO, n.no)
tot.yes.v <- c(tot.yes.v, tot.yes)
tot.no.v <- c(tot.no.v, tot.no)
} else {
pct.of.CA.for.YES <- c(pct.of.CA.for.YES, NA)
pct.of.CA.for.NO <- c(pct.of.CA.for.NO, NA)
abs.of.CA.for.YES <- c(abs.of.CA.for.YES, NA)
abs.of.CA.for.NO <- c(abs.of.CA.for.NO, NA)
tot.yes.v <- c(tot.yes.v, NA)
tot.no.v <- c(tot.no.v, NA)
}

if(doROCfilt==TRUE){ res <- res.filt; mypred <- mypred.filt }

# PLot ROC curve of the prediction
ALL.classes <- res$true.class
ALL.classes <- as.factor(ALL.classes)
names(ALL.classes) <- rownames(res)
# Discrimination thresholds for the computation of the ROC curve
cuts <- seq(from=0, to=1, length.out=50)
cuts <- c(cuts, 0.5); cuts <- cuts[order(cuts, decreasing=FALSE)]
TPR.v <- FPR.v <- TP.v <- FP.v <- TN.v <- FN.v <- rep(NA, length(cuts))
names(TPR.v) <- names(FPR.v) <- names(TP.v) <- names(FP.v) <- names(TN.v) <- names(FN.v) <- as.character(cuts)
# Assess prediction for each cut
for(mycut in cuts){
	res <- assess.RF(pred.test=mypred[names(ALL.classes), yes.group], prob.cut.positive=mycut,
			positive.level=yes.group, ALL.classes=ALL.classes)
	
	if(!is.na(res["VPR"])&!is.na(res["FPR"])&res["VPR"]!=Inf&res["FPR"]!=Inf){
		TPR.v[as.character(mycut)] <- res["VPR"]
		FPR.v[as.character(mycut)] <- res["FPR"]
		TP.v[as.character(mycut)] <- res["VP"]
		FP.v[as.character(mycut)] <- res["FP"]
		TN.v[as.character(mycut)] <- res["VN"]
		FN.v[as.character(mycut)] <- res["FN"]
	}
}

cuts <- format(cuts, digits=1)

ind <- which(!is.na(TPR.v))
TPR.v <- TPR.v[ind]; FPR.v <- FPR.v[ind]; TP.v <- TP.v[ind]; FP.v <- FP.v[ind]; TN.v <- TN.v[ind]; FN.v <- FN.v[ind]; cuts <- cuts[ind]

AUC <- abs(trapz(FPR.v, TPR.v))

if(is.nan(AUC)){
	warning("AUC is NaN")
	return()
} else {
	AUC <- format(AUC, digits=2)
}

mymain <- paste("AUC:", AUC, "numvars:", numvars, sep=" ")

plot(x=FPR.v, y=TPR.v, type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), main=mymain)
abline(a=0, b=1, col="red")
ind <- c(which(cuts==0.53), which(cuts==0.63), which(cuts==0.73))
if(length(ind)>0) text(x=FPR.v[ind], y=TPR.v[ind], labels=cuts[ind], col="blue", cex=0.6)
lines(x=FPR.v, y=TPR.v, col="blue")

points.df <- data.frame(FPR=FPR.v, TPR=TPR.v, TP=TP.v, FP=FP.v, TN=TN.v, FN=FN.v, FNR=1-TPR.v, discrimination.threshold=cuts)
write.table(points.df, file=file.path("./prediction.extra.tables", paste(study.tem, "prediction.points.ROC.for.threshold", threshold, "tsv", sep=".")),
		sep="\t", quote=FALSE, row.names=FALSE)

numvars.v <- c(numvars.v, numvars)
AUC.v <- c(AUC.v, AUC)

}
dev.off()

pct.true.pos.for.thr <- data.frame(numvars=numvars.v, abs.of.CA.for.YES=abs.of.CA.for.YES, abs.of.CA.for.NO=abs.of.CA.for.NO,
						tot.yes.v=tot.yes.v, tot.no.v=tot.no.v,
						pct.of.CA.for.YES=pct.of.CA.for.YES, pct.of.CA.for.NO=pct.of.CA.for.NO,
						tot.pct.of.CA=(abs.of.CA.for.YES+abs.of.CA.for.NO)/(tot.yes.v+tot.no.v),
						AUC=AUC.v)

colnames(pct.true.pos.for.thr) <- c("numvars", paste("num.of.CA.for", yes.group, sep="."), paste("num.of.CA.for", no.group, sep="."),
								paste("total", yes.group, sep="."), paste("total", no.group, sep="."),
								paste("pct.of.CA.for", yes.group, sep="."), paste("pct.of.CA.for", no.group, sep="."),
								"total.pct.of.CA", "AUC")

write.table(pct.true.pos.for.thr, file=paste("pct.true.pos.for.threshold", threshold, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)


if(doROCfilt==FALSE){
write.table(data.frame(numvars=numvars.v, AUC=AUC.v),
			file=paste("ROC.RF", study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
} else {
write.table(data.frame(numvars=numvars.v, AUC=AUC.v),
			file=paste("ROC.RF.for.threshold", threshold, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}







