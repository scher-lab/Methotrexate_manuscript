suppressPackageStartupMessages(library(labdsv))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(beeswarm))
suppressPackageStartupMessages(library(vegan))

# -------------------------------------------------------------
#
# MATRIX:
#   Abundances matrix where the columns represent samples (observations)
#   and the rows represent features to study (taxons, genes, ...).
#
# GROUPS:
#   Group assigned to each sample.
#   Named vector where the values represent distinct groups/conditions and 
#   the names represent the samples.
#
# COLOURS:
#   Colour assigned to each sample.
#   Named vector where the values represent distinct colours and 
#   the names represent the samples.
#
# LEGENDS:
#   Legend/Label assigned to each colour.  
#   Named vector where the values represent distinct labels and 
#   the names represent the colours.
#
#
# Colours are codified by name ("red", "blue", "green", ...).
# Usually colours matches to groups.
#
# -------------------------------------------------------------


doFiltering <- function(dat, cutoff, group){
	good <- c()
	for(mygroup in unique(as.character(group))){
		dat.temp <- dat[, which(group==mygroup)]
		cond <- apply(dat.temp, MARGIN=1, FUN=function(x) if(length(which(x>=cutoff))>=3) return(TRUE) else return(FALSE))
		good <- c(good, rownames(dat.temp)[which(cond==TRUE)])
	}
	return(unique(good))
}

#################################################
# Test differences between groups of observations (columns)
# for each feature (row) in the matrix.
#
# TEST: {t.test, wilcox.test, DESeq2.test}. For this particular study we just applied DESeq2.test
#################################################
differential.tests <- function(MATRIX, GROUPS, GROUP1, GROUP2, TEST, PAIRED, CUTOFF, COLNAME, FILE="", DESCRIPTION=NULL,
						myfitType="parametric", meta=NULL)
{

	if(length(apropos("mynorm"))==0) mynorm <- TRUE
  # Observations in each group
  items.group1 <- names(GROUPS)[GROUPS==GROUP1]
  items.group2 <- names(GROUPS)[GROUPS==GROUP2]
  
  # Differences between both groups
  if(TEST=="t.test"|TEST=="wilcox.test"){
  p.values <- apply(MATRIX, MARGIN=1, function(values) 
  {
    do.call(TEST, list(x=values[items.group1], y=values[items.group2], paired=PAIRED))$p.value
  })
  }
  
	if(TEST=="DESeq.test"){
		countTable <- MATRIX
		if(any(countTable==0)){
			countTable <- countTable + 1
		}
		condition <- as.factor(GROUPS)
		cds <- newCountDataSet(countTable, condition)
		cds <- estimateSizeFactors(cds)
		if(PAIRED==FALSE){
			cds <- estimateDispersions(cds, fitType=myfitType)
			res <- nbinomTest(cds, GROUP1, GROUP2)
			p.values <- res$pval
		} else {
			cds <- estimateDispersions(cds, method="blind", fitType=myfitType)
			vsd <- getVarianceStabilizedData(cds)
			if(mynorm==TRUE) vsd <- prop.table(vsd, margin=2)*100
			p.values <- apply(vsd, MARGIN=1, function(values) 
  			{
    				do.call("wilcox.test", list(x=values[items.group1], y=values[items.group2], paired=PAIRED))$p.value
  			})
		}
	}

	if(TEST=="DESeq2.test"){
		if(PAIRED==FALSE){
			coldata <- data.frame(sample=names(GROUPS), condition=as.factor(GROUPS))
			rownames(coldata) <- names(GROUPS)
			dds <- DESeqDataSetFromMatrix(countData = MATRIX, colData = coldata, design = ~ condition)
		} else {
			coldata <- data.frame(sample=names(GROUPS), condition=as.factor(GROUPS), donor=as.factor(meta$donor))
			rownames(coldata) <- names(GROUPS)
			dds <- DESeqDataSetFromMatrix(countData = MATRIX, colData = coldata, design = ~ donor + condition)
		}
		dds <- DESeq(dds, fitType=myfitType)
		res <- lfcShrink(dds, coef=2, type="apeglm")
		res <- as.data.frame(res@listData)
		p.values <- res$pvalue
		p.values.adj <- res$padj
	}

  if(mynorm==TRUE) MATRIX <- prop.table(MATRIX, margin=2)*100
  if(TEST!="DESeq2.test"){
	p.values.adj <- p.adjust(p.values, method="fdr")
  }
  
  # Mean value for each group
  means.group1 <- apply(MATRIX, MARGIN=1, function(values) {mean(values[items.group1])})
  means.group2 <- apply(MATRIX, MARGIN=1, function(values) {mean(values[items.group2])})

  mymin1 <- min(means.group1[which(means.group1!=0)])
  mymin2 <- min(means.group2[which(means.group2!=0)])
  mymin <- mean(c(mymin1, mymin2))
  FC <- log2((means.group1+mymin)/(means.group2+mymin))

  df <- data.frame(p.value=p.values, adj.p.value=p.values.adj, 
                   mean1=means.group1, mean2=means.group2, log2FC=FC, stringsAsFactors=FALSE)

  # Lower p.value --> higher difference
  df <- df[order(df$p.value, decreasing=FALSE),]
  
  # Get only significant differences
  df2plot <- df[which(df$adj.p.value < CUTOFF),]
  if(nrow(df2plot)<5) df2plot <- df[which(df$p.value < CUTOFF),]
  
  # Write p.values into a text file
  df.aux <- cbind(rownames(df), df)
  colnames(df.aux)[1] <- COLNAME
  if(!is.null(DESCRIPTION)) df.aux$DESCRIPTION <- DESCRIPTION[as.character(df.aux[, 1])]

  if(FILE!="") write.table(df.aux, paste(FILE,"tsv",sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Boxplot values of rows with significant differences
  if(FILE!=""){

  pdf(file=paste(FILE,"pdf",sep="."))
  for(sign.row in rownames(df2plot))
  {
    values <- list(MATRIX[sign.row,items.group1], MATRIX[sign.row,items.group2])
    names(values) <- c(GROUP1, GROUP2)
    if(!is.null(DESCRIPTION)) mymain <- DESCRIPTION[sign.row] else mymain <- ""
    beeswarm(values, pch=16, col=c("red","blue"), las=1, xlab="Condition", ylab=sign.row, main=mymain, cex.lab=0.6, cex.main=0.5)
    boxplot (values, las=1,  add=TRUE)
    legend("topright", legend=c(paste("p.value:"         , format(df2plot[sign.row,"p.value"    ],digits=2)),
                                paste("adjusted p.value:", format(df2plot[sign.row,"adj.p.value"],digits=2))), cex=0.7)
  }
  temp <- dev.off()
  
  }

  return(df.aux)
}

