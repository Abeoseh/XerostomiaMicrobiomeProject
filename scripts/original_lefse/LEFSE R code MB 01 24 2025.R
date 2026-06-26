library(lefser)
library(dplyr)
library('SummarizedExperiment')
countdata=read.csv('Xero18v24_counts3.csv',row.names=1,sep=',')
coldata=read.csv('Xero18v24_ColData.csv',row.names='Sample',sep=',')

radata=as.matrix(countdata)
se=SummarizedExperiment(assays=list(counts=radata),colData=as.matrix(coldata))
se

#determine differential features of BL vs. 18m with site as sub-class
lefser_resultsTP=lefser(se,classCol='dex',subclassCol='Site',trim.names=FALSE,checkAbundances=FALSE, lda.threshold=0.001,kruskal.threshold=0.05,wilcox.threshold=0.05)
View(lefser_resultsTP)
lefserPlot(lefser_resultsTP,trim.names=FALSE)

lefser_resultsSite=lefser(se,groupCol='Site',trim.names=FALSE,checkAbundances=FALSE, lda.threshold=0,kruskal.threshold=0.05,wilcox.threshold=0.05)
View(lefser_resultsSite)
lefserPlot(lefser_resultsSite)

#data must be bimodal
lefser_resultsXero=lefser(se,groupCol='Xero',trim.names=FALSE,checkAbundances=FALSE, lda.threshold=0,kruskal.threshold=0.05,wilcox.threshold=0.05)
View(lefser_resultsXero)
lefserPlot(lefser_resultsXero,trim.names=FALSE)
View(lefser_resultsXero$scores[[]])

#determine differntial features of BL vs. 18m 2's vs. 3s with site as sub-class
countdata2v3=read.csv('Xero18v24 2s vs 3s counts.csv',row.names=1,sep=',')
coldata2v3=read.csv('Xero18v24 2s vs 3s ColData.csv',row.names='Sample',sep=',')

radata2v3=as.matrix(countdata2v3)
se2v3=SummarizedExperiment(assays=list(counts=radata2v3),colData=as.matrix(coldata2v3))
se2v3

lefser_results2v3=lefser(se2v3,groupCol='Xero',blockCol='Site',trim.names=FALSE,checkAbundances=FALSE, lda.threshold=0.001,kruskal.threshold=0.05,wilcox.threshold=0.05)
View(lefser_results2v3)
lefserPlot(lefser_results2v3,trim.names=FALSE)

#determine differntial features of BL vs. 18m 2's vs. 4's with site as sub-class
countdata2v4=read.csv('Xero18v24 2s vs 4s counts.csv',row.names=1,sep=',')
coldata2v4=read.csv('Xero18v24 2s vs 4s ColData.csv',row.names='Sample',sep=',')

radata2v4=as.matrix(countdata2v4)
se2v4=SummarizedExperiment(assays=list(counts=radata2v4),colData=as.matrix(coldata2v4))
se2v4

lefser_results2v4=lefser(se2v4,groupCol='Xero',blockCol='Site',trim.names=FALSE,checkAbundances=FALSE, lda.threshold=0.001,kruskal.threshold=0.05,wilcox.threshold=0.05)
View(lefser_results2v4)
lefserPlot(lefser_results2v4,trim.names=FALSE)')