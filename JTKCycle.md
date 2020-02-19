install.packages("/Users/lopeze/Downloads/MetaCycle", repos = NULL, type="source")

========
Prepare data file for wt53 analysis
========

```R
require(DESeq2)

colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

#Eliminate Frq08_DD24_rep3 sample from colData and countData
colData <- colData[!(colData$Sample=="Frq08_DD24_rep1"),]
countData <- subset(countData, select=-Frq08_DD24_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep2"),]
countData <- subset(countData, select=-Frq08_DD24_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]
countData <- subset(countData, select=-Frq08_DD24_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep1"),]
countData <- subset(countData, select=-Frq08_DD18_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep2"),]
countData <- subset(countData, select=-Frq08_DD18_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep3"),]
countData <- subset(countData, select=-Frq08_DD18_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep1"),]
countData <- subset(countData, select=-Frq08_DD12_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep2"),]
countData <- subset(countData, select=-Frq08_DD12_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep3"),]
countData <- subset(countData, select=-Frq08_DD12_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep1"),]
countData <- subset(countData, select=-Frq08_DD6_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep2"),]
countData <- subset(countData, select=-Frq08_DD6_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep3"),]
countData <- subset(countData, select=-Frq08_DD6_rep3)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep1"),]
countData <- subset(countData, select=-Frq08_LL6_rep1)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep2"),]
countData <- subset(countData, select=-Frq08_LL6_rep2)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep3"),]
countData <- subset(countData, select=-Frq08_LL6_rep3)



colData <- colData[!(colData$Sample=="Wc153_DD6_rep1"),]
countData <- subset(countData, select=-Wc153_DD6_rep1)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep2"),]
countData <- subset(countData, select=-Wc153_DD6_rep2)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep3"),]
countData <- subset(countData, select=-Wc153_DD6_rep3)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep1"),]
countData <- subset(countData, select=-Wc153_LL6_rep1)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep2"),]
countData <- subset(countData, select=-Wc153_LL6_rep2)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep3"),]
countData <- subset(countData, select=-Wc153_LL6_rep3)

colData <- colData[!(colData$Sample=="WT53_LL6_rep1"),]
countData <- subset(countData, select=-WT53_LL6_rep1)
colData <- colData[!(colData$Sample=="WT53_LL6_rep2"),]
countData <- subset(countData, select=-WT53_LL6_rep2)
colData <- colData[!(colData$Sample=="WT53_LL6_rep3"),]
countData <- subset(countData, select=-WT53_LL6_rep3)

write.table(countData,"countData_WT",sep="\t",na="",quote=F)
```
```R
===============
require(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData_WT",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)


#To extract the first column of a file:
as.data.frame( countData[,1], drop=false)

#To extract the rownames of the countData file:
rownames(countData)
write.table(rownames(countData,"countData_annot"))

#Do it in Bash
cut -f1 countData_WT > counData_annot
```


=============
Run programme
=============
```R
source("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/countData_annot")
data <- read.delim("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/countData_WT")

#rownames(data) <- data[,1]
#data <- data[,-1]
jtkdist(4, 3)       # 13 total time points, 2 replicates per time point

periods <- 4      # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,6)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/JTK_WT53_2.txt"),row.names=F,col.names=T,quote=F,sep="\t")

```
=============
Prepare data file for Frq KO analysis   
=============
```R
require(DESeq2)

colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

#Eliminate Frq08_DD24_rep3 sample from colData and countData
colData <- colData[!(colData$Sample=="WT53_DD24_rep1"),]
countData <- subset(countData, select=-WT53_DD24_rep1)
colData <- colData[!(colData$Sample=="WT53_DD24_rep2"),]
countData <- subset(countData, select=-WT53_DD24_rep2)
colData <- colData[!(colData$Sample=="WT53_DD24_rep3"),]
countData <- subset(countData, select=-WT53_DD24_rep3)
colData <- colData[!(colData$Sample=="WT53_DD18_rep1"),]
countData <- subset(countData, select=-WT53_DD18_rep1)
colData <- colData[!(colData$Sample=="WT53_DD18_rep2"),]
countData <- subset(countData, select=-WT53_DD18_rep2)
colData <- colData[!(colData$Sample=="WT53_DD18_rep3"),]
countData <- subset(countData, select=-WT53_DD18_rep3)
colData <- colData[!(colData$Sample=="WT53_DD12_rep1"),]
countData <- subset(countData, select=-WT53_DD12_rep1)
colData <- colData[!(colData$Sample=="WT53_DD12_rep2"),]
countData <- subset(countData, select=-WT53_DD12_rep2)
colData <- colData[!(colData$Sample=="WT53_DD12_rep3"),]
countData <- subset(countData, select=-WT53_DD12_rep3)
colData <- colData[!(colData$Sample=="WT53_DD6_rep1"),]
countData <- subset(countData, select=-WT53_DD6_rep1)
colData <- colData[!(colData$Sample=="WT53_DD6_rep2"),]
countData <- subset(countData, select=-WT53_DD6_rep2)
colData <- colData[!(colData$Sample=="WT53_DD6_rep3"),]
countData <- subset(countData, select=-WT53_DD6_rep3)
colData <- colData[!(colData$Sample=="WT53_LL6_rep1"),]
countData <- subset(countData, select=-WT53_LL6_rep1)
colData <- colData[!(colData$Sample=="WT53_LL6_rep2"),]
countData <- subset(countData, select=-WT53_LL6_rep2)
colData <- colData[!(colData$Sample=="WT53_LL6_rep3"),]
countData <- subset(countData, select=-WT53_LL6_rep3)



colData <- colData[!(colData$Sample=="Wc153_DD6_rep1"),]
countData <- subset(countData, select=-Wc153_DD6_rep1)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep2"),]
countData <- subset(countData, select=-Wc153_DD6_rep2)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep3"),]
countData <- subset(countData, select=-Wc153_DD6_rep3)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep1"),]
countData <- subset(countData, select=-Wc153_LL6_rep1)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep2"),]
countData <- subset(countData, select=-Wc153_LL6_rep2)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep3"),]
countData <- subset(countData, select=-Wc153_LL6_rep3)

colData <- colData[!(colData$Sample=="Frq08_LL6_rep1"),]
countData <- subset(countData, select=-Frq08_LL6_rep1)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep2"),]
countData <- subset(countData, select=-Frq08_LL6_rep2)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep3"),]
countData <- subset(countData, select=-Frq08_LL6_rep3)

write.table(countData,"countData_FrqKO",sep="\t",na="",quote=F)
```
```R
===============
require(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData_WT",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)


#To extract the first column of a file:
as.data.frame( countData[,1], drop=false)

#To extract the rownames of the countData file:
rownames(countData)
write.table(rownames(countData,"countData_annot"))

#Do it in Bash
cut -f1 countData_WT > counData_annot

```
=============
Run programme
=============
```R
source("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/JTK_CYCLEv3.1.R")
project <- "Frq08KO"
options(stringsAsFactors=FALSE)
annot <- read.delim("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/countData_annot")
data <- read.delim("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/countData_FrqKO")

#rownames(data) <- data[,1]
#data <- data[,-1]


jtkdist(4, 3)      # 4 total time points, 2 replicates per time point

periods <- 4      # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,6)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/JTK_WT53_2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
