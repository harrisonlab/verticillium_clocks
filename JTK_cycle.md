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
annot <- read.delim("/Users/lopeze/Desktop/Bioinformatics/RNA_Seq/JTK_Cycle/countData_annot")
data <- read.delim("/Users/lopeze/Desktop/Bioinformatics/RNA_Seq/JTK_Cycle/countData_WT_duplicate.txt")

#rownames(data) <- data[,1]
#data <- data[,-1]
jtkdist(8, 3)       # 13 total time points, 2 replicates per time point

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
write.table(results,file=paste("/Users/lopeze/Desktop/Bioinformatics/JTK_Cycle/JTK_WT53_duplicate.txt"),row.names=F,col.names=T,quote=F,sep="\t")
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
write.table(results,file=paste("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/JTK_cycle53"),row.names=F,col.names=T,quote=F,sep="\t")
```

Heatmap_JTK

```R

count=read.csv("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/JTK.csv",header=T)
RM_dd6<-log(rowMeans(count[,7:9] +0.000001))
RM_dd12<-log(rowMeans(count[,10:12]+0.000001))
RM_dd18<-log(rowMeans(count[,13:15]+0.000001))
RM_dd24<-log(rowMeans(count[,16:18]+0.000001))
count2=cbind(count[,1:6],RM_dd6,RM_dd12,RM_dd18,RM_dd24)
count2[,7]=count2[,7]-median(count2[,7])
count2[,8]=count2[,8]-median(count2[,8])
count2[,9]=count2[,9]-median(count2[,9])
count2[,10]=count2[,10]-median(count2[,10])
count3=count2[which(count2[,3]<0.01),]
#count4=count3[which(count3[,5]==6),]
#count4=count3[which(count3[,5]==3),]
count4=count3

library("gplots")
orthologs_matrix<-data.matrix(count4[,c(7:10)])
breaks = seq(min(orthologs_matrix),max(orthologs_matrix),length.out=1000)
gradient1 = colorpanel( sum( breaks[-1]<=1 ),"blue", "black" )
gradient2 = colorpanel( sum( breaks[-1]>1 ), "black", "yellow")
hm.colors = c(gradient1,gradient2)

#cluster by pattern, not phase offset
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/1.pdf",width=7,height=7)
heatmap.2(as.matrix(orthologs_matrix),scale="none",breaks=breaks,col=hm.colors,
          Colv=FALSE,dendrogram="row",trace="none",
          margin=c(5,10),lwid=c(1.5,2.0))

#Order by phase offset
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/2.pdf",width=7,height=7)
count4=count3[order (count3[,5]),]
heatmap.2(as.matrix(orthologs_matrix),scale="none",breaks=breaks,col=hm.colors,
 Colv=FALSE,Rowv=FALSE,dendrogram="none",trace="none")




######PLOT THE LOG WITH CENTERING
count=read.csv("JTK.csv",header=T)
RM_dd6<-log(rowMeans(count[,7:9] +0.000001))
RM_dd12<-log(rowMeans(count[,10:12]+0.000001))
RM_dd18<-log(rowMeans(count[,13:15]+0.000001))
RM_dd24<-log(rowMeans(count[,16:18]+0.000001))
count2=cbind(count[,1:6],RM_dd6,RM_dd12,RM_dd18,RM_dd24)
#count2[,7]=count2[,7]
#count2[,8]=count2[,8]
#count2[,9]=count2[,9]
#count2[,10]=count2[,10]
count2[,7]=count2[,7]-median(count2[,7])
count2[,8]=count2[,8]-median(count2[,8])
count2[,9]=count2[,9]-median(count2[,9])
count2[,10]=count2[,10]-median(count2[,10])
count3=count2[which(count2[,3]<0.05),]
#count4=count3[which(count3[,5]==6),]
#count4=count3[which(count3[,5]==3),]
count4=count3
library(lattice)
count5=melt(count4[,c(1,7:10)], id.vars = c("gene"))
pdf(file = "cycle.pdf",width=70,height=70)
xyplot(value ~ variable| gene, data = count5,
type = "l",
lty = c(1, 2, 2, 1),
lwd = c(1, 1, 1, 3),
col.line = c(rep("black",3), "red"))
dev.off()



##PLOT THE RPKM NO CENTERING
count=read.csv("JTK.csv",header=T)
RM_dd6<-(rowMeans(count[,7:9] ))
RM_dd12<-(rowMeans(count[,10:12]))
RM_dd18<-(rowMeans(count[,13:15]))
RM_dd24<-(rowMeans(count[,16:18]))
count2=cbind(count[,1:6],RM_dd6,RM_dd12,RM_dd18,RM_dd24)
count2[,7]=count2[,7]
count2[,8]=count2[,8]
count2[,9]=count2[,9]
count2[,10]=count2[,10]
count3=count2[which(count2[,3]<0.05),]
#count4=count3[which(count3[,5]==6),]
#count4=count3[which(count3[,5]==3),]
count4=count3
library(lattice)
count5=melt(count4[,c(1,7:10)], id.vars = c("gene"))
pdf(file = "cycle2.pdf",width=70,height=70)
xyplot(value ~ variable| gene, data = count5,
type = "l",
lty = c(1, 2, 2, 1),
lwd = c(1, 1, 1, 3),
col.line = c(rep("black",3), "red"))
dev.off()


######PLOT THE LOG WITHOUT CENTERING
count=read.csv("JTK.csv",header=T)
RM_dd6<-log(rowMeans(count[,7:9] +0.000001))
RM_dd12<-log(rowMeans(count[,10:12]+0.000001))
RM_dd18<-log(rowMeans(count[,13:15]+0.000001))
RM_dd24<-log(rowMeans(count[,16:18]+0.000001))
count2=cbind(count[,1:6],RM_dd6,RM_dd12,RM_dd18,RM_dd24)
count2[,7]=count2[,7]
count2[,8]=count2[,8]
count2[,9]=count2[,9]
count2[,10]=count2[,10]
#count2[,7]=count2[,7]-median(count2[,7])
#count2[,8]=count2[,8]-median(count2[,8])
#count2[,9]=count2[,9]-median(count2[,9])
#count2[,10]=count2[,10]-median(count2[,10])
count3=count2[which(count2[,3]<0.05),]
#count4=count3[which(count3[,5]==6),]
#count4=count3[which(count3[,5]==3),]
count4=count3
library(lattice)
count5=melt(count4[,c(1,7:10)], id.vars = c("gene"))
pdf(file = "cycle3.pdf",width=70,height=70)
xyplot(value ~ variable| gene, data = count5,
type = "l",
lty = c(1, 2, 2, 1),
lwd = c(1, 1, 1, 3),
col.line = c(rep("black",3), "red"))
dev.off()


#TimecourseLD

source("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/JTK_CYCLEv3.1.R")

project <- "Timecourse_LD"

options(stringsAsFactors=FALSE)
annot <- read.delim("/Users/lopeze/Desktop/Statistics_R/Timecourse_LD/Timecourse_LD_annot.txt")
data <- read.delim("/Users/lopeze/Desktop/Statistics_R/Timecourse_LD/Timecourse_LD_data.txt")

rownames(data) <- data[,1]
data <- data[,-1]
jtkdist(6, 3)       # 13 total time points, 2 replicates per time point

periods <- 5:6       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

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

save(results,file=paste("/Users/lopeze/Desktop/Statistics_R/Timecourse_LD/JTK_TimecourseLD",project,"rda",sep="."))
write.table(results,file=paste("/Users/lopeze/Desktop/Statistics_R/Timecourse_LD/JTK_TimecourseLD.txt"),row.names=F,col.names=T,quote=F,sep="\t")


#Timecourse_T

source("/Users/lopeze/Desktop/Bioinformatics/RNA_seq/JTK_Cycle/JTK_CYCLEv3.1.R")

project <- "Timecourse_LD"

options(stringsAsFactors=FALSE)
annot <- read.delim("/Users/lopeze/Desktop/Statistics_R/Timecourse_T/Timecourse_T_annot.txt")
data <- read.delim("/Users/lopeze/Desktop/Statistics_R/Timecourse_T/Timecourse_T_data.txt")

rownames(data) <- data[,1]
data <- data[,-1]
jtkdist(6, 2)       # 13 total time points, 2 replicates per time point

periods <- 5:6       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

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

save(results,file=paste("/Users/lopeze/Desktop/Statistics_R/Timecourse_LD/JTK_TimecourseLD",project,"rda",sep="."))
write.table(results,file=paste("/Users/lopeze/Desktop/Statistics_R/Timecourse_T/JTK_TimecourseT.txt"),row.names=F,col.names=T,quote=F,sep="\t")
