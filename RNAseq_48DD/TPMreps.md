Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep1/TPM4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep1/JTK_WT53_TPMrep1.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep2/TPM4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep2/JTK_WT53_TPMrep2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep3/TPM4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep3/JTK_WT53_TPMrep3.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```
Q1<-read.delim("CorrectTPM/Rep1/JTK_WT53_TPMrep1Pval.txt",header=T,sep="\t")
82 genes
Q2<-read.table("CorrectTPM/Rep2/JTK_WT53_TPMrep2Pval.txt",header=T,sep="\t")
375 genes
Q3<-read.table("CorrectTPM/Rep3/JTK_WT53_TPMrep3Pval.txt",header=T,sep="\t")
308 genes


Q4<-merge(Q1,Q2, by.x="ID",by.y="ID")
4 genes
Q5<-merge(Q1,Q3, by.x="ID",by.y="ID")
1 genes
Q6<-merge(Q2,Q3, by.x="ID",by.y="ID")
10 genes
```



## Normalized


Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep1/TPM_quantile4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep1/JTK_WT53_quantileTPMrep1.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep2/TPM_quantile4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep2/JTK_WT53_quantileTPMrep2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/Rep3/TPM_quantile4jtk.txt")

jtkdist(12, 1)       # 12 total time points, 3 replicates per time point

periods <- 5:7       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(count,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,count)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/Rep3/JTK_WT53_quantileTPMrep3.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```


```
T1<-read.delim("CorrectTPM/Rep1/JTK_WT53_quantileTPMrep1Pval.txt",header=T,sep="\t")
130 genes
T2<-read.table("CorrectTPM/Rep2/JTK_WT53_quantileTPMrep2Pval.txt",header=T,sep="\t")
412 genes
T3<-read.table("CorrectTPM/Rep3/JTK_WT53_quantileTPMrep3Pval.txt",header=T,sep="\t")
201 genes


T4<-merge(T1,T2, by.x="ID",by.y="ID")
5 genes
T5<-merge(T1,T3, by.x="ID",by.y="ID")
6 genes
T6<-merge(T2,T3, by.x="ID",by.y="ID")
8 genes

rawvsquantileTPM
P1<-merge(T1,Q1, by.x="ID",by.y="ID")
14
P2<-merge(T2,Q2, by.x="ID",by.y="ID")
186
P3<-merge(T3,Q3, by.x="ID",by.y="ID")
91
```


----------------------------------------------------------