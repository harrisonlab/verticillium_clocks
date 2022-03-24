Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep1and2/TPM_id2.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep1and2/JTK_WT53_TPMrep12.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep1and3/TPM_id2.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep1and3/JTK_WT53_TPMrep13.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep2and3/TPM_id2.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep2and3/JTK_WT53_TPMrep23.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```


## Normalized

Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep1and2/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep1and2/JTK_WT53_TPMquantilerep12.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep1and3/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep1and3/JTK_WT53_TPMquantilerep13.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG/Rep2and3/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

periods <- 2:12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/Rep2and3/JTK_WT53_TPMquantilerep23.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
