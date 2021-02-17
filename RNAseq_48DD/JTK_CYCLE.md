JTK_CYCLE 

```r
# Load libraries

source("JTKversion3/JTK_CYCLEv3.1.R")
library(docker4seq)
```

```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/fpkmnorm4jtk.txt")

jtkdist(12, 3)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/fpkm_norm_quantile4jtk.txt")

jtkdist(12, 3)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#Rep1
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep1/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep1/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep1/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep1/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#Rep2
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep2/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep2/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep2/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep2/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#Rep3
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep3/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep3/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/Rep3/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/Rep3/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
```r
#RepAB
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepAB/fpkmnorm4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepAB/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepAB/fpkm_norm_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepAB/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#RepAC
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepAC/fpkmnorm4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepAC/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepAC/fpkm_norm_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepAC/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
```r
#RepBC
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepBC/fpkmnorm4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepBC/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("fpkm_norm/RepBC/fpkm_norm_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/fpkm_norm/RepBC/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

-----------------------
```r
#JTK_TPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/TPM_id2.txt")

jtkdist(12, 3)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTK_TPM_quantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/TPM_quantile4jtk.txt")

jtkdist(12, 3)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/JTK_WT53_tpmquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#Rep1
#JTK_TPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep1/TPM4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep1/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTK_TPM_quantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep1/TPM_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep1/JTK_WT53_tpmquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#Rep2
#JTK_TPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep2/TPM4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep2/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTK_TPM_quantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep2/TPM_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep2/JTK_WT53_tpmquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
```r
#Rep3
#JTK_TPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep3/TPM4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep3/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTK_TPM_quantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/Rep3/TPM_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/Rep3/JTK_WT53_tpmquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#RepAB
#JTKTPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepAB/TPM4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepAB/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepAB/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepAB/JTK_WT53_TPMquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#RepAC
#JTKTPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepAC/TPM4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepAC/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepAC/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepAC/JTK_WT53_TPMquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#RepBC
#JTKTPM

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepBC/TPM4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepBC/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

#JTKfpkmnormquantiled

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("TPM/RepBC/TPM_quantile4jtk.txt")

jtkdist(12, 2)       # 12 total time points, 3 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/CleanRNA/TPM/RepBC/JTK_WT53_TPMquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```