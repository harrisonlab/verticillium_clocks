```
ECHO has replaced JTK in most of my analysis. Thefore the commands here are deprecated and not used for the final analysis. This file contains multiple test runs in JTK_Cycle using different datasets. 
```
# JTK_CYCLE 

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

Explore results

A1<-read.delim("fpkm_norm/JTK_WT53_fpkmnorm_pval.txt",header=T)
129 genes
A2<-read.table("fpkm_norm/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
130
A3<-merge(A1,A2, by.x="ID",by.y="ID")
87

B1<-read.delim("TPM/JTK_WT53_TPM_pval.txt",header=T)
108
B2<-read.delim("TPM/JTK_WT53_tpmquantiled_pval.txt",header=T)
157
B3<-merge(B1,B2, by.x="ID",by.y="ID")
58

C1<-merge(A1,B1, by.x="ID",by.y="ID")
23
C2<-merge(A2,B2, by.x="ID",by.y="ID")
59
C3<-merge(A3,B3, by.x="ID",by.y="ID")
20

D1<-read.delim("fpkm_norm/Rep1/JTK_WT53_fpkmnorm_pval.txt",header=T)
143 genes
D2<-read.table("fpkm_norm/Rep1/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
142
D3<-merge(D1,D2, by.x="ID",by.y="ID")
72

E1<-read.delim("fpkm_norm/Rep2/JTK_WT53_fpkmnorm_pval.txt",header=T)
143 genes
E2<-read.table("fpkm_norm/Rep2/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
142
E3<-merge(E1,E2, by.x="ID",by.y="ID")
292

F1<-read.delim("fpkm_norm/Rep3/JTK_WT53_fpkmnorm_pval.txt",header=T)
174 genes
F2<-read.table("fpkm_norm/Rep3/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
197
F3<-merge(F1,F2, by.x="ID",by.y="ID")
98

Z1<-merge(D1,E1, by.x="ID",by.y="ID")
6
Z2<-merge(D1,F1, by.x="ID",by.y="ID")
3
Z3<-merge(E1,F1, by.x="ID",by.y="ID")
6
Z4<-merge(D2,E2, by.x="ID",by.y="ID")
10
Z5<-merge(D2,F2, by.x="ID",by.y="ID")
3
Z6<-merge(E2,F2, by.x="ID",by.y="ID")
8

G1<-read.delim("fpkm_norm/RepAB/JTK_WT53_fpkmnorm_pval.txt",header=T)
157 genes
G2<-read.table("fpkm_norm/RepAB/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
176
G3<-merge(G1,G2, by.x="ID",by.y="ID")
104

H1<-read.delim("fpkm_norm/RepAC/JTK_WT53_fpkmnorm_pval.txt",header=T)
95 genes
H2<-read.table("fpkm_norm/RepAC/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
97
H3<-merge(H1,H2, by.x="ID",by.y="ID")
54

I1<-read.delim("fpkm_norm/RepBC/JTK_WT53_fpkmnorm_pval.txt",header=T)
176 genes
I2<-read.table("fpkm_norm/RepBC/JTK_WT53_fpkmnormquantiled_pval.txt",header=T,sep="\t")
178
I3<-merge(I1,I2, by.x="ID",by.y="ID")
107

Z1<-merge(G1,H1, by.x="ID",by.y="ID")
7
Z2<-merge(G1,I1, by.x="ID",by.y="ID")
3
Z3<-merge(H1,I1, by.x="ID",by.y="ID")
11
X4<-merge(G2,H2, by.x="ID",by.y="ID")
9
X5<-merge(G2,I2, by.x="ID",by.y="ID")
5
X6<-merge(H2,I2, by.x="ID",by.y="ID")
8

fpkm1<-read.delim("fpkm_norm/JTK_WT53_fpkmnorm.txt",header=T)
Inter<-read.delim("interproscan_gene.txt",header=T)
Inter2<-merge(fpkm1,Inter, by.x="ID",by.y="ID")
write.table(Inter2,"fpkm_robust_interpro.txt",sep="\t",na="",quote=F)

fpkm2<-read.delim("fpkm_norm/JTK_WT53_fpkmnormquantiled.txt",header=T)
Inter<-read.delim("interproscan_gene.txt",header=T)
Inter3<-merge(fpkm2,Inter, by.x="ID",by.y="ID")
write.table(Inter3,"fpkm_robust_quantiled_interpro.txt",sep="\t",na="",quote=F)

# JTK using single reps

fpkm

## Raw

Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep1/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep1/JTK_WT53_fpkmrep1.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep2/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep2/JTK_WT53_fpkmrep2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep3/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep3/JTK_WT53_fpkmrep3.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```


```
N1<-read.delim("Correctfpkm/Rep1/JTK_WT53_fpkmrep1Pval.txt",header=T,sep="\t")
143 genes
N2<-read.table("Correctfpkm/Rep2/JTK_WT53_fpkmrep2Pval.txt",header=T,sep="\t")
405 genes
N3<-read.table("Correctfpkm/Rep3/JTK_WT53_fpkmrep3Pval.txt",header=T,sep="\t")
174 genes


N4<-merge(N1,N2, by.x="ID",by.y="ID")
#6 genes
N5<-merge(N1,N3, by.x="ID",by.y="ID")
3 genes
N6<-merge(N2,N3, by.x="ID",by.y="ID")
6 genes
```


## Normalized

Rep1
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep1/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep1/JTK_WT53_fpkmquantilerep1.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep2
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep2/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep2/JTK_WT53_fpkmquantilerep2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

Rep3
```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/Rep3/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/Rep3/JTK_WT53_fpkmquantilerep3.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```
M1<-read.delim("Correctfpkm/Rep1/JTK_WT53_fpkmquantilerep1Pval.txt",header=T,sep="\t")
142 genes
M2<-read.table("Correctfpkm/Rep2/JTK_WT53_fpkmquantilerep2Pval.txt",header=T,sep="\t")
438 genes
M3<-read.table("Correctfpkm/Rep3/JTK_WT53_fpkmquantilerep3Pval.txt",header=T,sep="\t")
197 genes


M4<-merge(M1,M2, by.x="ID",by.y="ID")
10 genes
M5<-merge(M1,M3, by.x="ID",by.y="ID")
3 genes
M6<-merge(M2,M3, by.x="ID",by.y="ID")
8 genes
```


rawvsquantileTPM
O1<-merge(N1,M1, by.x="ID",by.y="ID")
72
O2<-merge(N2,M2, by.x="ID",by.y="ID")
292
O3<-merge(N3,M3, by.x="ID",by.y="ID")
98

# New period

```r
#JTKfpkmnorm


project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/fpkmnorm4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/JTK_WT53_fpkmnorm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("Correctfpkm/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/Correctfpkm/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```


```r
#TPM


project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/TPM_id2.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("CorrectTPM/TPM_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/CorrectTPM/JTK_WT53_TPMquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```