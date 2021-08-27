

## Quantile normalization

```r
# Load libraries

source("JTKversion3/JTK_CYCLEv3.1.R")
library(docker4seq)
```

```r
#fpkm norm

df <- read.table("fpkm_norm/quantile/toquantile.txt",header=T)

head(data)

rownames(df)<-df$ID

df_rank <- apply(df,2,rank,ties.method="min")
df_sorted <- data.frame(apply(df, 2, sort))


quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

new_data <- quantile_normalisation(df[,2:37])
new_data

write.table(new_data,"fpkm_norm/quantile/fpkm_norm_quantile.txt",sep="\t",na="",quote=F)
```

```r
#JTKfpkmnorm

project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")
count <- read.delim("vAG2/fpkm_norm_quantile4jtk.txt")

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG2/JTK_WT53_fpkmnormquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```


```r
#TPM

df <- read.table("TPM/TPM_id.txt",header=T)

head(df)

rownames(df)<-df$ID

df_rank <- apply(df,2,rank,ties.method="min")
df_sorted <- data.frame(apply(df, 2, sort))


quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

new_data <- quantile_normalisation(df[,2:37])
new_data

write.table(new_data,"TPM/TPM_quantile.txt",sep="\t",na="",quote=F)
```

```r
#JTKTPMquantile

count <- read.delim("vAG/TPM_quantile4jtk.txt")

jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/JTK_WT53_TPMquantiled.txt"),row.names=F,col.names=T,quote=F,sep="\t")




Raw
```r
#JTKTPMquantile

count <- read.delim("vAG2/fpkmnorm4jtk.txt")

jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG2/JTK_WT53_fpkm.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```

```r
#JTKTPMquantile

count <- read.delim("vAG/TPM_id2.txt")

jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/vAG/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```




C1<-read.delim("vAG/JTK_WT53_TPMPval.txt",header=T)
56 genes
C2<-read.table("vAG2/JTK_WT53_fpkmPval.txt",header=T,sep="\t")
79
C3<-merge(C1,C2, by.x="ID",by.y="ID")
19

#fpkmvsTPM quantiled

A1<-read.delim("vAG2/JTK_WT53_fpkmnormquantiledPval.txt",header=T)
74
A2<-read.delim("vAG/JTK_WT53_TPMquantiledPval.txt",header=T)
89
A3<-merge(A1,A2, by.x="ID",by.y="ID")
43
write.table(A3,"quantvs.txt",sep="\t",na="",quote=F)

C4<-merge(A1,C2, by.x="ID",by.y="ID")
31
C5<-merge(A2,C1, by.x="ID",by.y="ID")
27

C6<-merge(C4,C5, by.x="ID",by.y="ID")
C7<-merge(C3,A3, by.x="ID",by.y="ID")

C8<-merge(C6,C7, by.x="ID",by.y="ID")
write.table(C8,"common.txt",sep="\t",na="",quote=F)



Q4<-merge(Q1,T1, by.x="ID",by.y="ID")
#109 genes
Q5<-merge(Q2,T2, by.x="ID",by.y="ID")
#35 genes






2 reps

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")

count <- read.delim("fpkm_norm/Rep2and3/fpkmnorm4jtk.txt")

jtkdist(12, 2)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/fpkm_norm/Rep2and3/JTK_WT53_norm_fpkm.txt"),row.names=F,col.names=T,quote=F,sep="\t")

2 reps

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")

count <- read.delim("TPM/Rep2and3/TPM4jtk.txt")

jtkdist(12, 2)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/TPM/Rep2and3/JTK_WT53_TPM.txt"),row.names=F,col.names=T,quote=F,sep="\t")

1 reps

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")

count <- read.delim("TPM/Rep1/TPM4jtk.txt")

jtkdist(12, 1)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/TPM/Rep1/JTK_WT53_TPM_1rep.txt"),row.names=F,col.names=T,quote=F,sep="\t")
1 reps

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")

count <- read.delim("TPM/Rep1/TPM4jtk.txt")

jtkdist(12, 1)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/TPM/Rep1/JTK_WT53_TPM_1rep.txt"),row.names=F,col.names=T,quote=F,sep="\t")
'''

#fpkmvsTPM quantiled

A1<-read.delim("vAG2/JTK_WT53_fpkmnormquantiledPval.txt",header=T)
74
A2<-read.delim("vAG/JTK_WT53_TPMquantiledPval.txt",header=T)
89
A3<-merge(A1,A2, by.x="ID",by.y="ID")
43



cat interproscan.txt | cut -f1 | awk -F'.t1' '{print $1}' > table.tsv
DD<-read.delim("interproscan2.txt",header=T)
W2<-read.delim("vAG2/JTK_WT53_fpkmnormquantiled.txt",header=T)
W1<-merge(W2,DD, by.x="ID",by.y="ID")
write.table(W1,"withinterproscan.txt",sep="\t",na="",quote=F)

P2<-merge(Q1,DD, by.x="ID",by.y="ID")
write.table(P2,"JTKCycle_ResultsQuantileTPM.txt",sep="\t",na="",quote=F)



A1<-read.delim("vAG2/JTK_WT53_fpkmnormquantiledPval.txt",header=T)
74
A2<-read.delim("vAG/JTK_WT53_TPMquantiledPval.txt",header=T)
89
A3<-merge(A1,A2, by.x="ID",by.y="ID")
43


write.table(T3,"fpkmquantilevsfpkm_vZZ.txt",sep="\t",na="",quote=F)


T4<-read.delim("fpkm_norm/Rep2and3/JTK_Results.txt",header=T)
#442 genes
T5<-read.table("fpkm_norm/Rep1/JTK_Results.txt",header=T,sep="\t")
#77 genes


T6<-merge(T4,T5, by.x="ID",by.y="ID")
#2 genes

T7<-merge(T4,T1, by.x="ID",by.y="ID")
#65 genes

T8<-merge(T4,T2, by.x="ID",by.y="ID")
#63 genes

T9<-merge(T7,T8, by.x="ID",by.y="ID")
#47 genes




#TPM

Q1<-read.delim("TPM/quantiled_vZZ_Results.txt",header=T)
#110 genes ... 89 genes
Q2<-read.table("TPM/vZZ_Results.txt",header=T,sep="\t")
#78 genes ... 56 genes
Q3<-merge(Q1,Q2, by.x="ID",by.y="ID")
#41 genes ... 27 genes
write.table(Q3,"TPMquantilevsTPM_vZZ.txt",sep="\t",na="",quote=F)

Q4<-merge(Q1,T1, by.x="ID",by.y="ID")
#109 genes
Q5<-merge(Q2,T2, by.x="ID",by.y="ID")
#35 genes




Q6<-read.delim("TPM/Rep2and3/JTK_Results.txt",header=T)
#734 genes
Q7<-read.table("TPM/Rep1/JTK_Results.txt",header=T,sep="\t")
#53 genes


Q8<-merge(Q6,Q7, by.x="ID",by.y="ID")
#1 genes

Q9<-merge(Q6,Q1, by.x="ID",by.y="ID")
#56 genes

Q10<-merge(Q6,Q2, by.x="ID",by.y="ID")
#48 genes

Q11<-merge(T4,Q6, by.x="ID",by.y="ID")
245

Q12<-merge(T1,Q1, by.x="ID",by.y="ID")

T1 or Q1

cat interproscan.txt | cut -f1 | awk -F'.t1' '{print $1}' > table.tsv
DD<-read.delim("interproscan2.txt",header=T)
P1<-merge(T1,DD, by.x="ID",by.y="ID")
write.table(P1,"JTKCycle_ResultsQuantile.txt",sep="\t",na="",quote=F)

P2<-merge(Q1,DD, by.x="ID",by.y="ID")
write.table(P2,"JTKCycle_ResultsQuantileTPM.txt",sep="\t",na="",quote=F)



source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("countData_annot.txt")

count <- read.delim("vAG2/fpkmnorm4jtk.txt")

jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

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
write.table(results,file=paste("/Users/antoniogomez/Desktop/RNAseqVERT/TPM/Rep1/JTK_WT53_TPM_1rep.txt"),row.names=F,col.names=T,quote=F,sep="\t")