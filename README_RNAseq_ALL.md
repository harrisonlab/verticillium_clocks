

#Experiment all
```R
#Gene clustering of different group of samples
#Eliminate Frq08_DD24_rep3 sample from colData and countData
#Eliminate Frq53_LL06_rep2 sample from colData and countData
#Eliminate WT08_DD18_rep3 sample from colData and countData

library(DESeq2)
colData <- read.table("colData_all",header=T,sep="\t")
countData <- read.table("countData_all",header=T,sep="\t")
colData$Group <- factor(paste0(colData$Strain,colData$Light,colData$Time))

colData <- colData[!(colData$Sample=="WT08_DD18_rep3"),]      
countData <- subset(countData, select=-WT08_DD18_rep3)

design <- ~Group
dds <- 	DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- DESeq(dds)
rld <- rlog( dds, blind=FALSE )
resultsNames(dds)

=================
PCA plots elipses
=================
library("ggplot2")
library("ggrepel")

data <- plotPCA(rld, intgroup=c("Strain", "Light"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Strain)) +
 geom_point(size=3) +
 theme_minimal()+
 theme(axis.title.x = element_text(size=30),
 axis.title.y = element_text(size=30),
 legend.text=element_text(size=20),
 legend.title=element_text(size=20),
 axis.text.x = element_text(colour="black",size=20),
 axis.text.y = element_text(colour="black",size=20),)+
 stat_ellipse(aes(x=PC1,y=PC2, fill=Strain),
              geom="polygon", level=0.95, alpha=0.2) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) +
 geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()
ggsave("PCA_sample_names_elipses4.pdf", pca_plot, dpi=300, height=15, width=20)

=================
PCA plots
=================
library("ggplot2")
library("ggrepel")

data <- plotPCA(rld, intgroup=c("Strain", "Light"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, shape=Light, color=Strain)) +
 geom_point(size=5) +
 theme_minimal()+
 theme(axis.title.x = element_text(size=30),
 axis.title.y = element_text(size=30),
 legend.text=element_text(size=25),
 legend.title=element_text(size=30),
 axis.text.x = element_text(colour="black",size=25),
 axis.text.y = element_text(colour="black",size=25),)+
 #stat_ellipse(aes(x=PC1,y=PC2, fill=Strain),
              #geom="polygon", level=0.95, alpha=0.2) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) +
 geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("PCA_sample_names_thesis.pdf", pca_plot, dpi=300, height=15, width=20)

=================
Sample distanes
=================

sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sample)
colnames(sampleDistMatrix) <- paste(rld$Sample)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins = c(11, 11))

dev.off()

===============
Analysis of gene expression
===============
https://support.bioconductor.org/p/84366/
https://support.bioconductor.org/p/62355/

#Relevel to WT53 D
dds$Light <- relevel(dds$Light, "d")
dds$Strain <- relevel(dds$Strain, "WT53")

## 53 L vs D (treatment effect for genotype WT53)
resultsNames(dds)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "WT53l06h", "WT53d06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < 0, ]

write.table(sig.res.upregulated,"WT53_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT53_LD_down.txt",sep="\t",quote=F)
write.table(res,"WT53_LD.txt",sep="\t",quote=F)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]

write.table(sig.res.upregulated,"WT53_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT53_LD_down_LFC.txt",sep="\t",quote=F)

summary(res)

##08 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "WT08l6h", "WT08d6h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"WT08_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT08_LD_down.txt",sep="\t",quote=F)
write.table(res,"WT08_LD.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"WT08_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT08_LD_down_LFC.txt",sep="\t",quote=F)


##Wc1 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc153l06h", "Wc153d06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Wc1_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1_LD_down.txt",sep="\t",quote=F)
write.table(res,"Wc1_LD.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Wc1_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1_LD_down_LFC.txt",sep="\t",quote=F)


##WC1 vs WT53 D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc153d06h", "WT53d06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Wc1vsWT53_D_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1vsWT53_D_down.txt",sep="\t",quote=F)
write.table(res,"Wc1vsWT53_D.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Wc1vsWT53_D_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1vsWT53_D_down_LFC.txt",sep="\t",quote=F)


##WC1 vs WT53 L
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc153l06h", "WT53l06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Wc1vsWT53_L_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1vsWT53_L_down.txt",sep="\t",quote=F)
write.table(res,"Wc1vsWT53_L.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Wc1vsWT53_L_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wc1vsWT53_L_down_LFC.txt",sep="\t",quote=F)

##Frq53 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq53l6h", "Frq53d6h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq53_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53_LD_down.txt",sep="\t",quote=F)
write.table(res,"Frq53_LD.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq53_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53_LD_down_LFC.txt",sep="\t",quote=F)

##Frq53 vs WT53 D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq53d6h", "WT53d06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq53vsWT53_D_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53vsWT53_D_down.txt",sep="\t",quote=F)
write.table(res,"Frq53vsWT53_D.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq53vsWT53_D_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53vsWT53_D_down_LFC.txt",sep="\t",quote=F)

##Frq53 vs WT53 L
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq53l6h", "WT53l06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq53vsWT53_L_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53vsWT53_L_down.txt",sep="\t",quote=F)
write.table(res,"Frq53vsWT53_L.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq53vsWT53_L_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq53vsWT53_L_down_LFC.txt",sep="\t",quote=F)

##Frq08 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq08l06h", "Frq08d06h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq08_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08_LD_down.txt",sep="\t",quote=F)
write.table(res,"Frq09_LD.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq08_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08_LD_down_LFC.txt",sep="\t",quote=F)


##Frq08 vs WT53 D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq08d06h", "WT08d6h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq08vsWT08_D_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08vsWT08_D_down.txt",sep="\t",quote=F)
write.table(res,"Frq08vsWT08_D.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq08vsWT08_D_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08vsWT08_D_down_LFC.txt",sep="\t",quote=F)

##Frq08 vs WT08 L
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Frq08l06h", "WT08l6h"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"Frq08vsWT08_L_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08vsWT08_L_down.txt",sep="\t",quote=F)
write.table(res,"Frq08vsWT08_L.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"Frq08vsWT08_L_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08vsWT08_L_down_LFC.txt",sep="\t",quote=F)
```

===============
Dispersion
==============
```R
plotDispEsts(dds)
dev.off()
```

===============
MA plot and histogram and independent filtering
==============

```R
plotMA(res,ylim=c(-5,5))
dev.off()
```
#histogram

```R
hist(res$pvalue, breaks=20, col="grey50", border="white")
dev.off()
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()
```

#Independent filtering
```R
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:6/6))
bins <- cut(res$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(res$pvalue, bins, function(p)
                          mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
                     ylab = "fraction of small p values")

```

================
Preparing a file
================
#Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"normalised_counts.txt",sep="\t",na="",quote=F)

## When the raw_data is created, the column names are shifted, due to the lack of column name for the genes. To fx it, do nano filex, and then added the word "Gene", press tab and save ctr+x

install.packages("XVector", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

library(XVector)
library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)



# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
#Maria produced the FPKM files
input=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
#CDS lengths to use for FPKM calculations
python /home/sobczm/bin/popgen/renseq/write_seq_length.py Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa
#Remove the .t1 etc. transcript suffix in the file with CDS lengths info
#There is only one transcript per gene, so no need for further manipulation
sed -i 's/.t[0-9]//' Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt
#Add headers to the length table and fix the header in the count table.
sed -i '1s/^/Gene\tLength\n/' Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt
#Calculate FPKM with a new script
scripts=/home/sobczm/bin/popgen/rnaseq

#countData file may need to be changed and add a Gene column at the beginning
python $scripts/calculate_fpkm.py countData_all Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt

```

#Produce a detailed table of analyses
'''python
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

Run with commands:

#For WC1
```bash
for GeneGff in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3); do
    Strain=JR2
    Organism=V.dahliae
    Assembly=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
    TFs=$(ls analysis/transcription_factors/V.dahliae/JR2/JR2_TF_domains.tsv)
    InterPro=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
    Antismash=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/secondary_metabolites/antismash/JR2/fungi-38c21e2b-ce4f-4026-8f65-536412b28aee/geneclusters.txt)
    SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa)
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_all/ALL)
    DEG_Files=$(ls \
        $Dir1/WT53_LD.txt \
        $Dir1/Wc1_LD.txt \
        $Dir1/Wc1vsWT53_D.txt \
        $Dir1/Wc1vsWT53_L.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts.txt)
    FPKM=$(ls $Dir1/countData_all.fpkm)
    #ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_WC1_.tsv
done
```
#For Frq53 file

```bash
for GeneGff in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3); do
    Strain=JR2
    Organism=V.dahliae
    Assembly=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
    TFs=$(ls analysis/transcription_factors/V.dahliae/JR2/JR2_TF_domains.tsv)
    InterPro=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
    Antismash=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/secondary_metabolites/antismash/JR2/fungi-38c21e2b-ce4f-4026-8f65-536412b28aee/geneclusters.txt)
    SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa)
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_all/ALL)
    DEG_Files=$(ls \
        $Dir1/WT53_LD.txt \
        $Dir1/Frq53_LD.txt \
        $Dir1/Frq53vsWT53_D.txt \
        $Dir1/Frq53vsWT53_L.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts.txt)
    FPKM=$(ls $Dir1/countData_all.fpkm)
    #ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_FRQ53.tsv
done
```
#For FRQ08

```bash
for GeneGff in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3); do
    Strain=JR2
    Organism=V.dahliae
    Assembly=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
    TFs=$(ls analysis/transcription_factors/V.dahliae/JR2/JR2_TF_domains.tsv)
    InterPro=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
    Antismash=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/secondary_metabolites/antismash/JR2/fungi-38c21e2b-ce4f-4026-8f65-536412b28aee/geneclusters.txt)
    SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa)
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_all/ALL)
    DEG_Files=$(ls \
    $Dir1/WT08_LD.txt \
    $Dir1/Frq08_LD.txt \
    $Dir1/Frq08vsWT08_D.txt \
    $Dir1/Frq08vsWT08_L.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts.txt)
    FPKM=$(ls $Dir1/countData_all.fpkm)
    #ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_FRQ08.tsv
done
```



========================================================================
Merging data tables for TFs
========================================================================
```bash
DataSet.x<-read.csv("gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_WC1_.csv", header=T)
DataSet.y<-read.csv("analysis/transcription_factors/V.dahliae/JR2/JR2_TF_domains_v2.csv", header = T)

total<-merge(DataSet.x,DataSet.y, by.x="transcript_id", by.y="transcript_id",sort=T,all.x = T)
write.csv(total, "JR2_gene_table_WC1_TFs.csv")
```

=========================================
Creating a subset table from a gene list
=========================================
```bash
for num in 1
do
RNASeqData=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_WC1_TFs.csv
GeneList=analysis/transcription_factors/V.dahliae/JR2/JR2_TF.txt
echo "TFs:"
New=${RNASeqData/.csv/_extract.csv}
cat $RNASeqData | head -n 1 > $New
cat $RNASeqData | grep -w -f $GeneList >> $New
cat $New | tail -n +2 | wc -l
done
```

========================================================================
Merging data tables for pigments
========================================================================
```bash
DataSet.x<-read.csv("gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_WC1_.csv", header=T)
DataSet.y<-read.csv("analysis/secondary_metabolism/pigments/pigments.csv", header = T)

total<-merge(DataSet.x,DataSet.y, by.x="transcript_id", by.y="transcript_id",sort=T,all.x = T)
write.csv(total, "JR2_gene_table_WC1_pigments.csv")
```

=========================================
Creating a subset table from a gene list
=========================================
```bash
for num in 1
do
RNASeqData=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_WC1_TFs.csv
GeneList=analysis/transcription_factors/V.dahliae/JR2/JR2_TF.txt
echo "TFs:"
New=${RNASeqData/.csv/_extract.csv}
cat $RNASeqData | head -n 1 > $New
cat $RNASeqData | grep -w -f $GeneList >> $New
cat $New | tail -n +2 | wc -l
done
```

#FUNCTIONAL annotations
##Analysis of DEGs vs all genes
```bash
OutDir=analysis/enrichment/experiment_WT/
InterProTSV=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv
ProgDir=/home/adamst/git_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
```

#Functional annotation of WT bl vs d

```bash
#Extract the first column of one file and save it as a *_name.file
cut -f1 file > file2

#Remove .pi from column and save to other file
sed -i.bak 's/.p1//' experiment_all_gene_GO_annots.tsv
```

#Wc1_LD up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_ALL
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1_LD/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1_LD_up_LFC_names.txt
AllGenes=$OutDir/Wc1_LD_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1_LD_up_DEGs.txt
Set2Genes=$OutDir/Wc1_LD_up2.txt
AllGenes=$OutDir/Wc1_LD_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1_LD/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1_LD_down_LFC_names.txt
AllGenes=$OutDir/Wc1_LD_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1_LD_down_DEGs.txt
Set2Genes=$OutDir/Wc1_LD_down2.txt
AllGenes=$OutDir/Wc1_LD_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

#Wc1/WT53 in D up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1vsWT53_D/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1vsWT53_D_up_LFC_names.txt
AllGenes=$OutDir/Wc1vsWT53_D_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1vsWT53_D_up_DEGs.txt
Set2Genes=$OutDir/Wc1vsWT53_D_up2.txt
AllGenes=$OutDir/Wc1vsWT53_D_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1vsWT53_D/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1vsWT53_D_down_LFC_names.txt
AllGenes=$OutDir/Wc1vsWT53_D_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1vsWT53_D_down_DEGs.txt
Set2Genes=$OutDir/Wc1vsWT53_D_down2.txt
AllGenes=$OutDir/Wc1vsWT53_D_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Wc1/WT53 in L up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1vsWT53_L/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1vsWT53_L_up_LFC_names.txt
AllGenes=$OutDir/Wc1vsWT53_L_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1vsWT53_L_up_DEGs.txt
Set2Genes=$OutDir/Wc1vsWT53_L_up2.txt
AllGenes=$OutDir/Wc1vsWT53_L_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Wc1/Wc1vsWT53_L/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Wc1vsWT53_L_down_LFC_names.txt
AllGenes=$OutDir/Wc1vsWT53_L_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Wc1vsWT53_L_down_DEGs.txt
Set2Genes=$OutDir/Wc1vsWT53_L_down2.txt
AllGenes=$OutDir/Wc1vsWT53_L_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

#Frq53_LD up and down regulated

```bash
WorkDir=analysis/enrichment/experiment_ALL
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53_LD/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53_LD_up_LFC_names.txt
AllGenes=$OutDir/Frq53_LD_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53_LD_up_DEGs.txt
Set2Genes=$OutDir/Frq53_LD_up2.txt
AllGenes=$OutDir/Frq53_LD_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53_LD/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53_LD_down_LFC_names.txt
AllGenes=$OutDir/Frq53_LD_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53_LD_down_DEGs.txt
Set2Genes=$OutDir/Frq53_LD_down2.txt
AllGenes=$OutDir/Frq53_LD_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

#Frq53/WT53 in D up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53vsWT53_D/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53vsWT53_D_up_LFC_names.txt
AllGenes=$OutDir/Frq53vsWT53_D_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53vsWT53_D_up_DEGs.txt
Set2Genes=$OutDir/Frq53vsWT53_D_up2.txt
AllGenes=$OutDir/Frq53vsWT53_D_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53vsWT53_D/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53vsWT53_D_down_LFC_names.txt
AllGenes=$OutDir/Frq53vsWT53_D_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53vsWT53_D_down_DEGs.txt
Set2Genes=$OutDir/Frq53vsWT53_D_down2.txt
AllGenes=$OutDir/Frq53vsWT53_D_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Frq53/WT53 in L up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53vsWT53_L/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53vsWT53_L_up_LFC_names.txt
AllGenes=$OutDir/Frq53vsWT53_L_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53vsWT53_L_up_DEGs.txt
Set2Genes=$OutDir/Frq53vsWT53_L_up2.txt
AllGenes=$OutDir/Frq53vsWT53_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq53/Frq53vsWT53_L/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq53vsWT53_L_down_LFC_names.txt
AllGenes=$OutDir/Frq53vsWT53_L_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq53vsWT53_L_down_DEGs.txt
Set2Genes=$OutDir/Frq53vsWT53_L_down2.txt
AllGenes=$OutDir/Frq53vsWT53_L_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

#Frq08_LD up and down regulated

```bash
WorkDir=analysis/enrichment/experiment_ALL
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08_LD/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08_LD_up_LFC_names.txt
AllGenes=$OutDir/Frq08_LD_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08_LD_up_DEGs.txt
Set2Genes=$OutDir/Frq08_LD_up2.txt
AllGenes=$OutDir/Frq08_LD_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08_LD/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08_LD_down_LFC_names.txt
AllGenes=$OutDir/Frq08_LD_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08_LD_down_DEGs.txt
Set2Genes=$OutDir/Frq08_LD_down2.txt
AllGenes=$OutDir/Frq08_LD_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Frq08/WT08 in D up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08vsWT08_D/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08vsWT08_D_up_LFC_names.txt
AllGenes=$OutDir/Frq08vsWT08_D_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08vsWT08_D_up_DEGs.txt
Set2Genes=$OutDir/Frq08vsWT08_D_up2.txt
AllGenes=$OutDir/Frq08vsWT08_D_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08vsWT08_D/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08vsWT08_D_down_LFC_names.txt
AllGenes=$OutDir/Frq08vsWT08_D_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08vsWT08_D_down_DEGs.txt
Set2Genes=$OutDir/Frq08vsWT08_D_down2.txt
AllGenes=$OutDir/Frq08vsWT08_D_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Frq08/WT08 in L up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08vsWT08_L/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08vsWT08_L_up_LFC_names.txt
AllGenes=$OutDir/Frq08vsWT08_L_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08vsWT08_L_up_DEGs.txt
Set2Genes=$OutDir/Frq08vsWT08_L_up2.txt
AllGenes=$OutDir/Frq08vsWT08_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_ALL/Frq08/Frq08vsWT08_L/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL/Frq08vsWT08_L_down_LFC_names.txt
AllGenes=$OutDir/Frq08vsWT08_L_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Frq08vsWT08_L_down_DEGs.txt
Set2Genes=$OutDir/Frq08vsWT08_L_down2.txt
AllGenes=$OutDir/Frq08vsWT08_L_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Draw venn diagrams of differenitally expressed genes

##WT53_LD vs  WC1_LD

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_up_LFC.txt
inp2=Wc1_LD_up_LFC.txt
OutDir=Wc1vsWT53_LD_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_down_LFC.txt
inp2=Wc1_LD_down_LFC.txt
OutDir=Wc1vsWT53_LD_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Wc1vsWT53_LD_up_LFC_DEGs.tsv --out $WorkDir/Wc1vsWT53_LD_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Wc1vsWT53_LD_down_LFC_DEGs.tsv --out $WorkDir/Wc1vsWT53_LD_down_LFC_DEGs.pdf
```

##WT53_LD vs  Frq53_LD

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_up_LFC.txt
inp2=Frq53_LD_up_LFC.txt
OutDir=Frq53vsWT53_LD_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_down_LFC.txt
inp2=Frq53_LD_down_LFC.txt
OutDir=Frq53vsWT53_LD_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsWT53_LD_up_LFC_DEGs.tsv --out $WorkDir/Frq53vsWT53_LD_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsWT53_LD_down_LFC_DEGs.tsv --out $WorkDir/Frq53vsWT53_LD_down_LFC_DEGs.pdf
```

##WT08_LD vs  Frq08_LD

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT08_LD_up_LFC.txt
inp2=Frq08_LD_up_LFC.txt
OutDir=Frq08vsWT08_LD_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT08_LD_down_LFC.txt
inp2=Frq08_LD_down_LFC.txt
OutDir=Frq08vsWT08_LD_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq08vsWT08_LD_up_LFC_DEGs.tsv --out $WorkDir/Frq08vsWT08_LD_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq08vsWT08_LD_down_LFC_DEGs.tsv --out $WorkDir/Frq08vsWT08_LD_down_LFC_DEGs.pdf
```

##Frq53/WT53 vs Frq08/WT08 in D

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=Frq53vsWT53_D_up_LFC.txt
inp2=Frq08vsWT08_D_up_LFC.txt
OutDir=Frq53vsFrq08_D_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=Frq53vsWT53_D_down_LFC.txt
inp2=Frq08vsWT08_D_down_LFC.txt
OutDir=Frq53vsFrq08_D_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsFrq08_D_up_LFC_DEGs.tsv --out $WorkDir/Frq53vsFrq08_D_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsFrq08_D_down_LFC_DEGs.tsv --out $WorkDir/Frq53vsFrq08_D_down_LFC_DEGs.pdf
```

##Frq53/WT53 vs Frq08/WT08 in L

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=Frq53vsWT53_L_up_LFC.txt
inp2=Frq08vsWT08_L_up_LFC.txt
OutDir=Frq53vsFrq08_L_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=Frq53vsWT53_L_down_LFC.txt
inp2=Frq08vsWT08_L_down_LFC.txt
OutDir=Frq53vsFrq08_L_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/ALL
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsFrq08_L_up_LFC_DEGs.tsv --out $WorkDir/Frq53vsFrq08_L_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/Frq53vsFrq08_L_down_LFC_DEGs.tsv --out $WorkDir/Frq53vsFrq08_L_down_LFC_DEGs.pdf
```

#Substract list of genes from countData_all

sed -i "s/.t1//g" TFs_names.csv

```R
library(DESeq2)
#make sure that TFs and countData have a first column named "genes"
countData <- read.table("countData_all",header=T,sep="\t")
TFs <- read.csv("TFs_names.csv")
results=merge(TFs, countData, all.x=T, by="genes")

results2= results[,c(1,11,12,13,14,15,16,31,32,33,34,35,36)]

#write.table(results2,"countData_Wc53_TFs.txt",sep="\t",quote=F)
```

```R
library("RColorBrewer")
library("gplots")
library( "genefilter" )

#Heatmap with ColSideColors
TFs<-data.matrix(results2)
my_palette <- colorRampPalette(c("green","green2", "black","red2", "red"))(n = 100)
heatmap.2( TFs, scale="row",
trace="none", dendrogram="row", margins = c(8, 20), col = my_palette, cexCol=0.8,cexRow=0.4)
dev.off()
```
