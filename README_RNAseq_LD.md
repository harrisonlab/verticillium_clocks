http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
https://static-content.springer.com/esm/art%3A10.1007%2Fs13258-016-0479-2/MediaObjects/13258_2016_479_MOESM1_ESM.pdf

#Experiment all
```R
#Gene clustering of different group of samples
library(DESeq2)
colData <- read.table("colData_WT_LD",header=T,sep="\t")
countData <- read.table("countData_WT_LD",header=T,sep="\t")
colData$Group <- factor(paste0(colData$Strain,colData$Light))

design <- ~Group
dds <- 	DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- DESeq(dds)
rld <- rlog( dds, blind=FALSE )
resultsNames(dds)



=================
PCA plots
=================
library("ggplot2")
library("ggrepel")

data <- plotPCA(rld, intgroup=c("Strain", "Light"), returnData=TRUE)
data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, shape=Strain, color=Light)) +
 geom_point(size=8) +
 theme_minimal()+
 theme(axis.title.x = element_text(size=30),
 axis.title.y = element_text(size=30),
 legend.text=element_text(size=30),
 legend.title=element_text(size=30),
 axis.text.x = element_text(colour="black",size=25),
 axis.text.y = element_text(colour="black",size=25),)+
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) +
 geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("PCA_sample_names.pdf", pca_plot, dpi=300, height=15, width=20)

=================
Sample distanes
=================

rld <- rlog( dds )
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sample)
colnames(sampleDistMatrix) <- paste(rld$Sample)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins = c(10, 10))

dev.off()

===============
Analysis of gene expression
===============
https://support.bioconductor.org/p/84366/
https://support.bioconductor.org/p/62355/

## 53 L vs D (treatment effect for genotype WT53)
resultsNames(dds)
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "WT53l", "WT53d"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"WT53_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53_LD_down.txt",sep="\t",quote=F)
write.table(sig.res,"WT53_LD_sigres.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]

write.table(sig.res.upregulated,"WT53_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53_LD_down_LFC.txt",sep="\t",quote=F)


##08 L vs D
resultsNames(dds)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "WT08l", "WT08d"))
mcols(res, use.names=TRUE)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"WT08_LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT08_LD_down.txt",sep="\t",quote=F)
write.table(sig.res,"WT08_LD_sigres.txt",sep="\t",quote=F)
summary(res)

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange < -1, ]
write.table(sig.res.upregulated,"WT08_LD_up_LFC.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WT08_LD_down_LFC.txt",sep="\t",quote=F)

```
===============
Dispersion
==============
```R
plotDispEsts(dds)
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

===============
TopVarGenes
===============

```R
library("RColorBrewer")
library("gplots")
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)

#Heatmap with ColSideColors
StrainCols <- brewer.pal(11, "RdGy")[c(7, 11)]
my_palette <- colorRampPalette(c("orange", "lightgreen", "darkgreen"))(n = 299)
heatmap.2( assay(rld)[ topVarGenes,], scale="row",
trace="none", dendrogram="row", margins = c(5, 30), col = my_palette, cexCol=0.8,cexRow=0.4,ColSideColors=StrainCols[unclass(dds$Light)])
legend("topright", levels(dds$Light), col = StrainCols, lty = 1, lwd = 5, cex = 0.5)
dev.off()
```
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Light ]
mat <- assay(rld)[ topVarGenes, ],scale="row",
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$Strain,"-",rld$Light)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")


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
input=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT
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
python $scripts/calculate_fpkm.py countData_WT_LD Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt

```

#Produce a detailed table of analyses
'''python
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

Run with commands:
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
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_all/WT)
    DEG_Files=$(ls \
        $Dir1/WT08_LD.txt \
        $Dir1/WT53_LD.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts.txt)
    FPKM=$(ls $Dir1/countData_WT_LD.fpkm)
    #ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --TFs $TFs --InterPro $InterPro --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_WT_allgene.tsv
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

#WT53_LD up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_WT/LFC/WT53_LD/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT/WT53_LD_up_LFC_names.txt
AllGenes=$OutDir/WT53_LD_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/WT53_LD_up_DEGs.txt
Set2Genes=$OutDir/WT53_LD_up2.txt
AllGenes=$OutDir/WT53_LD_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_WT/LFC/WT53_LD/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT/WT53_LD_down_LFC_names.txt
AllGenes=$OutDir/WT53_LD_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/WT53_LD_down_DEGs.txt
Set2Genes=$OutDir/WT53_LD_down2.txt
AllGenes=$OutDir/WT53_LD_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

#WT08_LD up and down regulated
```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_WT/LFC/WT08_LD/UP
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT/WT08_LD_up_LFC_names.txt
AllGenes=$OutDir/WT08_LD_up_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/WT08_LD_up_DEGs.txt
Set2Genes=$OutDir/WT08_LD_up2.txt
AllGenes=$OutDir/WT08_LD_up.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

```bash
WorkDir=analysis/enrichment/experiment_WT
OutDir=analysis/enrichment/experiment_WT/LFC/WT08_LD/DOWN
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
AnnotTable=gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_3.tsv
DEGs=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT/WT08_LD_down_LFC_names.txt
AllGenes=$OutDir/WT08_LD_down_allgenes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/WT08_LD_down_DEGs.txt
Set2Genes=$OutDir/WT08_LD_down2.txt
AllGenes=$OutDir/WT08_LD_down.txt
cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
cat $AnnotTable | tail -n+2 | cut -f1 | cut -d'.' -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $WorkDir/experiment_all_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


#Draw venn diagrams of differenitally expressed genes

##All genes

###All DEGs
```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_sigres.txt
inp2=WT08_LD_sigres.txt
OutDir=WT53vsWT08_LD_all_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Upregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=WT53_LD_up_LFC.txt
inp2=WT08_LD_up_LFC.txt
OutDir=WT53vsWT08_LD_up_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
inp1=Wt53_LD_down_LFC.txt
inp2=WT08_LD_down_LFC.txt
OutDir=WT53vsWT08_LD_down_LFC_DEGs.tsv
$ProgDir/parse_RNASeq_2-way.py --input_1 $inp1 --input_2 $inp2 --out_file $OutDir
```

###Venn diagrams
#Mutants in light

```bash
ProgDir=/home/lopeze/git_repos/scripts/verticillium_clocks
WorkDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment_all/WT
$ProgDir/All_DEGs_venn_diag_2.r  --inp $WorkDir/WT53vsWT08_LD_all_DEGs.tsv --out $WorkDir/WT53vsWT08_LD_all_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/WT53vsWT08_LD_up_LFC_DEGs.tsv --out $WorkDir/WT53vsWT08_LD_up_LFC_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/WT53vsWT08_LD_down_LFC_DEGs.tsv --out $WorkDir/WT53vsWT08_LD_LFC_down_DEGs.pdf
