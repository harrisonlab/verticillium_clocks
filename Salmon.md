#Salmon tool for transcript quantification from RNA-seq data.
http://salmon.readthedocs.io/en/latest/salmon.html
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.

Note that it is designed to align to predicted transcripts and not to the whole genomeDir

V. dahliae transcripts


cat Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa | sed 's/cds.*//g' > Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa
#cat  | grep '>' | cut -f1 | sed 's/>//g'

```bash
for Transcriptome in $(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
echo "$Strain"
for RNADir in $(ls -d qc_rna/*/*/*); do
FileNum=$(ls $RNADir/F/*trim.fq.gz | wc -l)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
printf "\n"
FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f3 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
#Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
echo "$Timepoint"
OutDir=RNA_alignment/salmon/$Strain/$Prefix/$Timepoint
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p RNA_alignment/salmon/JR2/experiment_all

#This command creates a two column file with transcript_id and gene_id.
for File in $(ls RNA_alignment/salmon/JR2/experiment1/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > RNA_alignment/salmon/JR2/experiment_all/trans2gene.txt
done

# Put files in a convenient location for DeSeq.
#For experiment 1
for File in $(ls RNA_alignment/salmon/JR2/experiment1/*/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  mkdir -p RNA_alignment/salmon/JR2/experiment_all/$Prefix
  cp $PWD/$File RNA_alignment/salmon/JR2/experiment_all/$Prefix/quant.sf
done

#For experiment 2
for File in $(ls RNA_alignment/salmon/JR2/experiment2/*/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  mkdir -p RNA_alignment/salmon/JR2/experiment_all/$Prefix
  cp $PWD/$File RNA_alignment/salmon/JR2/experiment_all/$Prefix/quant.sf
done
```

#GENE EXPRESSION IN DESEQ2

This analysis was done repeating the salmon alignment with the option --keepduplicates.


/home/deakig/R3.4/bin/R

setwd("/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/salmon/JR2/experiment_all")

#===============================================================================
#       Load libraries (R 3.4 required)
#===============================================================================

library("DESeq2")
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
library(tximport)
library(rjson)
library(readr)
require("pheatmap")
require(data.table)
library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggrepel")

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# Create a file for all the samples in experiment 1 and experiment 2

# import transcript to gene mapping info
tx2gene <- read.table("trans2gene.txt",header=T,sep="\t")
head(tx2gene)

# import quantification files
txi <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene)
names(txi)
head(txi$counts)

# get the sample names from the folders
mysamples <- list.dirs("./",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.tx <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""), type = "salmon", txOut = TRUE)
txi.sum <- summarizeToGene(txi.tx,tx2gene)
head(txi.sum$counts)
all.equal(txi$counts, txi.sum$counts)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.sum[[i]])<<-mysamples}))

write.table(txi.sum$counts,"countData_all",sep="\t",na="",quote=F)

#===============================================================================
#       Eliminate samples from countData for different analysis
#===============================================================================

colData <- read.table("colData.txt",header=T,sep="\t")


#Eliminate sample from colData and countData

colData <- colData[!(colData$Sample=="Frq08_DD24_rep1"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD24_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep2"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD24_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD24_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep1"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD18_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep2"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD18_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep3"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD18_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep1"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD12_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep2"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD12_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep3"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD12_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep1"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD6_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep2"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD6_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep3"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_DD6_rep3)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep1"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_LL6_rep1)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep2"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_LL6_rep2)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep3"),]      
txi.sum <- subset(txi.sum$counts, select=-Frq08_LL6_rep3)



#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

ColData <- read.table("colData.txt",header=T,sep="\t")
rownames(ColData) <- colnames(txi.sum$counts)

# Group column
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

# 1st design
design <- ~ Timepoint + Cultivar
dds <- DESeqDataSetFromTximport(txi.reps,colData,design)

# Library normalisation
#dds <- estimateSizeFactors(dds)

# Set reference factor level
#dds$Cultivar<-factor(dds$Cultivar, levels=c("mycelium","GD","M9"))

# Deseq
dds<-DESeq(dds)
resultsNames(dds)
###
[1] "Intercept"          "Cultivar_M9_vs_GD"  "Timepoint_t2_vs_t1"
###
