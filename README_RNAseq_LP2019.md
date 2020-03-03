
# V. dahliae LP RNAseq experiment Dec2019

#To download the RNSseq data from a link using the command ling, use wget in the preferred directory:

```bash
wget link
```
#In order to open and extract the files from the .tar document:

```bash
tar -xvf X201SC19112783-Z01-F001_1_20200114_91chX0.tar
tar -xvf X201SC19112783-Z01-F001_2_20200114_CtM3Ta.tar
tar -xvf X201SC19112783-Z01-F002_20200206_y35fds.tar
```

#Copy the directory structure from the RawDatDir to the ProjectDir (only directories)

```bash
cd destination/directory
find /projects/vertclock/RNAseq_data/LP_experiment2019/X201SC19112783-Z01-F001_1_20200114/Rawdata  -type d -printf "%P\n" | xargs mkdir -p
```

#To create folders F and R inside each strain folder:

```bash
ls > t
for folder in $(cat t); do
	echo $folder
	mkdir -p ./"$folder"/F
	mkdir -p ./"$folder"/R
done
rm t
```
#T o create folders F and R inside each strain folder (simpler version):

```bash
for folder /Wc1_15_2/; do
	mkdir -p ./"$folder"/F
	mkdir -p ./"$folder"/R
done
```

#Move sequence data to the correct directory

```bash
RawDatDir=/projects/vertclock/RNAseq_data/LP_experiment2019/X201SC19112783-Z01-F001_1_20200114/Rawdata
ProjectDir=/projects/vertclock/raw_rna/V.dahliae/LP_experiment2019

mv $RawDatDir/W53_15_1/W53_15_1_1.fq.gz $ProjectDir/W53_15_1/F/.
mv $RawDatDir/W53_15_1/W53_15_1_2.fq.gz $ProjectDir/W53_15_1/R/.
mv $RawDatDir/W53_15_2/W53_15_2_1.fq.gz $ProjectDir/W53_15_2/F/.
mv $RawDatDir/W53_15_2/W53_15_2_2.fq.gz $ProjectDir/W53_15_2/R/.
mv $RawDatDir/W53_15_3/W53_15_3_1.fq.gz $ProjectDir/W53_15_3/F/.
mv $RawDatDir/W53_15_3/W53_15_3_2.fq.gz $ProjectDir/W53_15_3/R/.
mv $RawDatDir/W53_30_1/W53_30_1_1.fq.gz $ProjectDir/W53_30_1/F/.
mv $RawDatDir/W53_30_1/W53_30_1_2.fq.gz $ProjectDir/W53_30_1/R/.
mv $RawDatDir/W53_30_2/W53_30_2_1.fq.gz $ProjectDir/W53_30_2/F/.
mv $RawDatDir/W53_30_2/W53_30_2_2.fq.gz $ProjectDir/W53_30_2/R/.
mv $RawDatDir/W53_30_3/W53_30_3_1.fq.gz $ProjectDir/W53_30_3/F/.
mv $RawDatDir/W53_30_3/W53_30_3_2.fq.gz $ProjectDir/W53_30_3/R/.
mv $RawDatDir/W53_D_1/W53_D_1_1.fq.gz $ProjectDir/W53_D_1/F/.
mv $RawDatDir/W53_D_1/W53_D_1_2.fq.gz $ProjectDir/W53_D_1/R/.
mv $RawDatDir/W53_D_2/W53_D_2_1.fq.gz $ProjectDir/W53_D_2/F/.
mv $RawDatDir/W53_D_2/W53_D_2_2.fq.gz $ProjectDir/W53_D_2/R/.
mv $RawDatDir/W53_D_3/W53_D_3_1.fq.gz $ProjectDir/W53_D_3/F/.
mv $RawDatDir/W53_D_3/W53_D_3_2.fq.gz $ProjectDir/W53_D_3/R/.
mv $RawDatDir/Wc1_15_1/Wc1_15_1_1.fq.gz $ProjectDir/Wc1_15_1/F/.
mv $RawDatDir/Wc1_15_1/Wc1_15_1_2.fq.gz $ProjectDir/Wc1_15_1/R/.
mv $RawDatDir/Wc1_30_1/Wc1_30_1_1.fq.gz $ProjectDir/Wc1_30_1/F/.
mv $RawDatDir/Wc1_30_1/Wc1_30_1_2.fq.gz $ProjectDir/Wc1_30_1/R/.
mv $RawDatDir/Wc1_30_3/Wc1_30_3_1.fq.gz $ProjectDir/Wc1_30_3/F/.
mv $RawDatDir/Wc1_30_3/Wc1_30_3_2.fq.gz $ProjectDir/Wc1_30_3/R/.
mv $RawDatDir/Wc1_D_3/Wc1_D_3_1.fq.gz $ProjectDir/Wc1_D_3/F/.
mv $RawDatDir/Wc1_D_3/Wc1_D_3_2.fq.gz $ProjectDir/Wc1_D_3/R/.

RawDatDir=/projects/vertclock/RNAseq_data/LP_experiment2019/X201SC19112783-Z01-F002/raw_data
ProjectDir=/projects/vertclock/raw_rna/V.dahliae/LP_experiment2019
mv $RawDatDir/Wc1_15_2_1.fq.gz $ProjectDir/Wc1_15_2/F/.
mv $RawDatDir/Wc1_15_2_2.fq.gz $ProjectDir/Wc1_15_2/R/.
mv $RawDatDir/Wc1_15_3_1.fq.gz $ProjectDir/Wc1_15_3/F/.
mv $RawDatDir/Wc1_15_3_2.fq.gz $ProjectDir/Wc1_15_3/R/.
mv $RawDatDir/Wc1_30_2_1.fq.gz $ProjectDir/Wc1_30_2/F/.
mv $RawDatDir/Wc1_30_2_2.fq.gz $ProjectDir/Wc1_30_2/R/.
mv $RawDatDir/Wc1_D_1_1.fq.gz $ProjectDir/Wc1_D_1/F/.
mv $RawDatDir/Wc1_D_1_2.fq.gz $ProjectDir/Wc1_D_1/R/.
mv $RawDatDir/Wc1_D_2_1.fq.gz $ProjectDir/Wc1_D_2/F/.
mv $RawDatDir/Wc1_D_2_2.fq.gz $ProjectDir/Wc1_D_2/R/.
```


## Data QC Perform qc of RNAseq timecourse data

```bash
screen -a
for FilePath in $(ls -d raw_rna/V.dahliae/LP_experiment2019/*); do
#Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
echo $FilePath
FileF=$(ls $FilePath/F/*.fq.gz);
FileR=$(ls $FilePath/R/*.fq.gz);
IlluminaAdapters=/projects/vertclock/git_repos/tools/seq_tools/rna_qc/ncbi_adapters.fa
ProgDir=/projects/vertclock/git_repos/tools/seq_tools/rna_qc/
sbatch $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
#sleep 10m
done
```

#Data quality was visualised using fastqc:

```bash
for RawData in $(ls qc_rna/V.dahliae/LP_experiment2019/*/*/*.fq.gz); do
				ProgDir=/projects/vertclock/git_repos/tools/seq_tools/dna_qc
				echo $RawData;
				sbatch $ProgDir/run_fastqc.sh $RawData
		done
```

## Trim data

Trimming was performed on data to trim adapters from sequences and remove poor quality data. This was done with fastq-mcf
Trimming was first performed on the strain that had a single run of data:

```bash
for StrainPath in $(ls -d raw_rna/V.dahliae/LP_experiment2019/*); do
				ProgDir=/projects/vertclock/git_repos/tools/seq_tools/rna_qc/
				IlluminaAdapters=/projects/vertclock/git_repos/tools/seq_tools/rna_qc/ncbi_adapters.fa
				ReadsF=$(ls $StrainPath/F/*.fq*)
				ReadsR=$(ls $StrainPath/R/*.fq*)
				echo $ReadsF
				echo $ReadsR
				sbatch $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
				#sleep 5m
		done
```

#Data quality was visualised once again following trimming:
 ```bash
for RawData in $(ls qc_rna/V.dahliae/LP_experiment2019/Wc1_D_3/R/*.fq.gz); do
				ProgDir=/projects/vertclock/git_repos/tools/seq_tools/dna_qc
				echo $RawData;
				sbatch $ProgDir/run_fastqc.sh $RawData
		done
```

## bbduk

Filter out rRNA Data using BBTools BBduk: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/

```bash
for StrainPath in $(ls -d qc_rna/V.dahliae/LP_experiment2019/*); do
    Strain=$(sed 's/.*\///' <<< $StrainPath)
    ProgDir=/projects/oldhome/deakig/pipelines/common/bbtools
    in1=$(ls $StrainPath/F/*trim.fq.gz)
    in2=$(ls $StrainPath/R/*trim.fq.gz)
    ref=/projects/oldhome/deakig/pipelines/common/resources/contaminants/ribokmers.fa.gz
		echo $StrainPath;
    sbatch -p long --mem 60000 -c 10 ./git_repos/tools/seq_tools/bbtools/duck_wrapper2.sh $ref $StrainPath $in1 $in2 $ProgDir $Strain
done
```

## Set up Bioconda environment

To install the new version of Salmon we must set up a Bioconda environment.
Follow: http://149.155.34.104/bioconda.html

#Download and install Anaconda2

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```
PATH=/home/lopeze/miniconda3/bin:$PATH

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#Creating a Conda environment called SALMON

```bash
conda create --name salmon
conda activate salmon
```

#Install salmon https://bioconda.github.io/recipes/salmon/README.html

```bash
conda install salmon
```

#have salmon print its help message to check it runs

```bash
salmon
```

#deactivate the environment once you are done with it

```bash
conda deactivate
```

## SALMON

Salmon tool for transcript quantification from RNA-seq data.
http://salmon.readthedocs.io/en/latest/salmon.html
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
Note that it is designed to align to predicted transcripts and not to the whole genomeDir.
Quantifying in mapping-based mode

#V. dahliae transcripts

```bash
cat Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa | sed 's/cds.*//g' > Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa
#cat  | grep '>' | cut -f1 | sed 's/>//g'
```

<!-- ```bash
for Transcriptome in $(ls /projects/oldhome/groups/harrisonlab/project_files/verticillium_dahliae/clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
echo "$Strain"
for RNADir in $(ls -d qc_rna/V.dahliae/LP_experiment2019/W53_15_1); do
FileNum=$(ls $RNADir/F/*cleaned.fq.gz | wc -l)
#Jobs=$(sbatch | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
printf "\n"
FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Experiment=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
#Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
echo "$Timepoint"
OutDir=RNA_alignment/salmon/$Experiment/$Strain/$Timepoint
ProgDir=/projects/vertclock/git_repos/tools/seq_tools/RNAseq
sbatch $ProgDir/sub_salmon_1.1.sh $Transcriptome $FileF $FileR $OutDir
done
done
done
``` -->

Activate salmon environment

```bash
conda activate salmon
```

```bash
salmon index \
-t /projects/oldhome/groups/harrisonlab/project_files/verticillium_dahliae/clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa \
-i transcripts_index \
--keepDuplicates -k 27

salmon quant \
-i transcripts_index \
-l IU \
-1 /projects/vertclock/qc_rna/V.dahliae/LP_experiment2019/Wc1_D_3/F/Wc1_D_3_1_cleaned.fq.gz \
-2 /projects/vertclock/qc_rna/V.dahliae/LP_experiment2019/Wc1_D_3/R/Wc1_D_3_2_cleaned.fq.gz \
--validateMappings
-p 4
--numBootstraps 1000
--dumpEq
--seqBias
--gcBias \
-o /projects/vertclock/RNA_alignment/salmon/LP_experiment2019/JR2/Wc1_D_3/transcripts_quant
```

#Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p projects/vertclock/RNA_alignment/salmon/LP_experiment2019/DeSeq2
```

#This command creates a two column file with transcript_id and gene_id.

```bash
for File in $(ls RNA_alignment/salmon/LP_experiment2019/JR2/*/transcripts_quant/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > RNA_alignment/salmon/LP_experiment2019/DeSeq2/trans2gene.txt
done
```
#Put files in a convenient location for DeSeq.

```bash
for File in $(ls RNA_alignment/salmon/LP_experiment2019/JR2/*/transcripts_quant/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  mkdir -p RNA_alignment/salmon/LP_experiment2019/DeSeq2/$Prefix
  cp $PWD/$File RNA_alignment/salmon/LP_experiment2019/DeSeq2/$Prefix/quant.sf
done
```

# GENE EXPRESSION IN DESEQ2

#This analysis was done repeating the salmon alignment with the option --keepduplicates.

```bash
setwd("/projects/vertclock/RNA_alignment/salmon/LP_experiment2019/DeSeq2")
```
#Launch R

```bash
/projects/software/R-3.6.1/bin/R
```
```R
#===============================================================================
#       Load libraries (R 3.4 required)
#===============================================================================

#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("Biostrings")
BiocManager::install("gplots")

# Libraries
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(12))
library("ggplot2")
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

# Create a file for all the samples
# import transcript to gene mapping info following https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
tx2gene <- read.table("trans2gene.txt",header=T,sep="\t")
head(tx2gene)

# import quantification files
txi <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene)
names(txi)
head(txi$counts)

# get the sample names from the folders
mysamples <- list.dirs("./",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
#txi.tx <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""), type = "salmon", txOut = TRUE)

#txi.sum <- summarizeToGene(txi.tx,tx2gene)
#head(txi.sum$counts)
#all.equal(txi$counts, txi.sum$counts)

#set the sample names for txi.genes
#invisible(sapply(seq(1,3), function(i) {colnames(txi.tx[[i]])<<-mysamples}))
invisible(sapply(seq(1,3), function(i) {colnames(txi[[i]])<<-mysamples}))

#write.table(txi.sum$counts,"countData_all",sep="\t",na="",quote=F)
```

```R
#===============================================================================
# DEseq2 Analysis
#===============================================================================
#Data is unordered as it is read in. This means data must be set into the same order as the samples were read into mysamples before integrating metadata and and read counts
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#countData <- read.table("countData",header=T,sep="\t")
colData <- read.table("colData",header=T,sep="\t")
colData$Group <- factor(paste0(colData$Strain, '_', colData$Time))

#Design
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi,colData,design)

#To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us. In a typical RNA-seq analysis this step is automatically performed by the DESeq() function
dds <- DESeq(dds)
resultsNames(dds)

#Transform values
rld <- rlog( dds, blind=FALSE )
vsd <- vst(dds, blind=FALSE)


#===============================================================================
#PCA plots
#===============================================================================
data <- plotPCA(vsd, intgroup="Strain", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Strain, shape=Strain)) +
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
 geom_text_repel(aes(label=colnames(vsd)))
 coord_fixed()
ggsave("PCA_sample_names_elipses.pdf", pca_plot, dpi=300, height=15, width=20)


#===============================================================================
#Sample distances
#===============================================================================
sampleDists <- dist( t( assay(vsd) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$Sample)
colnames(sampleDistMatrix) <- paste(vsd$Sample)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins = c(10, 10))

dev.off()

#===============================================================================
#Heatmap of the count matrix
#===============================================================================
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Strain","Time")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
				 cluster_cols=FALSE, annotation_col=df)

#===============================================================================
#Analysis of gene expression
#===============================================================================
#Relevel to W53 D
dds$Time <- relevel(dds$Time, "d")
dds$Strain <- relevel(dds$Strain, "W53")

#contrasts
resultsNames(dds)

## W53 15m vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "W53_15m", "W53_d"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"W53_15m_vs_d.txt",sep="\t",quote=F)


## W53 30m vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "W53_30m", "W53_d"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"W53_30m_vs_d.txt",sep="\t",quote=F)

## Wc1 15m vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc1_15m", "Wc1_d"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"Wc1_15m_vs_d.txt",sep="\t",quote=F)

## Wc1 30m vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc1_30m", "Wc1_d"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"Wc1_30m_vs_d.txt",sep="\t",quote=F)

## Wc1 D vs W53 D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc1_d", "W53_d"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"Wc1_d_vs_W53_d.txt",sep="\t",quote=F)

sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]


## Wc1 15m vs W53 15m
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc1_15m", "W53_15m"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"Wc1_15m_vs_W53_15m.txt",sep="\t",quote=F)

## Wc1 30m vs W53 30m
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group", "Wc1_30m", "W53_30m"))
mcols(res, use.names=TRUE)
summary(res)

write.table(res,"Wc1_30m_vs_W53_30m.txt",sep="\t",quote=F)
