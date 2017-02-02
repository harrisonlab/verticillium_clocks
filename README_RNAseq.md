#V. dahliae RNAseq

In order to open and extract the files from the .tar document:
tar -xvf C101HW16120207.tar

##Building a directory structure
RawDatDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNAseq_data/C101HW16120207/raw_data
ProjectDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/raw_RNA/experiment1/V.dahliae/$Strain

To create folders F and R inside each strain folder:

```bash
ls > t
for folder in $(cat t); do
  echo $folder
  mkdir -p ./"$folder"/F
  mkdir -p ./"$folder"/R
done
rm t
```

Sequence data was moved to the correct directory

```bash
RawDatDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNAseq_data/C101HW16120207/raw_data
ProjectDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/raw_rna
    mv $RawDatDir/EC_1_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep1/F/.
    mv $RawDatDir/EC_1_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep1/R/.
    mv $RawDatDir/EC_2_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep2/F/.
    mv $RawDatDir/EC_2_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep2/R/.
    mv $RawDatDir/EC_3_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep3/F/.
    mv $RawDatDir/EC_3_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD18_rep3/R/.
    mv $RawDatDir/EC_4_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep1/F/.
    mv $RawDatDir/EC_4_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep1/R/.
    mv $RawDatDir/EC_5_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep2/F/.
    mv $RawDatDir/EC_5_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep2/R/.
    mv $RawDatDir/EC_6_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep3/F/.
    mv $RawDatDir/EC_6_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD18_rep3/R/.
    mv $RawDatDir/EC_7_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep1/F/.
    mv $RawDatDir/EC_7_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep1/R/.
    mv $RawDatDir/EC_8_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep2/F/.
    mv $RawDatDir/EC_8_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep2/R/.
    mv $RawDatDir/EC_9_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep3/F/.
    mv $RawDatDir/EC_9_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD24_rep3/R/.
    mv $RawDatDir/EC_10_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep1/F/.
    mv $RawDatDir/EC_10_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep1/R/.
    mv $RawDatDir/EC_11_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep2/F/.
    mv $RawDatDir/EC_11_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep2/R/.
    mv $RawDatDir/EC_12_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep3/F/.
    mv $RawDatDir/EC_12_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD24_rep3/R/.
    mv $RawDatDir/EC_13_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep1/F/.
    mv $RawDatDir/EC_13_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep1/R/.
    mv $RawDatDir/EC_14_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep2/F/.
    mv $RawDatDir/EC_14_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep2/R/.
    mv $RawDatDir/EC_15_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep3/F/.
    mv $RawDatDir/EC_15_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD6_rep3/R/.
    mv $RawDatDir/EC_16_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep1/F/.
    mv $RawDatDir/EC_16_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep1/R/.
    mv $RawDatDir/EC_17_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep2/F/.
    mv $RawDatDir/EC_17_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep2/R/.
    mv $RawDatDir/EC_18_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep3/F/.
    mv $RawDatDir/EC_18_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD6_rep3/R/.
    mv $RawDatDir/EC_19_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep1/F/.
    mv $RawDatDir/EC_19_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep1/R/.
    mv $RawDatDir/EC_20_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep2/F/.
    mv $RawDatDir/EC_20_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep2/R/.
    mv $RawDatDir/EC_21_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep3/F/.
    mv $RawDatDir/EC_21_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_DD6_rep3/R/.
    mv $RawDatDir/EC_22_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep1/F/.
    mv $RawDatDir/EC_22_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep1/R/.
    mv $RawDatDir/EC_23_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep2/F/.
    mv $RawDatDir/EC_23_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep2/R/.
    mv $RawDatDir/EC_24_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep3/F/.
    mv $RawDatDir/EC_24_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_DD12_rep3/R/.
    mv $RawDatDir/EC_25_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep1/F/.
    mv $RawDatDir/EC_25_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep1/R/.
    mv $RawDatDir/EC_26_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep2/F/.
    mv $RawDatDir/EC_26_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep2/R/.
    mv $RawDatDir/EC_27_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep3/F/.
    mv $RawDatDir/EC_27_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_DD12_rep3/R/.
    mv $RawDatDir/EC_28_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep1/F/.
    mv $RawDatDir/EC_28_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep1/R/.
    mv $RawDatDir/EC_29_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep2/F/.
    mv $RawDatDir/EC_29_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep2/R/.
    mv $RawDatDir/EC_30_1.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep3/F/.
    mv $RawDatDir/EC_30_2.fq.gz $ProjectDir/experiment1/V.dahliae/53WT_LL6_rep3/R/.
    mv $RawDatDir/EC_31_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep1/F/.
    mv $RawDatDir/EC_31_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep1/R/.
    mv $RawDatDir/EC_32_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep2/F/.
    mv $RawDatDir/EC_32_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep2/R/.
    mv $RawDatDir/EC_33_1.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep3/F/.
    mv $RawDatDir/EC_33_2.fq.gz $ProjectDir/experiment1/V.dahliae/Frq08_LL6_rep3/R/.
    mv $RawDatDir/EC_34_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep1/F/.
    mv $RawDatDir/EC_34_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep1/R/.
    mv $RawDatDir/EC_35_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep2/F/.
    mv $RawDatDir/EC_35_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep2/R/.
    mv $RawDatDir/EC_36_1.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep3/F/.
    mv $RawDatDir/EC_36_2.fq.gz $ProjectDir/experiment1/V.dahliae/Wc153_LL6_rep3/R/.
```

#Data QC
Perform qc of RNAseq timecourse data

```bash
screen -a
for FilePath in $(ls -d raw_rna/experiment1/V.*/*); do
#Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
echo $FilePath
FileF=$(ls $FilePath/F/*.fq.gz);
FileR=$(ls $FilePath/R/*.fq.gz);
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
sleep 10m
done

  ```

Data quality was visualised using fastqc:

```bash
for RawData in $(ls qc_rna/experiment1/V.dahliae/*/*/*.fq.gz); do
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
        echo $RawData;
        qsub $ProgDir/run_fastqc.sh $RawData
    done
```

##Trim data

Trimming was performed on data to trim adapters from sequences and remove poor quality data. This was done with fastq-mcf

Trimming was first performed on the strain that had a single run of data:

```bash
for StrainPath in $(ls -d raw_rna/experiment1/V.dahliae/*); do
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
        IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
        ReadsF=$(ls $StrainPath/F/*.fq*)
        ReadsR=$(ls $StrainPath/R/*.fq*)
        echo $ReadsF
        echo $ReadsR
        qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
        sleep 5m
    done
    ```

    Data quality was visualised once again following trimming:

```bash
      for RawData in $(ls qc_rna/experiment1/V.dahliae/*/*/*.fq.gz); do
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
        echo $RawData;
        qsub $ProgDir/run_fastqc.sh $RawData
      done
      ```

  ##Filter data

```bash
for Strain in $(ls ./* -d) ; do
echo $Strain
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/bowtie/experiment1/V.dahliae/"$Strain"/F/
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/bowtie/experiment1/V.dahliae/"$Strain"/R/
done
```

```bash
for FilePath in $(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/*); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
echo $Strain
F_Read=$(ls $FilePath/F/*.fq.gz)
R_Read=$(ls $FilePath/R/*.fq.gz)
echo $F_Read
echo $R_read
/home/groups/harrisonlab/project_files/quorn/scripts/bowtie.sh \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/"$Strain"/F \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/"$Strain"/R \
/home/groups/harrisonlab/project_files/quorn/filtered/phix/phix \ /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/analysis/bowtie/experiment1/V.dahliae/"$Strain"/ \ "$Strain" 150 300 \
done
```

The results showed 0 reads 0.00% overall alignment rate which means there's no Phix in the samples.



#Align to reference genome using STAR

The program was copied to:
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks


```bash
for Strain in $(ls ./* -d) ; do
echo $Strain
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/"$Strain"
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/"$Strain"
done
```

```bash
STAR
for FilePath in $(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/*); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
echo $Strain
F_Read=$(ls $FilePath/F/*.fq.gz)
R_Read=$(ls $FilePath/R/*.fq.gz)
echo $F_Read
echo $R_read
qsub sub_star.sh
./STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae \
--outFileNamePrefix /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1 \
--readFilesCommand zcat \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/F/*.fq.gz \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/R/*.fq.gz \
--outSAMtype SAM
done
```



 ##To create the index file for STAR:

Remove the spaces from the Chromosomes names and subtitute them with _

```cat Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa | cut -f1 -d ' ' > Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa
 cat Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa |sed 's/ /_/g' > Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa
```

#-----
#First attempt
#-----

JR2 Genome


```
./STAR \
--runMode genomeGene rate \
--genomeDir /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae \
--genomeFastaFiles /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3 \
--outFileNamePrefix /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/ \ --readFilesCommand zcat \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/F/*.fq.gz \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/R/*.fq.gz \
--outSAMtype SAM \
--runThreadN 16
```

12008 Genome

```
./STAR \
--runMode genomeGenerate \
--genomeDir /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae \
--genomeFastaFiles /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_unmasked.fa \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/final_genes/V.dahliae/12008/final/final_genes_appended.gff3 \
--outFileNamePrefix /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/ \
--readFilesCommand zcat \
--readFilesIn /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/F/*.fq.gz \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/R/*.fq.gz \
--outSAMtype SAM \
--runThreadN 16
```



#------
# 2nd attempt
#------

12008

### To create the Index File

#GenomeDir=$WorkDir/index
#InGenome=../pathogenomics/repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_unmasked.fa
#InGff=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/final_genes/V.dahliae/12008/final/final_genes_appended.gff3
# ParentFeature="ID"

./STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff
--runThreadN 16 \
--sjdbOverhang ReadLength-1

### To align reads with STAR
./STAR \
--genomeDir $GenomeDir
--outFileNamePrefix /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/ \
--readFilesCommand zcat \
--readFilesIn /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/F/*.fq.gz \
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/R/*.fq.gz \
--outSAMtype SAM \
--runThreadN 16


JR2

### To create the Index File

GenomeDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1
InGenome=public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa
InGff=public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3
ParentFeature="Parent"

./STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff
--runThreadN 16 \
--sjdbOverhang ReadLength-1


### To align reads with STAR

screen -a
qlogin -pe smp 16 -l virtual_free=1.1G

GenomeDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1
OutDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/53WT_DD12_rep1
InReadF=qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/F/*.fq.gz
InReadR=qc_rna/experiment1/V.dahliae/53WT_DD12_rep1/R/*.fq.gz

./STAR \
--genomeDir $GenomeDir \
--outFileNamePrefix $OutDir \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 16


##Create sub_star.sh

#!/bin/bash
#Align RNAseq data with genome using STAR

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=1.2G

# ---------------
# Step 1
# Collect inputs
# ---------------

InGenome=$(basename $1)
InGff=$(basename $2)
InReadF=$(basename $3)
InReadR=$(basename $4)
# GenomeDir=$(basname $5)
OutDir=$5

CurDir=$PWD
WorkDir=$TMPDIR/star
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $InGenome
cp $CurDir/$2 $InGff
cp $CurDir/$3 $InReadF
cp $CurDir/$4 $InReadR
# cp $CurDir/$5 $GenomeDir
cp $CurDir/$5 $OutDir

GenomeDir=$WorkDir/index

# ---------------
# Step 2
# Create the Index File
# ---------------

ParentFeature="Parent"
# ParentFeature="ID"

./STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff
--runThreadN 8 \
--sjdbOverhang ReadLength-1

# ---------------
# Step 2=3
# Run STAR
# ---------------

./STAR \
--genomeDir $GenomeDir \
--outFileNamePrefix $OutDir \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 8


rm -r $GenomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.


## Run sub_star.sh in data set

```bash
InGenome=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3)
OutDir=$(ls RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep2/53WT_DD12_rep2)
InReadF=$(ls qc_rna/experiment1/V.dahliae/53WT_DD12_rep2/F/*.fq.gz)
InReadR=$(ls qc_rna/experiment1/V.dahliae/53WT_DD12_rep2/R/*.fq.gz)
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $InGenome $InGff $InReadF $InReadR $OutDir
```

## Run sub_star.sh in a loop

```bash
for FilePath in $(ls -d qc_rna/experiment1/V.dahliae/*); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
InGenome=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3)
InReadF=$(ls qc_rna/experiment1/V.dahliae/$Strain/F/*.fq.gz)
InReadR=$(ls qc_rna/experiment1/V.dahliae/$Strain/R/*.fq.gz)
OutDir=RNA_alignment/STAR/experiment1/V.dahliae/$Strain
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
echo "qsub $ProgDir/sub_star.sh $InGenome $InGff $InReadF $InReadR $OutDir"
done
```
