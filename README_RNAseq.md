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

GenomeDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep3
OutDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep3/53WT_DD12_rep3
InReadF=qc_rna/experiment1/V.dahliae/53WT_DD12_rep3/F/*.fq.gz
InReadR=qc_rna/experiment1/V.dahliae/53WT_DD12_rep3/R/*.fq.gz

./STAR \
--genomeDir $GenomeDir \
--outFileNamePrefix $OutDir \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 16


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
echo "Hello"
InGenome=$(basename $1)
InGff=$(basename $2)
InReadF=$(basename $3)
InReadR=$(basename $4)
# genomeDir=$(basname $5)
OutDir=$5

CurDir=$PWD
WorkDir=$TMPDIR/star
genomeDir=$WorkDir/index
echo "$WorkDir"
echo "$genomeDir"
mkdir -p $genomeDir
ls $genomeDir
cd $WorkDir

cp $CurDir/$1 $InGenome
cp $CurDir/$2 $InGff
cp $CurDir/$3 $InReadF
cp $CurDir/$4 $InReadR
# cp $CurDir/$5 $GenomeDir
cp $CurDir/$5 $OutDir



# ---------------
# Step 2
# Create the Index File
# ---------------
echo "Building index file"
ParentFeature="Parent"
# ParentFeature="ID"

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/STAR \
--runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff \
--runThreadN 8 \
--sjdbOverhang 99

# ---------------
# Step 2=3
# Run STAR
# ---------------

echo "Aligning RNAseq reads"

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/STAR \
--genomeDir $genomeDir \
--outFileNamePrefix star_aligment \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 8


rm -r $genomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.



## Run sub_star.sh in data set

```bash
InGenome=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
InReadF=$(ls qc_rna/experiment1/V.dahliae/53WT_DD18_rep1/F/*.fq.gz)
InReadR=$(ls qc_rna/experiment1/V.dahliae/53WT_DD18_rep1/R/*.fq.gz)
OutDir=RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD18_rep1/
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $InGenome $InGff $InReadF $InReadR $OutDir
```

## Run sub_star.sh in a loop
q
```bash
for FilePath in $(ls -d qc_rna/experiment1/V.dahliae/*); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
InGenome=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
InReadF=$(ls qc_rna/experiment1/V.dahliae/$Strain/F/*.fq.gz)
InReadR=$(ls qc_rna/experiment1/V.dahliae/$Strain/R/*.fq.gz)
OutDir=RNA_alignment/STAR/experiment1/V.dahliae/$Strain/$Strain
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $InGenome $InGff $InReadF $InReadR $OutDir
done
```


#FeatureCounts

Is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads.



```bash
for Strain in $(ls ./* -d) ; do
echo $Strain
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment1/V.dahliae/"$Strain"
mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/RNA_alignment/featureCounts/experiment1/V.dahliae/"$Strain"
done
```


```bash
InBam=$(ls RNA_alignment/STAR/experiment1/V.dahliae/Wc153_LL6_rep3/star/*.sam)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
#InGenome=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
OutDir=RNA_alignment/featureCounts/experiment1/V.dahliae/Wc153_LL6_rep3/featureCounts/Wc153_LL6_rep3_gene_featurecounts.txt
./subread-1.5.1-source/bin/featureCounts -p -B -M -R -a $InGff -t gene -g "gene_id" -o $OutDir $InBam

```

Status  RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/star/star_aligmentAligned.out.sam
Assigned        34203674
Unassigned_Ambiguity    10911560
Unassigned_MultiMapping 0
Unassigned_NoFeatures   1353797
Unassigned_Unmapped     35332
Unassigned_MappingQuality       0
Unassigned_FragmentLength       0
Unassigned_Chimera      0
Unassigned_Secondary    0
Unassigned_Nonjunction  0
Unassigned_Duplicate    0


Write .sh file for featureCounts program:

```bash
#!/bin/bash
#Align RNAseq data with genome using featureCounts

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1.2G

# ---------------
# Step 1
# Collect inputs
# ---------------

InBam=$(basename $1)
InGff=$(basename $2)
OutDir=$3
Prefix=$4

CurDir=$PWD
WorkDir=$TMPDIR/featureCounts
echo "$WorkDir"
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $InBam
cp $CurDir/$2 $InGff
cp $CurDir/$3 $OutDir


# ---------------
# Step 2
# Run featureCounts
# ---------------

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/subread-1.5.1-source/bin/featureCounts \
-p -B -M -R -a $InGff -t gene -g "gene_id" -o "$Prefix"_gene_featurecounts.txt $InBam


rm $InBam
rm $InGff
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.
```

## Run sub_featureCounts.sh in data set

```bash
InBam=$(ls RNA_alignment/STAR/experiment1/V.dahliae/53WT_DD12_rep1/star/*.sam)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
OutDir=RNA_alignment/featureCounts/experiment1/V.dahliae/53WT_DD12_rep1/
Prefix=53WT_DD12_rep1
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_featureCounts.sh $InBam $InGff $OutDir $Prefix
```

## Run sub_featureCounts.sh in a loop

```bash
for FilePath in $(ls -d qc_rna/experiment1/V.dahliae/*); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
InBam=$(ls RNA_alignment/STAR/experiment1/V.dahliae/$Strain/star/*.sam)
InGff=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
OutDir=RNA_alignment/featureCounts/experiment1/V.dahliae/$Strain/
Prefix=$Strain
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_featureCounts.sh $InBam $InGff $OutDir $Prefix
done
```



#Import data from featureCounts into R

Create a file with different columns about the experiment treatments:

Sample  Strain  Time    Light
53WT_LL6_rep1   53WT    6h      l
53WT_LL6_rep2   53WT    6h      l
53WT_LL6_rep3   53WT    6h      l
Frq08_LL6_rep1  Frq08   6h      l
Frq08_LL6_rep2  Frq08   6h      l
Frq08_LL6_rep3  Frq08   6h      l
Wc153_LL6_rep1  Wc153   6h      l
Wc153_LL6_rep2  Wc153   6h      l
Wc153_LL6_rep3  Wc153   6h      l
53WT_DD6_rep1   53WT    6h      d
53WT_DD6_rep2   53WT    6h      d
53WT_DD6_rep3   53WT    6h      d
Frq08_DD6_rep1  Frq08   6h      d
Frq08_DD6_rep2  Frq08   6h      d
Frq08_DD6_rep3  Frq08   6h      d
Wc153_DD6_rep1  Wc153   6h      d
Wc153_DD6_rep2  Wc153   6h      d
Wc153_DD6_rep3  Wc153   6h      d
53WT_DD12_rep1  53WT    12h     d
53WT_DD12_rep2  53WT    12h     d
53WT_DD12_rep3  53WT    12h     d
Frq08_DD12_rep1 Frq08   12h     d
Frq08_DD12_rep2 Frq08   12h     d
Frq08_DD12_rep3 Frq08   12h     d
53WT_DD18_rep1  53WT    DD18h   d
53WT_DD18_rep2  53WT    DD18h   d
53WT_DD18_rep3  53WT    DD18h   d
Frq08_DD18_rep1 Frq08   DD18h   d
Frq08_DD18_rep2 Frq08   DD18h   d
Frq08_DD18_rep3 Frq08   DD18h   d
53WT_DD24_rep1  53WT    24h     d
53WT_DD24_rep2  53WT    24h     d
53WT_DD24_rep3  53WT    24h     d
Frq08_DD24_rep1 Frq08   24h     d
Frq08_DD24_rep2 Frq08   24h     d
Frq08_DD24_rep3 Frq08   24h     d


R script to import and merge data into a formad good for DESeq2

install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

```R
#install and load libraries             
require(data.table)

#load tables into a "list of lists"
qq <- lapply(list.files(".",".*_gene_featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

#rename the samples column
lapply(qq,function(x) {names(x)[7]<-sub(".*\\/","",sub("\\/star.*","",names(x)[7]))})

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#output countData
write.table(countData,"countData",sep="\t",na="",quote=F)
#output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F)
```

colnames(countData) <- sub("X","",colnames(countData))
countData <- countData[,colData$Sample]

===========================
#To run DESeq2
===========================
```R
require(DESeq2)

colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

design <- ~Group

dds <- 	DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- DESeq(dds, fitType="local")


#Run DESeq2 removing an outlier

library(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")

colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)
#Eliminate Frq08_DD24_rep3 sample from colData and countData
colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
countData <- subset(countData, select=-Frq08_DD24_rep3)

design <- ~Group
dds <-  DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- DESeq(dds, fitType="local")
```

=================
Sample distanes
=================
```
vst<-varianceStabilizingTransformation(dds)
sampleDists<-dist(t(assay(vst)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap(sampleDistMatrix)

Sample distances measured with rlog transformation:

rld <- rlog( dds )
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours)

dev.off()
```
=================
PCA plots
=================
```
vst<-varianceStabilizingTransformation(dds)
plotPCA(vst,intgroup="Group")

Plot using rlog transformation:
plotPCA(rld,intgroup="Group")
plotPCA(rld,intgroup=c("Time","Group")),col=cols)


#Plot using rlog transformation, showing sample names:
plotPCA(rld, intgroup="Group")

library("ggplot2")
library("ggrepel")

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("PCA_sample_names.pdf", pca_plot, dpi=300, height=15, width=20)
```
===============
Analysis of gene expression
===============
```
#Example:
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08l6h","Frq08d6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]


##53WT L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTl6h","53WTd6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

summary(res)
write.table(sig.res.upregulated,"WtLDtxt_up",sep="\t",quote=F)
write.table(sig.res.downregulated,"WtLDtxt_down",sep="\t",quote=F)
write.table(sig.res.,"WtLDtxt",sep="\t",quote=F)

out of 3569 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1799, 50%
LFC < 0 (down)   : 1770, 50%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 8)


##Frq08 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08l6h","Frq08d6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08LD_down.txt",sep="\t",quote=F)

write.table(sig.res,"Frq08LD.txt",sep="\t",quote=F)

out of 3299 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1675, 51%
LFC < 0 (down)   : 1624, 49%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)


##Wc153 L vs D
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Wc153l6h","Wc153d6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"WC153LD_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"WC153LD_down.txt",sep="\t",quote=F)

write.table(res,"WC153LD.txt",sep="\t",quote=F)

out of 3512 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1756, 50%
LFC < 0 (down)   : 1756, 50%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 4)


##Light 53WT vs Wc153
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTl6h","Wc153l6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"LightWt53Wc153_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"LightWt53Wc153_down.text",sep="\t",quote=F)

write.table(res,"LightWt53Wc153.txt",sep="\t",quote=F)

out of 4793 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2391, 50%
LFC < 0 (down)   : 2402, 50%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)


##Dark 53WT vs Wc153
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTd6h","Wc153d6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"DarkWt53Wc153_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"DarkWt53Wc153_down.txt",sep="\t",quote=F)

write.table(res,"DarkWt53Wc153.txt",sep="\t",quote=F)

out of 4004 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2026, 51%
LFC < 0 (down)   : 1978, 49%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)


##53Wt 6h vs 12h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTd6h","53WTd12h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Wt53h6h12_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53h6h12_down.txt",sep="\t",quote=F)

write.table(res,"Wt53h6h12.txt",sep="\t",quote=F)

out of 138 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 102, 74%
LFC < 0 (down)   : 36, 26%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 195)

##53Wt 12h vs 18
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTd12h","53WTdDD18h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Wt53h12h18_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53h12h18_down.txt",sep="\t",quote=F)

write.table(res,"Wt53h12h18.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 306, 2.7%
LFC < 0 (down)   : 425, 3.7%
outliers [1]     : 3, 0.026%
low counts [2]   : 884, 7.7%
(mean count < 13)

##53Wt 18h vs 24h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTdDD18h","53WTd24h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Wt53h18h24_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53h18h24_down.txt",sep="\t",quote=F)

write.table(res,"Wt53h18h24.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2, 0.018%
LFC < 0 (down)   : 12, 0.11%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)

##Frq08 6h vs 12h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08d6h","Frq08d12h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h6h12_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h6h12_down.txt",sep="\t",quote=F)

write.table(res,"Frq08h6h12.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 288, 2.5%
LFC < 0 (down)   : 142, 1.2%
outliers [1]     : 3, 0.026%
low counts [2]   : 2431, 21%
(mean count < 96)

##Frq08 12h vs 18h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08d12h","Frq08dDD18h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h12h18_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h12h18_down.txt",sep="\t",quote=F)

write.table(res,"Frq08h12h18.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 71, 0.62%
LFC < 0 (down)   : 98, 0.86%
outliers [1]     : 3, 0.026%
low counts [2]   : 2873, 25%
(mean count < 141)

##Frq08 18h vs 24h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08dDD18h","Frq08d24h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h18h24_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h18h24_down.txt",sep="\t",quote=F)
write.table(res,"Frq08h18h24.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 20, 0.18%
LFC < 0 (down)   : 9, 0.079%
outliers [1]     : 3, 0.026%
low counts [2]   : 1325, 12%
(mean count < 26)

##53Wt 12h vs 24h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTd12h","53WTd24h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Wt53h12h24_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53h12h24_down.txt",sep="\t",quote=F)
write.table(res,"Wt53h12h24.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 0, 0%
LFC < 0 (down)   : 5, 0.044%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)

##Fr08 12h vs 24h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08d12h","Frq08d24h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h12h24_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h12h24_down.txt",sep="\t",quote=F)
write.table(res,"Frq08h12h24.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 30, 0.26%
LFC < 0 (down)   : 8, 0.07%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)

##53Wt 6h vs 18h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","53WTd6h","53WTdDD18h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Wt53h6h18_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Wt53h6h18_down.txt",sep="\t",quote=F)
write.table(res,"Wt53h6h18.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 59, 0.52%
LFC < 0 (down)   : 47, 0.41%
outliers [1]     : 3, 0.026%
low counts [2]   : 884, 7.7%
(mean count < 13)

##Fr08 6h vs 18h
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08d6h","Frq08dDD18h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h6h18_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h6h18_down.txt",sep="\t",quote=F)
write.table(res,"Frq08h6h18.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2, 0.018%
LFC < 0 (down)   : 2, 0.018%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)

##Fr08 6hL vs 53WT 6hL
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08l6h","53WTl6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h6lWT53h6l_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h6lWT53h6l_down.txt",sep="\t",quote=F)
write.table(res,"Frq08h6lWT53h6l.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2198, 19%
LFC < 0 (down)   : 1987, 17%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)

##Fr08 6hD vs 53WT 6hD
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Frq08d6h","53WTd6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res.upregulated,"Frq08h6dWT53h6d_up.txt",sep="\t",quote=F)
write.table(sig.res.downregulated,"Frq08h6dWT53h6d_down.txt",sep="\t",quote=F)
write.table(res,"Frq08h6dWT53h6d.txt",sep="\t",quote=F)

out of 11409 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2569, 23%
LFC < 0 (down)   : 2498, 22%
outliers [1]     : 3, 0.026%
low counts [2]   : 0, 0%
(mean count < 0)
```

==================
Gene clustering
==================
```
rld <- rlog( dds )
head( assay(rld) )

library("RColorBrewer")
library("gplots")
library( "genefilter" )

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
heatmap.2( assay(rld)[ topVarGenes,], scale="row", par(lend = 1),
trace="none", dendrogram="column", margins = c(5, 10), col = my_palette, cexCol=0.8,cexRow=0.8

#Gene clustering of different group of samples
library(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")

colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)
#Eliminate Frq08_DD24_rep3 sample from colData and countData
colData <- colData[!(colData$Sample=="Frq08_DD24_rep1"),]      
countData <- subset(countData, select=-Frq08_DD24_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep2"),]      
countData <- subset(countData, select=-Frq08_DD24_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
countData <- subset(countData, select=-Frq08_DD24_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep1"),]      
countData <- subset(countData, select=-Frq08_DD18_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep2"),]      
countData <- subset(countData, select=-Frq08_DD18_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD18_rep3"),]      
countData <- subset(countData, select=-Frq08_DD18_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep1"),]      
countData <- subset(countData, select=-Frq08_DD12_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep2"),]      
countData <- subset(countData, select=-Frq08_DD12_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD12_rep3"),]      
countData <- subset(countData, select=-Frq08_DD12_rep3)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep1"),]      
countData <- subset(countData, select=-Frq08_DD6_rep1)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep2"),]      
countData <- subset(countData, select=-Frq08_DD6_rep2)
colData <- colData[!(colData$Sample=="Frq08_DD6_rep3"),]      
countData <- subset(countData, select=-Frq08_DD6_rep3)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep1"),]      
countData <- subset(countData, select=-Frq08_LL6_rep1)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep2"),]      
countData <- subset(countData, select=-Frq08_LL6_rep2)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep3"),]      
countData <- subset(countData, select=-Frq08_LL6_rep3)

design <- ~Group
dds <-  DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- DESeq(dds, fitType="local")

rld <- rlog( dds )
head( assay(rld) )

library("RColorBrewer")
library("gplots")
library( "genefilter" )

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
heatmap.2( assay(rld)[ topVarGenes,], scale="row", par(lend = 1),
trace="none", dendrogram="column", margins = c(5, 10), col = my_palette, cexCol=0.8,cexRow=0.8

```
===============
Plotting results  
===============
```
##Light 53WT vs Wc153

res= results(dds, alpha=alpha,contrast=c("Group","53WTl6h","Wc153l6h"))
res= results(dds, alpha=alpha,contrast=c("Group","Frq08l6h","53WTl6h"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene="VDAG_JR2_Chr7g03830", intgroup=c("Strain"))


```
========
Time course experiment
========

Using 6hDD, 12hDD, 18hDD, 24hDD
```R
library(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

colData <- colData[!(colData$Sample=="Wc153_DD6_rep1"),]      
countData <- subset(countData, select=-Wc153_DD6_rep1)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep2"),]      
countData <- subset(countData, select=-Wc153_DD6_rep2)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep3"),]      
countData <- subset(countData, select=-Wc153_DD6_rep3)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep1"),]      
countData <- subset(countData, select=-Wc153_LL6_rep1)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep2"),]      
countData <- subset(countData, select=-Wc153_LL6_rep2)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep3"),]      
countData <- subset(countData, select=-Wc153_LL6_rep3)

colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
countData <- subset(countData, select=-Frq08_DD24_rep3)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep1"),]      
countData <- subset(countData, select=-Frq08_LL6_rep1)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep2"),]      
countData <- subset(countData, select=-Frq08_LL6_rep2)
colData <- colData[!(colData$Sample=="Frq08_LL6_rep3"),]      
countData <- subset(countData, select=-Frq08_LL6_rep3)

colData <- colData[!(colData$Sample=="WT53_LL6_rep1"),]
countData <- subset(countData, select=-WT53_LL6_rep1)
colData <- colData[!(colData$Sample=="WT53_LL6_rep2"),]
countData <- subset(countData, select=-WT53_LL6_rep2)
colData <- colData[!(colData$Sample=="WT53_LL6_rep3"),]
countData <- subset(countData, select=-WT53_LL6_rep3)

library(ggplot2)

design <- ~Strain + Time + Strain:Time
ddsTC <- DESeqDataSetFromMatrix(countData,colData,design)
ddsTC <- DESeq(ddsTC, fitType="local")
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

data <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("Time","Strain"), returnData=TRUE)
ggplot(data, aes(x=Time, y=count, color=Strain, group=Strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
dev.off()

#to check a specific gene over time
data <- plotCounts(ddsTC, gene="VDAG_JR2_Chr2g01990",
                   intgroup=c("Time","Strain"), returnData=TRUE)
ggplot(data, aes(x=Time, y=count, color=Strain, group=Strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
dev.off()

#To create a list of the values
write.table(resTC[order(resTC$padj),],"time_course1",sep="\t",quote=F)
```

Using 6hLL, 6hDD, 12hDD, 18hDD, 24hDD
```R
library(DESeq2)
library(ggplot2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

colData <- colData[!(colData$Sample=="Wc153_DD6_rep1"),]      
countData <- subset(countData, select=-Wc153_DD6_rep1)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep2"),]      
countData <- subset(countData, select=-Wc153_DD6_rep2)
colData <- colData[!(colData$Sample=="Wc153_DD6_rep3"),]      
countData <- subset(countData, select=-Wc153_DD6_rep3)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep1"),]      
countData <- subset(countData, select=-Wc153_LL6_rep1)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep2"),]      
countData <- subset(countData, select=-Wc153_LL6_rep2)
colData <- colData[!(colData$Sample=="Wc153_LL6_rep3"),]      
countData <- subset(countData, select=-Wc153_LL6_rep3)

colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
countData <- subset(countData, select=-Frq08_DD24_rep3)


colData$TimePoint <- paste0(colData$Time,colData$Light)
design <- ~Strain + TimePoint + Strain:TimePoint
ddsTC <- DESeqDataSetFromMatrix(countData,colData,design)
ddsTC <- DESeq(ddsTC, fitType="local")
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

data <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("Time","Light","Strain"), returnData=TRUE)
ggplot(data, aes(x=Time, y=count, color=Strain, group=Strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
dev.off()

#to check a specific gene over time
data <- plotCounts(ddsTC, gene="VDAG_JR2_Chr2g01990",
                   intgroup=c("Time","Strain"), returnData=TRUE)
ggplot(data, aes(x=Time, y=count, color=Strain, group=Strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
dev.off()

#To create a list of the values
write.table(resTC[order(resTC$padj),],"time_course2",sep="\t",quote=F)
```

#Functional annotation of JR2 protein files
##Interproscan
Interproscan was used to give gene models functional annotations. Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using 'screen' to allow jobs to be submitted and monitored in the background. This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs was redirected to a temporary output file named interproscan_submission.log .

```bash
screen -a
  cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Genes in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.pep.all.fa); do
  echo $Genes
  $ProgDir/sub_interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
  ```

Following interproscan annotation split files were combined using the following commands:

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
 for Proteins in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.pep.all.fa); do
   Strain=$(echo $Proteins | rev | cut -d '/' -f2 | rev)
   echo $Strain
   InterProRaw=gene_pred/interproscan/$Strain/raw: url "title"
   $ProgDir/append_interpro.sh $Proteins $InterProRaw
 done
 ```


 #===============================================================================
 #       FPKM
 #===============================================================================
```R
 mygenes <- read.table ("genes.txt",header=T,sep="\t")
countData <- read.table("countData2",header=T,sep="\t")
 library(naturalsort)
 t1 <- counts(dds)
 t1 <- mygenes[geneid(t1)]
 rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)
 myfpkm <- fpkm(dds,robust=T)
 myfpkm <- cbind(myfpkm, t1@ranges@width)
 rm(t1)
 myfpkm <- myfpkm[naturalsort(rownames(myfpkm)),]
