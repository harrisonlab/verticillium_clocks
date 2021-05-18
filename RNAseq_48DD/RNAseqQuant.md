# RNA-Seq analysis 

#### Reorganise raw data

```bash
mkdir -p /projects/vertclock2020
cd /projects/vertclock2020

mkdir -p RNAseq_data/48DD_experiment2020/WT53/T04/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T04/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T08/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T08/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T12/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T12/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T16/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T16/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T20/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T20/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T24/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T24/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T28/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T28/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T32/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T32/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T36/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T36/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T40/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T40/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T44/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T44/R/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T48/F/
mkdir -p RNAseq_data/48DD_experiment2020/WT53/T48/R/

cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_4C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T04/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_8C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T08/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_12C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T12/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_16C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T16/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_20C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T20/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_24C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T24/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_28C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T28/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_32C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T32/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_36C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T36/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_40C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T40/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_44C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T44/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_48A_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_48A_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_48B_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_48B_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/R/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/rawdata/WT53_48C_1.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/F/
cp /archives/2020_niab_general/20201117_VdahliaeClock_RNAseq48DD/upload1/WT53_48C_2.fq.gz RNAseq_data/48DD_experiment2020/WT53/T48/R/
```

## Perform qc on RNA-Seq data

```bash
# Run fastqc
    for RawData in $(ls RNAseq_data/48DD_experiment2020/WT53/*/*/*.fq.gz | grep 'WT53_1\|WT53_2\|WT53_3\|WT53_4\|WT53_8') 
    do
        echo $RawData
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch $ProgDir/fastqc.sh $RawData
    done
```


```bash
# Run fastq-mcf
    for RNADir in $(ls -d RNAseq_data/48DD_experiment2020/WT53/T08); do
    FileNum=$(ls $RNADir/F/*_1.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
            sbatch $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA
        done
    done
```

```bash
    # Run fastqc
    for RawData in $(ls qc_rna/48DD_experiment2020/WT53/T08/*/*.fq.gz | grep 'WT53_1\|WT53_2\|WT53_3\|WT53_4\|WT53_8') 
    do
        echo $RawData
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch $ProgDir/fastqc.sh $RawData
    done
```

## Decontamination of rRNA reads in RNAseq data

```bash
    for RNADir in $(ls -d qc_rna/48DD_experiment2020/WT53/*); do
    FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            echo $RNADir
            Strain=$(sed 's/.*\///' <<< $RNADir)
            Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            echo $Sample_Name
            sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned/$Timepoint/$Sample_Name $FileF $FileR $ProgDir $Strain
        done
    done
```

## Salmon 

```bash
    # I have installed the latest version in a new env
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda create -n salmon salmon
```

```bash
    for Transcriptome in $(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
        Organism=V.dahliae
        echo "$Organism - $Strain"
        for RNADir in $(ls -d qc_rna/48DD_experiment2020/WT53/*/cleaned/*/*); do
            FileNum=$(ls $RNADir/F/*_1_cleaned.fq.gz | wc -l)
            for num in $(seq 1 $FileNum); do
                printf "\n"
                FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
                FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
                echo $FileF
                echo $FileR
                Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
                Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
                Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
                echo "$Prefix"
                echo "$Timepoint"
                echo "$Sample_Name"
                OutDir=alignment/salmon/48DD_experiment/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
                ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
                sbatch -p himem $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
            done
        done
    done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/48DD_experiment/V.dahliae/JR2/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/48DD_experiment/V.dahliae/JR2/*/*/*/quant.sf); do
  Prefix=$(echo $File | cut -f8 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/$Prefix
  cp $PWD/$File RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/vertclock2020")

#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(BiocParallel)
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
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggrepel)

library(MetaCycle)
library(docker4seq)



#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

write.table(txi.genes,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txireps.txt",sep="\t",na="",quote=F)


pca(experiment.table="txigenes2.txt", type="TPM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="txigenes_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/Vd_clocks_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Define the DESeq 'GLM' model
design <- ~ Condition
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
###
[1] "Intercept"                    "Condition_WT53_16_vs_WT53_12"
 [3] "Condition_WT53_20_vs_WT53_12" "Condition_WT53_24_vs_WT53_12"
 [5] "Condition_WT53_28_vs_WT53_12" "Condition_WT53_32_vs_WT53_12"
 [7] "Condition_WT53_36_vs_WT53_12" "Condition_WT53_4_vs_WT53_12" 
 [9] "Condition_WT53_40_vs_WT53_12" "Condition_WT53_44_vs_WT53_12"
[11] "Condition_WT53_48_vs_WT53_12" "Condition_WT53_8_vs_WT53_12" 
###

#Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


write.csv(vst, file="vst.csv")

# Correcting the experimental desing

# Define the DESeq 'GLM' model
design <- ~ Condition + Experiment
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Library normalisation (this is done in the DeSeq function)
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
resultsNames(dds)

> resultsNames(dds)
 [1] "Intercept"                    "Condition_WT53_16_vs_WT53_12"
 [3] "Condition_WT53_20_vs_WT53_12" "Condition_WT53_24_vs_WT53_12"
 [5] "Condition_WT53_28_vs_WT53_12" "Condition_WT53_32_vs_WT53_12"
 [7] "Condition_WT53_36_vs_WT53_12" "Condition_WT53_4_vs_WT53_12" 
 [9] "Condition_WT53_40_vs_WT53_12" "Condition_WT53_44_vs_WT53_12"
[11] "Condition_WT53_48_vs_WT53_12" "Condition_WT53_8_vs_WT53_12" 
[13] "Experiment_B_vs_A"            "Experiment_C_vs_A" 

# Varianze stabilizing transformation

vst<-varianceStabilizingTransformation(dds)
write.csv(assay(vst), file="vst_true.csv")

pdf("PCA_vst_true.pdf")
plotPCA(vst,intgroup=c("Experiment"))
dev.off()

# Remove batch effect associated with the Experiment (replicas)
mat1 <- assay(vst)
mat1 <- limma::removeBatchEffect(mat1, vst$Experiment)
write.csv(mat1, file="vst_true_batchcorrected.csv")
assay(vst) <- mat1

pdf("PCA_vst_true_corrected.pdf")
plotPCA(vst,intgroup=c("Experiment"))
dev.off()

# Blind false
vst4<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst4), file="vst_false.csv")

pdf("PCA_vst_false.pdf")
plotPCA(vst4,intgroup=c("Experiment"))
dev.off()

# Remove batch effect associated with the Experiment (replicas)
mat <- assay(vst4)
mat <- limma::removeBatchEffect(mat, vst4$Experiment)
write.csv(mat, file="vst_false_batchcorrected.csv")
assay(vst4) <- mat

pdf("PCA_vst_false_corrected.pdf")
plotPCA(vst4,intgroup=c("Experiment"))
dev.off()

# Test vst function 
vst2<-vst(dds,blind=FALSE)
write.csv(assay(vst2), file="vst2.csv")
vst3<-vst(dds,blind=TRUE)
write.csv(assay(vst3), file="vst3.csv")

pdf("PCA_vst2.pdf")
plotPCA(vst2,intgroup=c("Experiment"))
dev.off()
pdf("PCA_vst3.pdf")
plotPCA(vst3,intgroup=c("Experiment"))
dev.off()

C2 <- assay(vst2)
C2 <- limma::removeBatchEffect(C2, vst2$Experiment)
write.csv(C2, file="vst2_batchcorrected.csv")
assay(vst2) <- C2
pdf("PCA_vst2_corrected.pdf")
plotPCA(vst2,intgroup=c("Experiment"))
dev.off()

C3 <- assay(vst3)
C3 <- limma::removeBatchEffect(C3, vst3$Experiment)
write.csv(C3, file="vst3_batchcorrected.csv")
assay(vst3) <- C3
pdf("PCA_vst3_corrected.pdf")
plotPCA(vst3,intgroup=c("Experiment"))
dev.off()


#Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"C+E/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"C+E/normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"C+E/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"C+E/fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="C+E/raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="C+E/normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="C+E/fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="C+E/fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

# Full PCA with sample names
data <- plotPCA(vst4, intgroup=c("Condition","Experiment"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Experiment)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(vst4))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
ggsave("PCA_sample_names1.pdf", pca_plot, dpi=300, height=10, width=12)


setwd("/Users/antoniogomez/Desktop/VertRNA")

require(DESeq2)
colData <- read.table("colData",header=T,sep="\t")
countData <- read.table("countData_WT",header=T,sep="\t")
colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)

colData <- read.table("/Users/antoniogomez/Desktop/VertRNA/Vd_clocks_RNAseq_design.txt",header=T,sep="\t")
countData <- read.table("/Users/antoniogomez/Desktop/VertRNA/C+E/fpkm_norm_counts.txt",header=T,sep="\t")

write.table(rownames(countData,"countData_annot"))

#Do it in Bash
cut -f1 countData_WT > counData_annot

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("counData_annot.txt")
data <- read.delim("TPM/TPM_id2.txt")


jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/VertRNA/TPM/JTK_WT53_2.txt"),row.names=F,col.names=T,quote=F,sep="\t")


source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("counData_annot.txt")
data <- read.delim("fpkm/fpkm4jtk_counts.txt")


jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/VertRNA/fpkm/JTK_WT53_fpkm_counts.txt"),row.names=F,col.names=T,quote=F,sep="\t")

```

Interproscan


```bash
awk -F'.p1' '{print $1}' public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.pep.all.fa > genes.pep.fa
```


```bash
# This command will split your gene fasta file and run multiple interproscan jobs.
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Genes in $(ls genes.pep.fa); do
    echo $Genes
    $ProgDir/interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```
```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Proteins in $(ls genes.pep.fa); do
    Strain=JR2
    Organism=V.dahliae
    InterProRaw=gene_pred/interproscan/V.dahliae/JR2/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done


for GeneGff in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3); do
Strain=JR2
Organism=V.dahliae
Assembly=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InterPro=$(ls /projects/oldhome/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
#SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
OutDir=gene_pred/annotation/$Organism/$Strain
mkdir -p $OutDir
GeneFasta=$(ls genes.pep.fa)
ProgDir=/home/gomeza/git_repos/scripts/verticillium_clocks/annotation
$ProgDir/Vd_annotation_tables_vAG.py --gene_gff $GeneGff --gene_fasta $GeneFasta --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp.tsv
done

cat $OutDir/"$Strain"_gene_table_incl_exp.tsv | cut -f1 | awk -F'.t1' '{print $1}' > table.tsv

```
```r
source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("counData_annot.txt")
data <- read.delim("fpkm/normalised_counts4jtk.txt")


jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

periods <- 12       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/VertRNA/fpkm/JTK_WT53_norm_counts.txt"),row.names=F,col.names=T,quote=F,sep="\t")

same results using norm counts or fpkm robus yes

source("JTKversion3/JTK_CYCLEv3.1.R")
project <- "WT53"
options(stringsAsFactors=FALSE)
annot <- read.delim("counData_annot.txt")
data <- read.delim("fpkm/fpkm_robust/fpkm4jtk.txt")


jtkdist(12, 3)       # 13 total time points, 2 replicates per time point

periods <- 2:11       # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,4)  # 4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("/Users/antoniogomez/Desktop/VertRNA/fpkm/fpkm_robust/last_fpkm_norm2.txt"),row.names=F,col.names=T,quote=F,sep="\t")
```
T1<-read.delim("inter.txt",header=T)
T2<-read.table("fpkm/fpkm_robust/01.txt",header=T,sep="\t")
T3<-merge(T2,T1, by.x="ID",by.y="ID",all.x=TRUE)
write.table(T3,"TableAG.txt",sep="\t",na="",quote=F)
T3<-read.table("Table3.txt",header=T,sep="\t")

write.table(T3,"Echointer.txt",sep="\t",na="",quote=F)
write.table(T5,"jtkinter.txt",sep="\t",na="",quote=F)
