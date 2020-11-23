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
for RNADir in $(ls -d RNAseq_data/48DD_experiment2020/WT53/T48); do
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

for RawData in $(ls RNAseq_data/48DD_experiment2020/WT53/*/*/*.fq.gz | grep 'WT53_1\|WT53_2\|WT53_3\|WT53_4\|WT53_8') 
do
echo $RawData
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
sbatch $ProgDir/fastqc.sh $RawData
done
```
```

```bash
# Run fastqc
for Strain in Strain1 Strain2; do
    RawData=$(ls raw_dna/paired/$Organism/$Strain/*/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
done
```


## Salmon 

```bash
for Transcriptome in $(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
    Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
    Organism=V.dahliae
    echo "$Organism - $Strain"
    for RNADir in $(ls -d RNAseq_data/48DD_experiment2020/WT53/*); do
    FileNum=$(ls $RNADir/F/*_1.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
            Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1.fq.gz//g')
            echo "$Prefix"
            echo "$Timepoint"
            echo "$Sample_Name"
            OutDir=alignment/salmon/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done
done
```
