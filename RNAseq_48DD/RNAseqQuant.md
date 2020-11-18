# RNA-Seq analysis 

```bash
mkdir -p /home/scratch/gomeza/rna_seq/
tar -xvf /home/scratch/gomeza/rna_seq/20171211/C101HW17030405_20180102_5_Yvad6z.tar
```

##Reorganise raw data

```bash
mkdir -p N.ditissima/Hg199/mycellium/F/
mkdir -p N.ditissima/Hg199/mycellium/R/
mkdir -p N.ditissima/GD/t0/F/
mkdir -p N.ditissima/GD/t1/F/
mkdir -p N.ditissima/GD/t2/F/
mkdir -p N.ditissima/GD/t0/R/
mkdir -p N.ditissima/GD/t1/R/
mkdir -p N.ditissima/GD/t2/R/
mkdir -p N.ditissima/M9/t0/F/
mkdir -p N.ditissima/M9/t1/F/
mkdir -p N.ditissima/M9/t2/F/
mkdir -p N.ditissima/M9/t0/R/
mkdir -p N.ditissima/M9/t1/R/
mkdir -p N.ditissima/M9/t2/R/

mv 20171211/C101HW17030405/raw_data/GD_C1_3_1.fq.gz N.ditissima/GD/t0/F/
mv 20171211/C101HW17030405/raw_data/GD_C1_3_2.fq.gz N.ditissima/GD/t0/R/
mv 20171211/C101HW17030405/raw_data/GD_6A3_1.fq.gz N.ditissima/GD/t1/F/
mv 20171211/C101HW17030405/raw_data/GD_6A3_2.fq.gz N.ditissima/GD/t1/R/
mv 20171211/C101HW17030405/raw_data/M9_C1_3_2.fq.gz N.ditissima/M9/t0/R/
mv 20171211/C101HW17030405/raw_data/M9_C1_3_1.fq.gz N.ditissima/M9/t0/F/
mv 20171211/C101HW17030405/raw_data/M9_6A3_1.fq.gz N.ditissima/M9/t1/F/
mv 20171211/C101HW17030405/raw_data/M9_6A3_2.fq.gz N.ditissima/M9/t1/R/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A4_1.fq.gz N.ditissima/GD/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A4_2.fq.gz N.ditissima/GD/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A3_1.fq.gz N.ditissima/GD/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A3_2.fq.gz N.ditissima/GD/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A5_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A5_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A2_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_6A2_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_6A2_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A2_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C2_3_1.fq.gz N.ditissima/GD/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C2_3_2.fq.gz N.ditissima/GD/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C3_3_2.fq.gz N.ditissima/GD/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C3_3_1.fq.gz N.ditissima/GD/t0/F/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C2_3_1.fq.gz N.ditissima/M9/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C3_3_1.fq.gz N.ditissima/M9/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C3_3_2.fq.gz N.ditissima/M9/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C2_3_2.fq.gz N.ditissima/M9/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A3_1.fq.gz N.ditissima/M9/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A3_2.fq.gz N.ditissima/M9/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_6A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_6A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A4_1.fq.gz N.ditissima/M9/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A4_2.fq.gz N.ditissima/M9/t1/R/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_1_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_2_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_3_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_3_2.fq.gz N.ditissima/Hg199/mycelium/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_2_2.fq.gz N.ditissima/Hg199/mycelium/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_1_2.fq.gz N.ditissima/Hg199/mycelium/R/
```

## Perform qc on RNA-Seq timecourse and mycelium data
```bash
# Run fastqc
for Strain in Strain1 Strain2; do
    RawData=$(ls raw_dna/paired/$Organism/$Strain/*/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
done
```


```bash
# Run fastq-mcf
for Strain in Strain1 Strain2; do
    Read_F=raw_dna/paired/*/*/F/*.fastq.gz
    Read_R=raw_dna/paired/*/*/R/*.fastq.gz
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
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
    for Transcriptome in $(ls path/to/predicted/transcriptome/final_genes_appended_renamed.cdna.fasta); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
        echo "$Organism - $Strain"
        for RNADir in $(ls -d path/to/RNAseq/reads); do
        FileF=$(ls $RNADir/*.1.fq) # grep -e 'Sample name'
        FileR=$(ls $RNADir/*.2.fq) # grep -e 'Sample name'
        echo $FileF
        echo $FileR
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/.1.fq//g')
        echo $Sample_Name
        OutDir=RNAseq_analysis/salmon/$Organism/$Strain/$Sample_Name
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
        sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done
```
