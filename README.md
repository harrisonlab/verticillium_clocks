# Verticillium_clocks
==========

Documentation of identification of clock genes in verticillium genomes


Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/verticclium/clocks

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:


#Building of directory structure

#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash

  # For original sequencing runs
  mkdir -p /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
	cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
	Species="V.dahliae"
	Strain="51"
	mkdir -p raw_dna/paired/$Species/$Strain/F
	mkdir -p raw_dna/paired/$Species/$Strain/R    
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/51/F/wilt_51_F_appended.fastq raw_dna/paired/$Species/$Strain/F/.
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/51/R/wilt_51_R_appended.fastq raw_dna/paired/$Species/$Strain/R/.
  Strain="53"
	mkdir -p raw_dna/paired/$Species/$Strain/F
	mkdir -p raw_dna/paired/$Species/$Strain/R    
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/53/F/wilt_53_F_appended.fastq raw_dna/paired/$Species/$Strain/F/.
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/53/R/wilt_53_R_appended.fastq raw_dna/paired/$Species/$Strain/R/.
  Strain="58"
	mkdir -p raw_dna/paired/$Species/$Strain/F
	mkdir -p raw_dna/paired/$Species/$Strain/R    
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/58/F/wilt_58_F_appended.fastq raw_dna/paired/$Species/$Strain/F/.
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/58/R/wilt_58_R_appended.fastq raw_dna/paired/$Species/$Strain/R/.
  Strain="61"
	mkdir -p raw_dna/paired/$Species/$Strain/F
	mkdir -p raw_dna/paired/$Species/$Strain/R    
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/61/F/wilt_61_F_appended.fastq raw_dna/paired/$Species/$Strain/F/.
  cp /home/groups/harrisonlab/project_files/verticillium_dahliae/wilt/raw_dna/paired/V.dahliae/61/R/wilt_61_R_appended.fastq raw_dna/paired/$Species/$Strain/R/.
  # For new sequencing run
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160412_M04465_0010-AMLCU    
  Species="V.dahliae"
  Strain="12008"
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDat/Vd12008_S1_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/.
  cp $RawDat/Vd12008_S1_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/.
```



This process was repeated for RNAseq data:

```bash

```




#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

<!-- Firstly, those strains with more than one run were identified:

```bash
	for Strain in $(ls -d raw_dna/paired/*/*); do
	NumReads=$(ls $Strain/F/*.gz | wc -l);
		if [ $NumReads -gt 1 ]; then
			echo "$Strain";
			echo "$NumReads";
		fi;
	done
```

```
	raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	2
	raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	2
``` -->

Trimming was first performed on all strains that had a single run of data:

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```
<!--
Trimming was then performed for strains with multiple runs of data

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
	echo "Fus2"
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	ReadsF=$(ls $StrainPath/F/s_6_1_sequence.fastq.gz)
	ReadsR=$(ls $StrainPath/R/s_6_2_sequence.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	ReadsF=$(ls $StrainPath/F/FUS2_S2_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/FUS2_S2_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	echo "HB6"
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	ReadsF=$(ls $StrainPath/F/HB6_S4_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/HB6_S4_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	ReadsF=$(ls $StrainPath/F/HB6_S5_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/HB6_S5_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
``` -->


Data quality was visualised once again following trimming:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz ); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
	for TrimPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```
<!--
and for strains with muiltiple runs of data:

```bash
	for TrimPath in $(ls -d raw_dna/paired/*/Fus2); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/s_6_1_sequence.fastq.gz)
		TrimR1=$(ls $TrimPath/R/s_6_2_sequence.fastq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/FUS2_S2_L001_R1_001.fastq.gz)
		TrimR2=$(ls $TrimPath/R/FUS2_S2_L001_R2_001.fastq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
	for TrimPath in $(ls -d raw_dna/paired/*/HB6); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/HB6_S4_L001_R1_001.fastq.gz)
		TrimR1=$(ls $TrimPath/R/HB6_S4_L001_R2_001.fastq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/HB6_S5_L001_R1_001.fastq.gz)
		TrimR2=$(ls $TrimPath/R/HB6_S5_L001_R2_001.fastq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
``` -->

mode kmer abundance prior to error correction was reported using the following
commands:

```bash
for File in $(ls qc_dna/kmc/*/*/*_true_kmer_summary.txt); do
basename $File;
tail -n3 $File | head -n1 ;
done
```

Results were as follows. Where notes next to the kmer coverage reports a fail,
the results were initally <5 but manual inspection of the kmer distribution graph
revealed the coverage to be greater, an that the <5 result was a result of poor
thresholding and did not represeent median coverage.

```
	PG8_true_kmer_summary.txt
	The mode kmer abundance is:  53 < fail
	125_true_kmer_summary.txt
	The mode kmer abundance is:  46 < pass
	55_true_kmer_summary.txt
	The mode kmer abundance is:  29 < pass
	A1-2_true_kmer_summary.txt
	The mode kmer abundance is:  16 < pass
	A13_true_kmer_summary.txt
	The mode kmer abundance is:  22 < pass
	A23_true_kmer_summary.txt
	The mode kmer abundance is:  33 < pass
	A28_true_kmer_summary.txt
	The mode kmer abundance is:  36 < pass
	CB3_true_kmer_summary.txt
	The mode kmer abundance is:  21 < pass
	D2_true_kmer_summary.txt
	The mode kmer abundance is:  11 < pass
	Fus2_true_kmer_summary.txt
	The mode kmer abundance is:  109 < 2lib
	HB17_true_kmer_summary.txt
	The mode kmer abundance is:  27 < pass
	HB6_true_kmer_summary.txt
	The mode kmer abundance is:  91 < 2 lib
	PG_true_kmer_summary.txt
	The mode kmer abundance is:  58 < pass
	N139_true_kmer_summary.txt
	The mode kmer abundance is:  26 < pass
	FOP1_true_kmer_summary.txt
	The mode kmer abundance is:  32 <- fail
	FOP5_true_kmer_summary.txt
	The mode kmer abundance is:  62 <- fail
	L5_true_kmer_summary.txt
	The mode kmer abundance is:  35 <- pass
	PG18_true_kmer_summary.txt
	The mode kmer abundance is:  24 <- pass
	PG3_true_kmer_summary.txt
	The mode kmer abundance is:  28 <- pass
	A8_true_kmer_summary.txt
	The mode kmer abundance is:  21 <- pass
```

#Assembly

Assembly was performed with:
* Spades

## Spades Assembly



```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
			sleep 5m
			printf "."
			Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
		done		
		printf "\n"
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
	done
```

Assembly for PG8, FOP1 and FOP5 failed due to a lack of memory, as such the assembly was
resubmitted with more RAM.

```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'FOP5' -e 'PG8' -e 'FOP1'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 10
	done
```
<!--
Assemblies were submitted for genomes with data from multiple sequencing runs:

```bash
for StrainPath in $(ls -d qc_dna/paired/F.*/HB6); do
  echo $StrainPath
    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/HB6_S4_L001_R1_001_trim.fq.gz);
    TrimR1_Read=$(ls $StrainPath/R/HB6_S4_L001_R2_001_trim.fq.gz);
    TrimF2_Read=$(ls $StrainPath/F/HB6_S5_L001_R1_001_trim.fq.gz);
    TrimR2_Read=$(ls $StrainPath/R/HB6_S5_L001_R2_001_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
  done
	for StrainPath in $(ls -d qc_dna/paired/F.*/Fus2); do
	  echo $StrainPath
	    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
	    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
	    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
	    echo $Strain
	    echo $Organism
	    TrimF1_Read=$(ls $StrainPath/F/s_6_1_sequence_trim.fq.gz);
	    TrimR1_Read=$(ls $StrainPath/R/s_6_2_sequence_trim.fq.gz);
	    TrimF2_Read=$(ls $StrainPath/F/FUS2_S2_L001_R1_001_trim.fq.gz);
	    TrimR2_Read=$(ls $StrainPath/R/FUS2_S2_L001_R2_001_trim.fq.gz);
	    echo $TrimF1_Read
	    echo $TrimR1_Read
	    echo $TrimF2_Read
	    echo $TrimR2_Read
	    OutDir=assembly/spades/$Organism/$Strain
	    qsub $ProgDir/subSpades_2lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
	  done
``` -->

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep 'FOP2'); do
	    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
