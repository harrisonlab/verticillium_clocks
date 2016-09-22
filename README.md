# Verticillium_clocks
==========

Documentation of identification of clock genes in verticillium genomes


Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks

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

To gzip files:

for File in $(ls raw_dna/paired/*/*/*/*.fastq); do
gzip $File > $File.gz
done

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
		ProgDir=/home/lopeze/git_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

<!-- Firstly, those strains with more than one run were identified:

```bash
	ls

```

```
	raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	2
	raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	2
``` -->

Trimming was first performed on all strains that had a single run of data:

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/*); do
		ProgDir=/home/lopeze/git_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/lopeze/git_repos/tools/seq_tools/ncbi_adapters.fa
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
		ProgDir=/home/lopeze/git_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data
```bash
		for TrimPath in $(ls -d qc_dna/paired/*/*); do
		ProgDir=/home/lopeze/git_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fq.gz)
		TrimR=$(ls $TrimPath/R/*.fq.gz)
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
51_true_kmer_summary.txt
The mode kmer abundance is:  52
53_true_kmer_summary.txt
The mode kmer abundance is:  37
58_true_kmer_summary.txt
The mode kmer abundance is:  71
61_true_kmer_summary.txt
The mode kmer abundance is:  53
```

#Assembly

Assembly was performed with:
* Spades

## Spades Assembly



```bash
    for StrainPath in $(ls -d qc_dna/paired/*/* ); do
    ProgDir=/home/lopeze/git_repos/tools/seq_tools/assemblers/spades
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

Quast --> quality assessment tool for genome assemblies.

```bash
ProgDir=/home/lopeze/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
51 reults:
Assembly                   contigs_min_500bp  contigs_min_500bp broken      
# contigs                  1482               1557                    
Largest contig             321888             321888                  
Total length               33065465           33064715                
GC (%)                     55.93              55.93                   
N50                        64178              59021                   
N75                        34220              32519                   

53 results:                
Assembly                   contigs_min_500bp  contigs_min_500bp broken         
# contigs                  1387               1479                    
Largest contig             237640             193510                  
Total length               32389379           32388393                
GC (%)                     56.40              56.40                   
N50                        47913              45145                   
N75                        26178              24006                   

58 results:
Assembly                   contigs_min_500bp  contigs_min_500bp broken             
# contigs                  1160               1262                    
Largest contig             294690             294690                  
Total length               32601233           32600213                
GC (%)                     55.87              55.87                   
N50                        79211              68028                   
N75                        41540              37617                   


61 results:
Assembly                   contigs_min_500bp  contigs_min_500bp broken                 
# contigs                  1239               1364               
Largest contig             359375             331650                  
Total length               32897558           32896308                
GC (%)                     55.64              55.64                   
N50                        81559              71457                   
N75                        41700              36159                   



#Repeatmasking

Repeat masking was performed with:
* Repeatmasker
* Repeatmodeler

The best assembly was used to perform Repeatmasking

```bash

for BestAssembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta);
do
ProgDir=/home/lopeze/git_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $BestAssembly
qsub $ProgDir/transposonPSI.sh $BestAssembly
done
```

Strain 51:
** % bases masked by hardmasked and softmasked: 3.45% (bases masked:1142147 bp)

** % bases masked by transposon psi: 1.88% (bases masked:622307 bp)

Strain 53:
** % bases masked by hardmasked and softmasked: 2.07% (bases masked:684231 bp)

** % bases masked by transposon psi: 0.58% (bases masked:189342 bp)

Strain 58:
** % bases masked by hardmasked and softmasked: 4% (bases masked:1302941 bp)

** % bases masked by transposon psi: 2.45% (bases masked:798225 bp)

Strain 61:
** % bases masked by hardmasked and softmasked: 4.78% (bases masked:1573657 bp)

** % bases masked by transposon psi: 3.24 % (bases masked:1065488 bp)


Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "$OutFile"
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
```
repeat_masked/V.dahliae/51/filtered_contigs_repmask/51_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1319234
repeat_masked/V.dahliae/53/filtered_contigs_repmask/53_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
857005
repeat_masked/V.dahliae/58/filtered_contigs_repmask/58_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1456508
repeat_masked/V.dahliae/61/filtered_contigs_repmask/61_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1726406



cp /home/armita/generic_profiles/2016-07-28/.generic_profile /home/lopeze/.profile
. ~/.profile
