A Bioproject and Biosample was made with NCBI genbank for submission of genomes. Following the creation of these submissions, the .fasta assembly was uploaded through the submission portal. A note was provided requesting that the assembly be run through the contamination screen to aid a more detailed resubmission in future. The returned FCSreport.txt was downloaded from the NCBI webportal and used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
for Assembly in $(ls assembly/spades/V.dahliae/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
  done
  ```

  These downloaded files were used to correct assemblies:

```bash
  for Assembly in $(ls assembly/spades/V.dahliae/*/filtered_contigs/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/FCSreport.txt)
  OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
  mkdir -p $OutDir
  ProgDir=/home/lopeze/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp.fasta --coord_file $NCBI_report > $OutDir/log.txt
  done
```

Strain 62 second round of correction

```bash
  for Assembly in $(ls assembly/spades/V.dahliae/61/ncbi_edits/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  NCBI_report=$(ls genome_submission/$Organism/$Strain/second_submission/FCSreport.txt)
  OutDir=assembly/spades/$Organism/$Strain/ncbi_edits_2
  mkdir -p $OutDir
  ProgDir=/home/lopeze/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp.fasta --coord_file $NCBI_report > $OutDir/log.txt
  done
```

Quast was used to collect details on these assemblies again

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Quast for 61 second round


```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/61/ncbi_edits_2/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits_2
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

12251
AAssembly                   contigs_min_500bp  contigs_min_500bp broken
contigs (>= 0 bp)        1483               1555                    
contigs (>= 1000 bp)     1096               1163                    
Total length (>= 0 bp)     33028358           33027638                
Total length (>= 1000 bp)  32766912           32763294                
contigs                  1483               1555                    
Largest contig             321888             321888                  
Total length               33028358           33027638                
GC (%)                     55.96              55.96                   
N50                        63849              59021                   
N75                        34220              32528                   
L50                        154                166                     
L75                        326                350                     
N's per 100 kbp          2.32               0.14    

12158
Assembly                   contigs_min_500bp  contigs_min_500bp broken
contigs (>= 0 bp)        1155               1257
contigs (>= 1000 bp)     898                997
Total length (>= 0 bp)     32560838           32559818
Total length (>= 1000 bp)  32385589           32382197
contigs                  1154               1256
Largest contig             294690             294690
Total length               32560355           32559335
GC (%)                     55.89              55.89
N50                        79211              68028
N75                        41647              37762
L50                        127                142
L75                        266                304
N's per 100 kbp          3.25               0.11

12253
Assembly                   contigs_min_500bp  contigs_min_500bp broken
contigs (>= 0 bp)        1382               1473
contigs (>= 1000 bp)     1263               1347
Total length (>= 0 bp)     32355537           32354561
Total length (>= 1000 bp)  32270360           32264719
contigs                  1382               1473
Largest contig             237640             193510
Total length               32355537           32354561
GC (%)                     56.42              56.42
N50                        47913              45158
N75                        26456              24013
L50                        206                218
L75                        427                458
 N's per 100 kbp          3.08               0.06

 Assembly                   contigs_min_500bp  contigs_min_500bp broken
contigs (>= 0 bp)        1238               1363                    
contigs (>= 1000 bp)     964                1085                    
Total length (>= 0 bp)     32862437           32861187                
Total length (>= 1000 bp)  32674024           32669801                
contigs                  1237               1362                    
Largest contig             359375             331650                  
Total length               32861946           32860696                
GC (%)                     55.66              55.66                   
N50                        81782              71457                   
N75                        41683              36159                   
L50                        126                146                     
L75                        266                314                     
N's per 100 kbp          3.91               0.11       


##Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/lopeze/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp.fasta | grep -v '61'); do
Organism=$(echo $BestAss | rev | cut -d "/" -f4 | rev)
Strain=$(echo $BestAss | rev | cut -d "/" -f3 | rev)
OutDir=repeat_masked/$Organism/$Strain/ncbi_filtered_contigs_repmask
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```

A differemt loaction was used for isolate 61 as this required multiple rounds of ncbi edits

```bash
ProgDir=/home/lopeze/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/61/ncbi_edits_2/contigs_min_500bp.fasta); do
Organism=$(echo $BestAss | rev | cut -d "/" -f4 | rev)
Strain=$(echo $BestAss | rev | cut -d "/" -f3 | rev)
OutDir=repeat_masked/$Organism/$Strain/ncbi_filtered_contigs_repmask
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
  ```  
The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/V.*/*/ncbi_filtered_contigs_repmask | grep '61'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
  ```
  V.dahliae       51
The number of bases masked by RepeatMasker:     1195031
The number of bases masked by TransposonPSI:    310921
The total number of masked bases are:   1377060

  V.dahliae       58
  The number of bases masked by RepeatMasker:     1258760
  The number of bases masked by TransposonPSI:    360475
  The total number of masked bases are:   1418679

  V.dahliae       53
  The number of bases masked by RepeatMasker:     689605
  The number of bases masked by TransposonPSI:    221954
  The total number of masked bases are:   863683

  V.dahliae       61
  The number of bases masked by RepeatMasker:	1574548
  The number of bases masked by TransposonPSI:	407135
  The total number of masked bases are:	1726438

  Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

  ```bash
  for File in $(ls repeat_masked/V.*/*/ncbi*/*_contigs_softmasked.fa); do
  OutDir=$(dirname $File)
  TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
  OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
  bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
  echo "$OutFile"
  echo "Number of masked bases:"
  cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
  done
  ```
  12151
Number of masked bases:
1377060
12253
Number of masked bases:
863683
12158
Number of masked bases:
1418679
12161
Number of masked bases:
1726438

#Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

  ```bash
ProgDir=/home/lopeze/git_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
for Genome in $(ls repeat_masked/V.*/*/ncbi*/*_contigs_softmasked.fa); do
echo $Genome;
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
  ```
12151
*** Number of cegma genes present and complete: 94.76% ** Number of cegma genes present and partial: 97.98%
12253
*** Number of cegma genes present and complete: 94.35% ** Number of cegma genes present and partial: 97.98%
12158
*** Number of cegma genes present and complete: 95.16% ** Number of cegma genes present and partial: 97.18%
12161
*** Number of cegma genes present and complete: 95.56% ** Number of cegma genes present and partial: 97.58%


  ```bash
ProgDir=/home/lopeze/git_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
for Genome in $(ls repeat_masked/V.*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
echo $Genome;
 qsub $ProgDir/sub_cegma.sh $Genome dna;
done
    ```
12151
*** Number of cegma genes present and complete: 94.76% ** Number of cegma genes present and partial: 97.98%
12253
*** Number of cegma genes present and complete: 94.35% ** Number of cegma genes present and partial: 97.98%
12158
*** Number of cegma genes present and complete: 95.16% ** Number of cegma genes present and partial: 97.18%
12161
*** Number of cegma genes present and complete: 95.56% ** Number of cegma genes present and partial: 97.58%

Outputs were summarised using the commands:

  ```bash
for File in $(ls gene_pred/cegma/V.*/*/*_dna_cegma.completeness_report); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
Species=$(echo $File | rev | cut -f3 -d '/' | rev);
printf "$Species\t$Strain\n";
cat $File | head -n18 | tail -n+4;printf "\n";
done > gene_pred/cegma/cegma_results_dna_summary.txt
```

#Gene prediction

##Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could be made. To do this tophat and cufflinks were run, aligning the reads against a single genome. The fragment length and stdev were printed to stdout while cufflinks was running.

  ```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/V.*/12008PDA); do
      Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
      echo "$Timepoint"
      FileF=$(ls $RNADir/F/*_trim.fq.gz)
      FileR=$(ls $RNADir/R/*_trim.fq.gz)
      OutDir=ncbi_alignment/$Organism/$Strain/$Timepoint
      ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
    done
  done
    ```
    51 --> 81.7% overall read mapping rate. 73.0% concordant pair alignment rate.
    53 --> 76.5% overall read mapping rate. 68.5% concordant pair alignment rate.
    58 --> 73.7% overall read mapping rate. 65.7% concordant pair alignment rate.
    61 --> 78.8% overall read mapping rate. 70.1% concordant pair alignment rate.

      ```bash
    for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/V.*/12008CD); do
      Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
      echo "$Timepoint"
      FileF=$(ls $RNADir/F/*_trim.fq.gz)
      FileR=$(ls $RNADir/R/*_trim.fq.gz)
      OutDir=ncbi_alignment/$Organism/$Strain/$Timepoint
      ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
    done
  done
    ```
    51 --> 80.4% overall read mapping rate. 70.8% concordant pair alignment rate.
    53 --> 77.4% overall read mapping rate. 68.3% concordant pair alignment rate.
    58 --> 70.2% overall read mapping rate. 61.2% concordant pair alignment rate.
    61 --> 73.2% overall read mapping rate. 63.8% concordant pair alignment rate.


Alignments were concatenated prior to running cufflinks: Cufflinks was run to produce the fragment length and stdev statistics:

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for AcceptedHits in $(ls ncbi_alignment/$Organism/$Strain/*/accepted_hits.bam); do
Timepoint=$(echo $AcceptedHits | rev | cut -f2 -d '/' | rev)
echo $Timepoint
OutDir=gene_pred/ncbi_cufflinks/$Organism/$Strain/"$Timepoint"_prelim
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
done
  ```

  12251
PDA
> Processed 18108 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 11141720.32
>       Raw Map Mass: 11141720.32
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 233.10
>                  Estimated Std Dev: 59.61
[17:57:30] Assembling transcripts and estimating abundances.
> Processed 18245 loci.                        [*************************] 100%

CDD
> Processed 17015 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 8276688.40
>       Raw Map Mass: 8276688.40
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 224.33
>                  Estimated Std Dev: 62.66
[17:43:10] Assembling transcripts and estimating abundances.
> Processed 17155 loci.                        [*************************] 100%

12253
PDA
> Processed 18194 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 10396226.10
>       Raw Map Mass: 10396226.10
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 232.92
>                  Estimated Std Dev: 59.58
[18:04:14] Assembling transcripts and estimating abundances.
> Processed 18333 loci.                        [*************************] 100%   7%

CDD
> Processed 17078 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 7964357.15
>       Raw Map Mass: 7964357.15
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 224.17
>                  Estimated Std Dev: 62.58
[17:43:03] Assembling transcripts and estimating abundances.
> Processed 17223 loci.                        [*************************] 100%

12158
PDA
> Processed 19123 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 10053530.97
>       Raw Map Mass: 10053530.97
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 230.36
>                  Estimated Std Dev: 59.22
[18:01:29] Assembling transcripts and estimating abundances.
> Processed 19284 loci.                        [*************************] 100%

CDD
> Processed 18576 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 7283829.83
>       Raw Map Mass: 7283829.83
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 221.45
>                  Estimated Std Dev: 62.04
[17:50:04] Assembling transcripts and estimating abundances.
> Processed 18745 loci.                        [*************************] 100%

12161
PDA
> Processed 19194 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 10828749.04
>	Raw Map Mass: 10828749.04
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 235.63
>	           Estimated Std Dev: 63.60
[10:50:39] Assembling transcripts and estimating abundances.
> Processed 19343 loci.                        [*************************] 100%

CDD
> Processed 18620 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 7639017.83
>	Raw Map Mass: 7639017.83
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 224.52
>	           Estimated Std Dev: 64.72
[10:36:08] Assembling transcripts and estimating abundances.
> Processed 18786 loci.                        [*************************] 100%


Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/*/12008PDA | grep -v -e '_rep'); do
Timepoint_PDA=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint_PDA"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=ncbi_alignment/$Organism/$Strain/$Timepoint_PDA
InsertGap='-250'
InsertStdDev='64'
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
done
printf "\n"
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
done
done
    ```



    51 --> 81.7% overall read mapping rate. 73.0% concordant pair alignment rate.
    53 --> 76.5% overall read mapping rate. 68.5% concordant pair alignment rate.
    58 --> 73.7% overall read mapping rate. 65.7% concordant pair alignment rate.
    61 --> 78.8% overall read mapping rate. 70.1% concordant pair alignment rate.

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/*/12008CD | grep -v -e '_rep'); do
Timepoint_CD=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint_CD"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=ncbi_alignment/$Organism/$Strain/$Timepoint_CD
InsertGap='-235'
InsertStdDev='64'
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
done
printf "\n"
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
done
done
  ```



  51 --> 80.4% overall read mapping rate. 70.8% concordant pair alignment rate.
  53 --> 77.4% overall read mapping rate. 68.3% concordant pair alignment rate.
  58 --> 70.2% overall read mapping rate. 61.2% concordant pair alignment rate.
  61 --> 73.2% overall read mapping rate. 63.8% concordant pair alignment rate.


#Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

ls ~/.gm_key
cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key

Braker predictiction was performed using softmasked genome, not unmasked one.


```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
done
printf "\n"
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker/$Organism/"$Strain"
AcceptedHits=ncbi_alignment/$Organism/$Strain/concatenated/concatenated.bam
mkdir -p ncbi_alignment/$Organism/$Strain/concatenated
samtools merge -f $AcceptedHits \ncbi_alignment/V.dahliae/$Strain/12008PDA/accepted_hits.bam \ncbi_alignment/V.dahliae/$Strain/12008CD/accepted_hits.bam
GeneModelName="$Organism"_"$Strain"_braker_two
rm -r /home/lopeze/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_two
ProgDir=/home/lopeze/git_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
  ```

##Amino acid sequences and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/V.dahliae/*/V.dahliae_*_braker_two/augustus.gff); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
  ```

The relationship between gene models and aligned reads was investigated. To do this aligned reads needed to be sorted and indexed:

Note - IGV was used to view aligned reads against the Fus2 genome on my local machine.

```bash
for
InBam=$(ls alignment/V.dahliae/*/concatenated/concatenated.bam)
ViewBam=alignment/V.dahliae/*/concatenated/concatenated_view.bam
SortBam=alignment/V.dahliae/*/concatenated/concatenated_sorted
samtools view -b $InBam > $ViewBam
samtools sort $ViewBam $SortBam
samtools index $SortBam.bam
```

#Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and therefore features can not be restricted by strand when they are intersected.


```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/ncbi_cufflinks/$Organism/$Strain/concatenated_prelim
mkdir -p $OutDir
AcceptedHits=ncbi_alignment/$Organism/$Strain/concatenated/concatenated.bam
ProgDir=/home/lopeze/git_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
  ```
  Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary1/$Organism/$Strain
    CufflinksGTF=gene_pred/ncbi_cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
    ProgDir=/home/lopeze/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
    done
    ```

#Identification of duplicated genes in additional CodingQuary gene models

```bash
for AddGenes in $(ls gene_pred/codingquary1/V.*/*/additional/additional_genes.gff); do
Strain=$(echo $AddGenes| rev | cut -d '/' -f3 | rev)
Organism=$(echo $AddGenes | rev | cut -d '/' -f4 | rev)
OutDir=$(dirname $AddGenes)
echo "$Organism - $Strain" > $OutDir/duplicated_genes.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $AddGenes >> $OutDir/duplicated_genes.txt
cat $OutDir/duplicated_genes.txt
echo ""
done
  ```

  V.dahliae - 51
Duplicate gene found:
NODE_279_length_40998_cov_22.1969	24380	24910
NODE_279_length_40998_cov_22.1969	CodingQuarry_v2.0	gene	24380	24910	.	+	.	ID=NS.01835;Name=;
NODE_279_length_40998_cov_22.1969	CodingQuarry_v2.0	gene	24380	24910	.	+	.	ID=CUFF.4517.1.10;Name=;

V.dahliae - 53
Duplicate gene found:
NODE_271_length_40835_cov_18.2922	16060	16590
NODE_271_length_40835_cov_18.2922	CodingQuarry_v2.0	gene	16060	16590	.	-	.	ID=NS.01607;Name=;
NODE_271_length_40835_cov_18.2922	CodingQuarry_v2.0	gene	16060	16590	.	-	.	ID=CUFF.4072.1.5;Name=;

V.dahliae - 58

V.dahliae - 61
Duplicate gene found:
NODE_223_length_50767_cov_26.2663	12238	12558
NODE_223_length_50767_cov_26.2663	CodingQuarry_v2.0	gene	12238	12558	.	+	.	ID=CUFF.3883.1.2;Name=;
NODE_223_length_50767_cov_26.2663	CodingQuarry_v2.0	gene	12238	12558	.	+	.	ID=CUFF.3884.1.3;Name=;


To check how many proteins/genes were predicted by coding quary for each strain:
gene_pred/codingquary1/V.dahliae/51/out

cat PredictedPass.gff3 | grep -w 'gene' | wc -l

12151: 10832
12158: 10387
12253: 10850
12161: 10360


##Remove those lines containing 'CUFF.9944.1.116' because we are not very confident for the cuff results

the command used was : cd /gene_pred/codingquary1/V.dahliae//additional: cat additional_genes.gff | grep -v -w 'CUFF.9944.1.116' > new_additional_genes.gff

Note that at this stage all codingquary genes contain . characters rather than _ characters


```bash
for Assembly in $(ls repeat_masked/*/61/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=gene_pred/final_genes/$Organism/$Strain/edited
mkdir -p $OutDir
BrakerGff=$(ls gene_pred/braker/$Organism/$Strain/V.dahliae_*_braker_two/augustus.gff3)
CodingQuaryGff=$(ls gene_pred/codingquary1/$Organism/$Strain/out/PredictedPass.gff3)
PGNGff=$(ls gene_pred/codingquary1/$Organism/$Strain/out/PGN_predictedPass.gff3)
cp -i $BrakerGff $OutDir/final_genes_Braker_ed.gff3
cat $CodingQuaryGff | grep -v -w -e 'CUFF.3884.1.3' -e 'CUFF.3883.1.2' > $OutDir/PredictedPass_ed.gff3
cp -i $PGNGff $OutDir/PGN_predictedPass_ed.gff3
echo ""
done
```

cat PredictedPass_ed.gff3 | grep -w 'gene' | wc -l

12251: 10831
12253: 10849
12158: 10387
12161: 10358

To make the loop: (not done)
    ```bash
for Assembly in $(ls repeat_masked/*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
if [ $Strain == 51 ]; then
echo ""
cat $CodingGff | grep -v 'CUFF.4517.1.10' > $CodingGffEd
elif [ $Strain == 53 ]; then
echo ""
cat $CodingGff | grep -v 'CCUFF.4072.1.5' > $CodingGffEd
fi [ $Strain == 61 ]; then
echo ""
cat $CodingGff | grep -v 'CUFF.3883.1.2' > $CodingGffEd
done
    ```

#Then additional transcripts were added to Braker gene models, when CodingQuary

genes were predicted in regions of the genome, not containing Braker gene models:

```bash
for EditDir in $(ls -d gene_pred/final_genes/*/*/edited); do
Strain=$(echo $EditDir | rev | cut -d '/' -f2 | rev)
Organism=$(echo $EditDir | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
BrakerGff=$EditDir/final_genes_Braker_ed.gff3
CodingQuaryGff=$EditDir/PredictedPass_ed.gff3
PGNGff=$EditDir/PGN_predictedPass_ed.gff3
# ManGff=$EditDir/manual_annotations.gff3
AddDir=$EditDir/additional
FinalDir=gene_pred/final_genes/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
# FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
# cp $ManGff $FinalDir/final_genes_manual.gff3
# $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_manual.gff3 $FinalDir/final_genes_manual
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
# GffManual=$FinalDir/final_genes_manual.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
# cat $GffBraker $GffQuary $GffManual > $GffAppended
done
  ```

The final number of genes per isolate was observed using:

```bash
for DirPath in $(ls -d  gene_pred/final_genes/V.dahliae/*/final); do
  echo $DirPath;
  cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
  cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
  cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
  echo "";
  done
  ```
gene_pred/final_genes/V.dahliae/51/final
9655
732
10387

gene_pred/final_genes/V.dahliae/53/final
9718
716
10434

gene_pred/final_genes/V.dahliae/58/final
9458
0
9458

gene_pred/final_genes/V.dahliae/61/final
9457
550
10007
