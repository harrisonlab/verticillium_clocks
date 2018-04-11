# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the F. oxysporum f. sp.
narcissi genome for submisson are shown below.

# Preliminary submission

A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov

A preliminary submission was made for the .fasta assembly to check if
any contigs needed to be split. This step was performed early in the annotation
process (prior to gene prediction) to ensure that annotation did not have to
be repeated at the end of the project.


The following note was provided in the WGS submission page on NCBI in the box
labeled "Private comments to NCBI staff":

```
I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'

# Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  ls genome_submission/SRA_metadata_acc.txt genome_submission/SRA_metadata_acc.xlsx
```

As these files included a file > 500 Mb, a presubmission folder was requested.
This aids submission of large data files. This file was created on the ftp server
at ftp-private.ncbi.nlm.nih.gov, with a private folder named
uploads/andrew.armitage@emr.ac.uk_6L2oakBI. Ncbi provided a username a password.
Files were uploaded into a folder created within my preload folder using ftp.

```bash
# Bioproject="PRJNA352681"
	SubFolder="Vd_PRJNA352681"
	mkdir $SubFolder
	for Read in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
	  echo $Read;
	  cp $Read $SubFolder/.
	done

ftp ftp-private.ncbi.nlm.nih.gov
cd uploads/emma.cascant-lopez@emr.ac.uk_F9JQ3w9M
mkdir Vd_PRJNA352681
cd Vd_PRJNA352681
# put Vd_PRJNA352681
prompt
mput *
bye
cd ../
rm -r $SubFolderssh
```

For WGS
# Calculate coverage using the count_nucl.pl script

## For one library

```bash
for DataDir in $(ls -d qc_dna/paired/*/*)
do
F_Read=$(ls $DataDir/F/*.gz)
R_Read=$(ls $DataDir/R/*.gz)
Strain=$(echo $DataDir | rev | cut -f1 -d '/' | rev)
Organism=$(echo $DataDir | rev | cut -f2 -d '/' | rev)
WorkDir=tmp_dir/$Strain
mkdir -p $WorkDir
cp -r $F_Read $WorkDir
cp -r $R_Read $WorkDir
cd $WorkDir
Read1=*F*
Read2=*R*
gunzip $Read1
gunzip $Read2
Sub1=*F*.fq
Sub2=*R*.fq
echo "$Organism - $Strain"
count_nucl.pl -i $Sub1 -i $Sub2 -g 35
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
done
```

WGS submission: SUB2105171 Vd12161
I would like to ask to NCBI to run the assembly thorough the contamination screen as this is not the final submission.
contigs_min_500bp.fasta


V.dahliae - 51
The estimated genome size is: 35000000 bp


The input file is: wilt_51_F_appended_trim.fq

Results for: wilt_51_F_appended_trim.fq
 Within this file of 1116178220 bp there were 6515637 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_51_R_appended_trim.fq

Results for: wilt_51_R_appended_trim.fq
 Within this file of 1111891504 bp there were 6515637 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 2228069724 nucleotides in this file.

 This equates to an estimated genome coverage of 63.66 .


V.dahliae - 53
The estimated genome size is: 35000000 bp


The input file is: wilt_53_F_appended_trim.fq

Results for: wilt_53_F_appended_trim.fq
 Within this file of 790141507 bp there were 3597693 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_53_R_appended_trim.fq

Results for: wilt_53_R_appended_trim.fq
 Within this file of 757410064 bp there were 3597693 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 1547551571 nucleotides in this file.

 This equates to an estimated genome coverage of 44.22 .


V.dahliae - 58
The estimated genome size is: 35000000 bp


The input file is: wilt_58_F_appended_trim.fq

Results for: wilt_58_F_appended_trim.fq
 Within this file of 1404169996 bp there were 6746472 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_58_R_appended_trim.fq

Results for: wilt_58_R_appended_trim.fq
 Within this file of 1385769457 bp there were 6746472 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 2789939453 nucleotides in this file.

 This equates to an estimated genome coverage of 79.71 .


V.dahliae - 61
The estimated genome size is: 35000000 bp


The input file is: wilt_61_F_appended_trim.fq

Results for: wilt_61_F_appended_trim.fq
 Within this file of 1087035735 bp there were 4990922 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_61_R_appended_trim.fq

Results for: wilt_61_R_appended_trim.fq
 Within this file of 1072541654 bp there were 4990922 fastq sequences
 of these 0 lines were empty.



Total results:

 There are a total of 2159577389 nucleotides in this file.

 This equates to an estimated genome coverage of 61.70 .


#Submission Commands

 Submisison of annotations with an assembly appears to be a complex process. If a genome is to be submitted without annotation then all that is needed is the fasta file containing the assembled contigs. If an annotated genome is to be submitted then a number of processing steps are required before submission. The fasta file of contigs and the gff file of annotations must be combined to form a .asn file. The program that does this conversion (tbl2asn) requires the fasta files and gff files to be formatted correctly. In the case of the gff file, this means parsing it to a .tbl file.

 The commands used to parse these files and prepare the Alternaria spp. genomes for submisson are shown below.

#Preliminary submission
 A Bioproject and biosample number was prepared for the genome submission at: https://submit.ncbi.nlm.nih.gov

A preliminary submission was made for the .fasta assembly to check if any contigs needed to be split. This step was performed early in the annotation process (prior to gene prediction) to ensure that annotation did not have to be repeated at the end of the project.

The following note was provided in the WGS submission page on NCBI in the box labeled "Private comments to NCBI staff":

```
I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'
```

#Making a table for locus tags:

locus tags were provided by ncbi when the bioproject was registered.

A table detailing their relationship to the strain was made manually. This could be searched later to extract locus tags for particular strains.

```bash
mkdir -p genome_submission/
printf \
"VD0001 SAMN05924275 12158
VD0002 SAMN05924795 12253
VD0003 SAMN05929039 12251
VD0004 SAMN05929042 12161" \
> genome_submission/Vd_PRJNA352681_locus_tags.txt
```

# Final Submission

 These commands were used in the final submission of Alternaria spp. genomes:


 ## Output directory
 An output and working directory was made for genome submission:

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
  echo "$Organism - $Strain"
  ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
  cd $ProjDir
  OutDir="genome_submission/$Organism/$Strain"
  mkdir -p $OutDir
done
```

## SbtFile
The genbank submission template tool was used at:
http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
This produce a template file detailing the submission.

## Setting varibales
Vairables containing locations of files and options for scripts were set:


# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/armita/git_repos/emr_repos/tools/genbank_submission"
# File locations:
SbtFile="genome_submission/template.sbt"
LabID="HarrisonEMR"

## Generating .tbl file (GAG)

The Genome Annotation Generator (GAG.py) can be used to convert gff files into
.tbl format, for use by tbl2asn.

It can also add annotations to features as provided by Annie the Annotation
extractor.

### Extracting annotations (Annie)

Interproscan and Swissprot annotations were extracted using annie, the
ANNotation Information Extractor. The output of Annie was filtered to
keep only annotations with references to ncbi approved databases.
Note - It is important that transcripts have been re-labelled as mRNA by this
point.

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  GffFile=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/final/$Organism/"$Strain"*/final/final_genes_appended_renamed.gff3)

  InterProTab=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/"$Strain"*/"$Strain"*_interproscan.tsv)
  SwissProtBlast=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/$Organism/"$Strain"*/swissprot_vJul2016_tophit_parsed.tbl)
  SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
  python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
  ProgDir=/home/armita/git_repos/emr_repos/tools/genbank_submission
  $ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv
done
```

#Running GAG

Gag was run using the modified gff file as well as the annie annotation file. Gag was noted to output database references incorrectly, so these were modified.

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir="genome_submission/$Organism/$Strain"
GffFile=$(ls gene_pred/final_genes/$Organism/"$Strain"*/final/final_genes_appended.gff3)
mkdir -p $OutDir/gag/round1_AA
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1_AA 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1_AA/genome.tbl
done
```

#tbl2asn round 1

tbl2asn was run an initial time to collect error reports on the current formatting of the .tbl file. Note - all input files for tbl2asn need to be in the same directory and have the same basename.

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir="genome_submission/$Organism/$Strain"

cp $Assembly $OutDir/gag/round1_AA/genome.fsa
SbtFile=$(ls genome_submission/template.sbt)
cp $SbtFile $OutDir/gag/round1_AA/genome.sbt
mkdir -p $OutDir/tbl2asn/round1_AA
tbl2asn -p $OutDir/gag/round1_AA/. -t $OutDir/gag/round1_AA/genome.sbt -r $OutDir/tbl2asn/round1_AA -M n -X E -Z $OutDir/gag/round1_AA/discrep.txt -j "[organism=$Organism] [strain=$Strain]"
done
```

#Editing .tbl file

The tbl2asn .val output files were observed and errors corrected. This was done with an in house script. The .val file indicated that some cds had premature stops, so these were marked as pseudogenes ('pseudo' - SEQ_FEAT.InternalStop) and that some genes had cds coordinates that did not match the end of the gene if the protein was hanging off a contig ('stop' - SEQ_FEAT.NoStop). Furthermore a number of other edits were made to bring the .tbl file in line with ncbi guidelines. This included: Marking the source of gene predictions and annotations ('add_inference'); Correcting locus_tags to use the given ncbi_id ('locus_tag'); Correcting the protein and transcript_ids to include the locus_tag and reference to submitter/lab id ('lab_id'), removal of annotated names of genes if you don't have high confidence in their validity (--gene_id 'remove'). If 5'-UTR and 3'-UTR were not predicted during gene annotation then genes, mRNA and exon features need to reflect this by marking them as incomplete ('unknown_UTR').

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  SubmissionID=$(cat genome_submission/Vd_PRJNA352681_locus_tags.txt | grep "$Strain" | cut -f1 -d ' ')
  echo $SubmissionID
  mkdir -p $OutDir/gag/edited_AA
  # ProgDir=/home/lopeze/git_repos/tools/genbank_submission
  ProgDir=/home/armita/git_repos/emr_repos/tools/genbank_submission
  $ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1_AA/genome.tbl --inp_val $OutDir/tbl2asn/round1_AA/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --rename_genes "g" --remove_product_locus_tags "True" --out_tbl $OutDir/gag/edited_AA/genome.tbl
done
```

#Generating a structured comment detailing annotation methods

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
  Annotation Provider\tHarrison Lab NIAB-EMR
  Annotation Date\tAUG-2017
  Annotation Version\tRelease 1.01
  Annotation Method\tAb initio gene prediction: Braker 1.9 and CodingQuary 2.0; Functional annotation: Swissprot (July 2016 release) and Interproscan 5.18-57.0" \
  > $OutDir/gag/edited_AA/annotation_methods.strcmt.txt
done
```

#Final run of tbl2asn

Following correction of the GAG .tbl file, tbl2asn was re-run to provide the final genbank submission file.

The options -l paired-ends -a r10k inform how to handle runs of Ns in the sequence, these options show that paired-ends have been used to estimate gaps and that runs of N's longer than 10 bp should be labelled as gaps.

```bash
for Assembly in $(ls repeat_masked/V.dahliae/*/ncbi_filtered_contigs_repmask/*_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  cp $Assembly $OutDir/gag/edited_AA/genome.fsa
  cp $SbtFile $OutDir/gag/edited_AA/genome.sbt
  mkdir $OutDir/tbl2asn/final
  tbl2asn -p $OutDir/gag/edited_AA/. -t $OutDir/gag/edited_AA/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$Organism] [strain=$Strain]" -l paired-ends -a r10k -w $OutDir/gag/edited_AA/annotation_methods.strcmt.txt
  cat $OutDir/tbl2asn/final/genome.sqn | sed 's/_pilon//g' | sed 's/\. subunit/kDa subunit/g' | sed 's/, mitochondrial//g' > $OutDir/tbl2asn/final/$FinalName.sqn
done
```
