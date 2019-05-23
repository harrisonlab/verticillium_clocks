# Optimise Codons

Commands used to generate a table for codon optimisation using RNAseq data
described in verticillium_clocks/README_RNAseq.md


## Locate input data

Work was performed in the directory:

```bash
  cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
```

```bash
# Light/dark experiment at 6hrs for wild type and Wc-1 KO line for isolate 120053  
ReadCounts=$(ls ../clocks/RNA_alignment/featureCounts/experiment_all/LD/countData_all_LD)
# Reference genome gene models
GeneModels=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
FastaCDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
```

Note - the Gff file includes ncrna features as "genes" that do not appear in the cds.




## Optimisation using all genes:




Remove genes with transposon annotations:
IPR000477 - reverse transriptase
IPR004875 - DDE superfamily endonuclease
IPR025476 - Helitron helicase-like domain
IPR012337 - Ribonuclease H-like superfamily

Remove genes with no annotations.
  (these may represent transposons)

```bash
OutDir=analysis/codon_usage/all_cds
mkdir -p $OutDir

CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
GeneList=$OutDir/$(basename ${CDS%.fa})._headers.txt
AnnotTab=$(ls ../clocks/gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_exp.tsv)

cat $CDS | grep '>' | tr -d '>' | sed 's/ $//g'> $GeneList

FilteredList=${GeneList%.txt}_no_transposons.txt
OutFasta=${GeneList%.txt}_no_transposons.fa

# make a list of genes with transposon annotations
cat $AnnotTab | grep -e 'IPR000477' -e 'IPR004875' -e 'IPR025476' -e 'IPR012337' -e 'transpos' > $OutDir/putative_transposon_IDs.txt
# make a list of genes with no interproscan annotations
cat $AnnotTab | cut -f1,29 | grep -P -v "\t\w" | cut -f1 > $OutDir/no_annotation_IDs.txt

cat $GeneList | grep ".t1$" | grep -v -f $OutDir/putative_transposon_IDs.txt | grep -v -f $OutDir/no_annotation_IDs.txt > $FilteredList
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $FilteredList > $OutFasta

cat $GeneList | wc -l
cat $OutFasta | grep '>' | wc -l
```

```
  11424
  8955
```

```bash
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir
OutDir=analysis/codon_usage/all_cds
FastaCDS=$(ls analysis/codon_usage/all_cds/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed._headers_no_transposons.fa)
mkdir $OutDir
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
$ProgDir/codonw.sh $(basename $FastaCDS)
cd $ProjDir
ls $OutDir/*/codon.coa
```

Make a codon usage table to be usable by optimizer:
http://genomes.urv.es/OPTIMIZER/tutorial.php

```bash
OutDir=analysis/codon_usage/all_cds
CodonFile=$(ls $OutDir/*/*_pass.cutot)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/codonw2optimizer.py --codonw $CodonFile
```

Vd codon optimisation table based on total genes:

```
UUU: 0.64; UCU: 0.79; UAU: 0.50; UGU: 0.53; UUC: 1.36; UCC: 1.30; UAC: 1.50; UGC: 1.47; UUA: 0.11; UCA: 0.61; UAA: 0.69; UGA: 1.35; UUG: 0.64; UCG: 1.36; UAG: 0.96; UGG: 1.00; CUU: 0.82; CCU: 0.84; CAU: 0.61; CGU: 0.78; CUC: 2.36; CCC: 1.48; CAC: 1.39; CGC: 2.27; CUA: 0.29; CCA: 0.65; CAA: 0.57; CGA: 0.79; CUG: 1.78; CCG: 1.03; CAG: 1.43; CGG: 0.86; AUU: 0.90; ACU: 0.54; AAU: 0.47; AGU: 0.44; AUC: 1.86; ACC: 1.42; AAC: 1.53; AGC: 1.50; AUA: 0.23; ACA: 0.69; AAA: 0.38; AGA: 0.46; AUG: 1.00; ACG: 1.34; AAG: 1.62; AGG: 0.84; GUU: 0.69; GCU: 0.73; GAU: 0.63; GGU: 0.76; GUC: 2.06; GCC: 1.86; GAC: 1.37; GGC: 2.24; GUA: 0.24; GCA: 0.53; GAA: 0.56; GGA: 0.53; GUG: 1.00; GCG: 0.88; GAG: 1.44; GGG: 0.47
```


## Analysis of highly expressed genes:

```bash
  mkdir analysis/codon_usage/by_expression
```

### Align RNAseq data to the Vd genome


# RNAseq

Perform RNAseq pseudo-alignment using Salmon

```bash
for Transcriptome in $(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../clocks/qc_rna/experiment*/V.dahliae/*); do
FileNum=$(ls $RNADir/F/*.fq.gz | wc -l)
for num in $(seq 1 $FileNum); do
Jobs=$(qstat | grep 'sub_salmon' | grep 'qw' | wc -l)
while [ $Jobs -gt 0 ]; do
sleep 2
printf "."
Jobs=$(qstat | grep 'sub_salmon' | grep 'qw' | wc -l)
done
printf "\n"
FileF=$(ls $RNADir/F/*.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed 's/_1_trim.fq.gz//g')
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint - $Prefix"
OutDir=alignment/salmon/$Organism/"$Strain"/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
done
```

Combine Salmon quasi-quantifications into a single file:


```bash
OutDir=alignment/salmon/TPM
mkdir -p $OutDir

for File in $(ls alignment/salmon/*/*/*/*/quant.sf | head -n1); do
  cat $File | cut -f1,2,3 > $OutDir/gene_info.tsv
done
# Put files in a convenient location for DeSeq. Analysis was not performed on
# Strawberry control samples.
for File in $(ls alignment/salmon/*/*/*/*/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  printf "${Prefix}\n" > $OutDir/${Prefix}_TPM.txt
  cat $File | cut -f4 | tail -n+2 >> $OutDir/${Prefix}_TPM.txt
done
paste $OutDir/gene_info.tsv $OutDir/*_TPM.txt > $OutDir/gene_TPM.tsv
```

The file of TPM was downloaded to my local machine (actually EMQA) for analysis using R studio.

```bash
EmqaDir="/Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation"
scp cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/alignment/salmon/TPM/gene_TPM.tsv $EmqaDir/.
```

```R

library("ggplot2")
gene_TPM <- read.delim("/Volumes/GGB/Harrison/Projects/BBSRC IPA Clocks/Science/WP2. Light, temperature, clock role in pathogenicity/Obj4. HIGS/codon_optimisation/gene_TPM.tsv")

summaryfunc <- function(prefix,x) {
  o1 <- paste(prefix, "summary.txt", sep="_")
  sink(o1)
  print("Summary of all genes")
  y <- summary(x)
  print(y)
  expressed_x <- x[x > 0]
  print("Summary of genes with TPM > 0")
  z <- summary(expressed_x)
  print(z)
  sink()
}
histfunc <- function(x){
  # hist(gene_TPM[x])
  # remove counts < 1
  sub_x <- subset(gene_TPM, x > 0)
  # calculate upper quantiles of subset genes:
  quantiles <- quantile(gene_TPM[,x], c(.25, .75, .90, .95))
  low <- round(as.numeric(quantiles[1]))
  quartile <- round(as.numeric(quantiles[2]))
  decile <- round(as.numeric(quantiles[3]))
  # p <- ggplot(data=sub_x, aes(sub_x[,x])) +
  #   geom_histogram(bins=100)
  # # p <- p + scale_x_log10()
  # p <- p + xlim(1, 500)
  # p <- p + geom_vline(xintercept=quartile, linetype="dashed", color = "black")
  # p <- p + geom_text(aes(x=quartile, label=paste(quartile, "top 25%", sep="\n"), y=0), colour="black", angle=0, text=element_text(size=6))
  # p <- p + geom_vline(xintercept=decile, linetype="dashed", color = "black")
  # p <- p + geom_text(aes(x=decile, label=paste(decile, "top 10%", sep="\n"), y=0), colour="black", angle=0, text=element_text(size=6))
  # o2 <- paste(x, "hist.png", sep="_")
  # ggsave(o2, plot = p)
  q_ids <- gene_TPM$Name[gene_TPM[,x] >= quartile]
  o3 <- paste(x, "quartile_ids.txt", sep="_")
  write(as.character(q_ids), file=o3)
  d_ids <- gene_TPM$Name[gene_TPM[,x] >= decile]
  o4 <- paste(x, "decile_ids.txt", sep="_")
  write(as.character(d_ids), file=o4)
  low_ids <- gene_TPM$Name[gene_TPM[,x] <= low]
  o4 <- paste(x, "low_ids.txt", sep="_")
  write(as.character(low_ids), file=o4)
}


lapply(colnames(gene_TPM[4:length(gene_TPM)]), function(x){
  summaryfunc(x, gene_TPM[x])
  })
lapply(colnames(gene_TPM[4:length(gene_TPM)]), function(x){
  histfunc(x)
  })
```

Collate this into a single list of highly expressed genes accross all timepoints:

```bash
EmqaDir='/Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio'
cat "$EmqaDir/R-studio/*decile_ids.txt" | sort | uniq -c | grep "60.*VDAG" | sed "s/*.VDAG/VDAG/g" > $EmqaDir/decile_genes_filtered.txt

cat "$EmqaDir/R-studio/*decile_ids.txt" | sort | uniq > $EmqaDir/decile_genes_combined.txt

cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*quartile_ids.txt | sort | uniq > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/quartile_genes_combined.txt

cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*decile_ids.txt | sort | uniq > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/decile_genes_combined.txt

cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*low_ids.txt | sort | uniq > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/low_genes_combined.txt

cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*quartile_ids.txt | sort | uniq -c | sort -nr | grep "57.*VDAG" | sed "s/*.VDAG/VDAG/g" > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/quartile_genes_filtered.txt

cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*decile_ids.txt | sort | uniq -c | sort -nr | grep "57.*VDAG" | sed "s/*.VDAG/VDAG/g" > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/decile_genes_filtered.txt


cat /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/R-studio/*low_ids.txt | sort | uniq -c | sort -nr | grep "57.*VDAG" | sed "s/*.VDAG/VDAG/g" > /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/low_genes_filtered.txt
```

```bash
scp /Volumes/GGB/Harrison/Projects/BBSRC\ IPA\ Clocks/Science/WP2.\ Light\,\ temperature\,\ clock\ role\ in\ pathogenicity/Obj4.\ HIGS/codon_optimisation/*.txt cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/codon_usage/by_expression/.
```

### Extract highly expressed genes

```bash
CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
GeneList=$(ls analysis/codon_usage/by_expression/quartile_genes_filtered.txt)
OutFasta=${GeneList%.txt}.fa
# A single transcript should be used from each gene. Check this assumption:
cat $GeneList | grep -v ".t1$" | wc -l
# cat $CDS | grep -w -f $OutDir/quorn_fpkm_148_IDs.txt > $OutDir/quorn_fpkm_148_IDs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $GeneList > $OutFasta
```

Remove genes with transposon annotations:
IPR000477 - reverse transriptase
IPR004875 - DDE superfamily endonuclease
IPR025476 - Helitron helicase-like domain
IPR012337 - Ribonuclease H-like superfamily

Remove genes with no annotations.
  (these may represent transposons)

```bash
OutDir=analysis/codon_usage/by_expression
GeneList=$(ls analysis/codon_usage/by_expression/quartile_genes_filtered.txt)
AnnotTab=$(ls ../clocks/gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_exp.tsv)
CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
FilteredList=${GeneList%.txt}_no_transposons.txt
OutFasta=${GeneList%.txt}_no_transposons.fa

# make a list of genes with transposon annotations
cat $AnnotTab | grep -e 'IPR000477' -e 'IPR004875' -e 'IPR025476' -e 'IPR012337' -e 'transpos' > $OutDir/putative_transposon_IDs.txt
# make a list of genes with no interproscan annotations
cat $AnnotTab | cut -f1,29 | grep -P -v "\t\w" | cut -f1 > $OutDir/no_annotation_IDs.txt

cat $GeneList | grep -v -f $OutDir/putative_transposon_IDs.txt | grep -v -f $OutDir/no_annotation_IDs.txt > $FilteredList
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $FilteredList > $OutFasta

cat $OutFasta | grep '>' | wc -l
```




```bash
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
OutDir=analysis/codon_usage/by_expression
cd $ProjDir
FastaCDS=$(ls $ProjDir/$OutDir/quartile_genes_filtered_no_transposons.fa)
mkdir -p $OutDir
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
$ProgDir/codonw.sh $(basename $FastaCDS)
cd $ProjDir

ls $ProjDir/$OutDir/*/codon.coa
```

Make a codon usage table to be usable by optimizer:
http://genomes.urv.es/OPTIMIZER/tutorial.php

```bash
CodonFile=$(ls $OutDir/*/*_pass.cutot)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/codonw2optimizer.py --codonw $CodonFile
```

Vd codon optimisation table based on highly expressed genes:

```
  UUU: 0.55; UCU: 0.85; UAU: 0.40; UGU: 0.44; UUC: 1.45; UCC: 1.41; UAC: 1.60; UGC: 1.56; UUA: 0.07; UCA: 0.51; UAA: 1.04; UGA: 1.06; UUG: 0.55; UCG: 1.36; UAG: 0.90; UGG: 1.00; CUU: 0.84; CCU: 0.93; CAU: 0.50; CGU: 0.97; CUC: 2.49; CCC: 1.66; CAC: 1.50; CGC: 2.45; CUA: 0.23; CCA: 0.53; CAA: 0.50; CGA: 0.75; CUG: 1.81; CCG: 0.89; CAG: 1.50; CGG: 0.71; AUU: 0.95; ACU: 0.58; AAU: 0.37; AGU: 0.39; AUC: 1.91; ACC: 1.62; AAC: 1.63; AGC: 1.49; AUA: 0.14; ACA: 0.58; AAA: 0.26; AGA: 0.36; AUG: 1.00; ACG: 1.22; AAG: 1.74; AGG: 0.76; GUU: 0.73; GCU: 0.80; GAU: 0.61; GGU: 0.91; GUC: 2.20; GCC: 1.98; GAC: 1.39; GGC: 2.32; GUA: 0.18; GCA: 0.44; GAA: 0.48; GGA: 0.46; GUG: 0.88; GCG: 0.77; GAG: 1.52; GGG: 0.31
```


## Calculate delta RSCU between codon sets:

High extracted genes were already extracted. Low extracted genes needed to be extracted:

### Extract low expressed genes


```bash
CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
cat $CDS >
GeneList=$(ls analysis/codon_usage/by_expression/low_genes_filtered.txt)
OutFasta=${GeneList%.txt}.fa
# A single transcript should be used from each gene. Check this assumption:
cat $GeneList | grep -v ".t1$" | wc -l
# cat $CDS | grep -w -f $OutDir/quorn_fpkm_148_IDs.txt > $OutDir/quorn_fpkm_148_IDs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $GeneList > $OutFasta
```

Remove genes with transposon annotations:
IPR000477 - reverse transriptase
IPR004875 - DDE superfamily endonuclease
IPR025476 - Helitron helicase-like domain
IPR012337 - Ribonuclease H-like superfamily

Remove genes with no annotations.
  (these may represent transposons)

```bash
OutDir=analysis/codon_usage/by_expression
GeneList=$(ls analysis/codon_usage/by_expression/low_genes_filtered.txt)
AnnotTab=$(ls ../clocks/gene_pred/annotation/V.dahliae/JR2/JR2_gene_table_incl_exp.tsv)
CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
FilteredList=${GeneList%.txt}_no_transposons.txt
OutFasta=${GeneList%.txt}_no_transposons.fa

# make a list of genes with transposon annotations
cat $AnnotTab | grep -e 'IPR000477' -e 'IPR004875' -e 'IPR025476' -e 'IPR012337' -e 'transpos' > $OutDir/putative_transposon_IDs.txt
# make a list of genes with no interproscan annotations
cat $AnnotTab | cut -f1,29 | grep -P -v "\t\w" | cut -f1 > $OutDir/no_annotation_IDs.txt

cat $GeneList | sed "s/.*VDAG/VDAG/g" | grep -v -f $OutDir/putative_transposon_IDs.txt | grep -v -f $OutDir/no_annotation_IDs.txt > $FilteredList
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $FilteredList > $OutFasta

cat $OutFasta | grep '>' | wc -l
```

```
1114
```

### Calculate delta RSCU

```bash
OutDir=analysis/codon_usage/by_expression/delta_rscu
mkdir -p $OutDir
CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
cat $CDS | cut -f1 -d ' ' > $OutDir/Vd_CDS_parsed.fa

Low=analysis/codon_usage/by_expression/low_genes_filtered_no_transposons.txt
High=analysis/codon_usage/by_expression/quartile_genes_filtered_no_transposons.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/codon_usage
$ProgDir/d_rscu.R --low $Low --high $High --cds $OutDir/Vd_CDS_parsed.fa | less
```


## Optimise codons for given proteins

### GFP

```bash
OutDir=analysis/codon_usage/by_expression/optimised
mkdir -p $OutDir
printf ">GFP\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSKLSKDPNEKRDHMVLLEFVTAAGITLGMDELYK\n" > $OutDir/GFP_unoptimised.aa

printf ">human_optimised_GFP\nATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG\n" > $OutDir/GFP_unoptimised.cds
printf ">Vd_HC_optimizer_GFP\nATGGTGAGCAAGGGGGAGGAGCTTTTCACGGGAGTCGTGCCAATTCTCGTGGAGCTTGACGGCGACGTGAACGGCCACAAATTCAGCGTGTCTGGAGAGGGAGAGGGCGACGCGACGTACGGAAAGCTGACGCTCAAGTTCATTTGCACCACCGGCAAGCTTCCTGTCCCGTGGCCGACACTGGTTACAACGCTGACGTACGGCGTCCAGTGCTTCTCCCGATACCCCGACCACATGAAGCAACATGATTTTTTTAAGTCGGCTATGCCCGAGGGCTACGTCCAAGAGCGTACCATTTTCTTCAAGGATGATGGAAACTACAAGACCCGAGCGGAAGTCAAATTCGAGGGCGACACGCTCGTCAACAGGATTGAGTTGAAGGGGATCGATTTCAAGGAAGATGGCAATATCCTCGGCCACAAGCTCGAATACAACTACAACAGCCACAACGTCTACATTATGGCGGACAAGCAGAAGAATGGAATCAAGGTTAACTTTAAAATCCGCCACAACATTGAAGATGGAAGTGTCCAGCTCGCGGACCATTACCAGCAGAACACGCCTATCGGCGATGGCCCCGTGCTCCTCCCGGACAACCACTATCTTTCGACCCAGTCAAAGTTGAGTAAGGACCCCAACGAGAAGCGCGACCATATGGTATTGCTGGAATTCGTGACCGCGGCAGGCATCACTCTGGGGATGGACGAGTTGTACAAG\n" >> $OutDir/GFP_unoptimised.cds

printf 'UUU: 0.55; UCU: 0.85; UAU: 0.40; UGU: 0.44; UUC: 1.45; UCC: 1.41; UAC: 1.60; UGC: 1.56; UUA: 0.07; UCA: 0.51; UAA: 1.04; UGA: 1.06; UUG: 0.55; UCG: 1.36; UAG: 0.90; UGG: 1.00; CUU: 0.84; CCU: 0.93; CAU: 0.50; CGU: 0.97; CUC: 2.49; CCC: 1.66; CAC: 1.50; CGC: 2.45; CUA: 0.23; CCA: 0.53; CAA: 0.50; CGA: 0.75; CUG: 1.81; CCG: 0.89; CAG: 1.50; CGG: 0.71; AUU: 0.95; ACU: 0.58; AAU: 0.37; AGU: 0.39; AUC: 1.91; ACC: 1.62; AAC: 1.63; AGC: 1.49; AUA: 0.14; ACA: 0.58; AAA: 0.26; AGA: 0.36; AUG: 1.00; ACG: 1.22; AAG: 1.74; AGG: 0.76; GUU: 0.73; GCU: 0.80; GAU: 0.61; GGU: 0.91; GUC: 2.20; GCC: 1.98; GAC: 1.39; GGC: 2.32; GUA: 0.18; GCA: 0.44; GAA: 0.48; GGA: 0.46; GUG: 0.88; GCG: 0.77; GAG: 1.52; GGG: 0.31' > $OutDir/Vd_codon_weights.cds

# CDS=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa)
# cat $CDS | cut -f1 -d ' ' > $OutDir/Vd_CDS_parsed.fa

ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/codon_optimisation
$ProgDir/score_codons.py --fasta_aa $OutDir/GFP_unoptimised.aa  --fasta_cds $OutDir/GFP_unoptimised.cds --codon_table $OutDir/Vd_codon_weights.cds --prefix $OutDir/GFP
```

```
randomised codons:
AUGGUGUCAAAGGGCGAGGAGCUAUUCACGGGGGUGGUACCAAUCCUAGUGGAACUAGACGGGGACGUGAACGGACACAAGUUCAGCGUGAGCGGGGAGGGGGAGGGAGAUGCGACGUAUGGCAAGCUUACAUUGAAGUUCAUCUGUACCACGGGGAAGCUACCGGUGCCGUGGCCUACUCUAGUGACAACGCUCACAUACGGAGUGCAGUGCUUCAGCCGGUAUCCGGACCACAUGAAGCAGCACGACUUCUUCAAGUCGGCCAUGCCGGAGGGGUACGUGCAAGAGCGAACCAUCUUUUUUAAGGACGACGGGAACUACAAGACGAGGGCCGAGGUGAAGUUCGAGGGGGACACGCUGGUAAACCGCAUCGAGCUUAAGGGAAUCGACUUCAAGGAAGAUGGAAAUAUACUUGGGCACAAGCUGGAGUACAACUACAACUCGCAUAACGUGUACAUAAUGGCCGACAAGCAGAAGAACGGAAUAAAGGUGAACUUCAAGAUAAGGCACAACAUAGAGGACGGGUCCGUGCAGCUAGCGGACCACUACCAGCAGAACACACCGAUAGGCGACGGACCCGUCCUGCUGCCCGACAACCACUACCUAAGUACGCAGUCAAAGCUUAGCAAGGAUCCGAACGAAAAGCGCGACCACAUGGUACUUCUUGAGUUCGUUACAGCGGCGGGAAUCACACUAGGUAUGGAUGAGCUCUACAAG
0.74 of 1.0
optimum sequence:
AUGGUCAGCAAGGGCGAGGAGCUCUUCACCGGCGUCGUCCCCAUCCUCGUCGAGCUCGACGGCGACGUCAACGGCCACAAGUUCAGCGUCAGCGGCGAGGGCGAGGGCGACGCCACCUACGGCAAGCUCACCCUCAAGUUCAUCUGCACCACCGGCAAGCUCCCCGUCCCCUGGCCCACCCUCGUCACCACCCUCACCUACGGCGUCCAGUGCUUCAGCCGCUACCCCGACCACAUGAAGCAGCACGACUUCUUCAAGAGCGCCAUGCCCGAGGGCUACGUCCAGGAGCGCACCAUCUUCUUCAAGGACGACGGCAACUACAAGACCCGCGCCGAGGUCAAGUUCGAGGGCGACACCCUCGUCAACCGCAUCGAGCUCAAGGGCAUCGACUUCAAGGAGGACGGCAACAUCCUCGGCCACAAGCUCGAGUACAACUACAACAGCCACAACGUCUACAUCAUGGCCGACAAGCAGAAGAACGGCAUCAAGGUCAACUUCAAGAUCCGCCACAACAUCGAGGACGGCAGCGUCCAGCUCGCCGACCACUACCAGCAGAACACCCCCAUCGGCGACGGCCCCGUCCUCCUCCCCGACAACCACUACCUCAGCACCCAGAGCAAGCUCAGCAAGGACCCCAACGAGAAGCGCGACCACAUGGUCCUCCUCGAGUUCGUCACCGCCGCCGGCAUCACCCUCGGCAUGGACGAGCUCUACAAG
1.0 of 1.0
worst sequence:
AUGGUAAGUAAAGGGGAAGAAUUAUUUACAGGGGUAGUACCAAUAUUAGUAGAAUUAGAUGGGGAUGUAAAUGGGCAUAAAUUUAGUGUAAGUGGGGAAGGGGAAGGGGAUGCAACAUAUGGGAAAUUAACAUUAAAAUUUAUAUGUACAACAGGGAAAUUACCAGUACCAUGGCCAACAUUAGUAACAACAUUAACAUAUGGGGUACAAUGUUUUAGUAGAUAUCCAGAUCAUAUGAAACAACAUGAUUUUUUUAAAAGUGCAAUGCCAGAAGGGUAUGUACAAGAAAGAACAAUAUUUUUUAAAGAUGAUGGGAAUUAUAAAACAAGAGCAGAAGUAAAAUUUGAAGGGGAUACAUUAGUAAAUAGAAUAGAAUUAAAAGGGAUAGAUUUUAAAGAAGAUGGGAAUAUAUUAGGGCAUAAAUUAGAAUAUAAUUAUAAUAGUCAUAAUGUAUAUAUAAUGGCAGAUAAACAAAAAAAUGGGAUAAAAGUAAAUUUUAAAAUAAGACAUAAUAUAGAAGAUGGGAGUGUACAAUUAGCAGAUCAUUAUCAACAAAAUACACCAAUAGGGGAUGGGCCAGUAUUAUUACCAGAUAAUCAUUAUUUAAGUACACAAAGUAAAUUAAGUAAAGAUCCAAAUGAAAAAAGAGAUCAUAUGGUAUUAUUAGAAUUUGUAACAGCAGCAGGGAUAACAUUAGGGAUGGAUGAAUUAUAUAAA
0.27 of 1.0
midpoint sequence:
0.67
AUGGUAUCGAAGGGGGAGGAGUUGUUUACAGGGGUCGUGCCGAUACUAGUAGAACUUGAUGGAGAUGUCAACGGGCAUAAGUUUAGUGUAAGUGGGGAAGGGGAAGGCGACGCCACAUAUGGAAAGCUUACGCUGAAGUUUAUAUGUACGACAGGGAAGCUGCCGGUCCCGUGGCCGACACUGGUCACGACGCUGACAUAUGGCGUGCAGUGCUUUAGCCGCUAUCCGGAUCACAUGAAGCAACACGACUUUUUUAAGAGCGCGAUGCCUGAGGGGUACGUACAGGAGAGGACUAUAUUCUUCAAGGACGACGGGAACUACAAGACGAGGGCUGAGGUAAAGUUCGAGGGAGACACAUUGGUCAACAGAAUAGAGCUGAAAGGGAUAGACUUCAAGGAGGAUGGGAACAUACUGGGUCAUAAGCUCGAGUACAACUAUAACUCGCAUAACGUAUACAUAAUGGCGGACAAGCAAAAGAACGGGAUAAAGGUGAAUUUCAAGAUCAGGCACAACAUAGAGGAUGGAAGUGUGCAGCUAGCGGAUCACUACCAGCAGAACACCCCAAUAGGCGACGGGCCCGUGUUGCUGCCUGACAACCAUUACCUAUCGACGCAGAGCAAGCUGAGUAAGGACCCGAACGAGAAGCGAGACCAUAUGGUAUUGCUCGAAUUCGUAACAGCGGCCGGUAUCACUCUAGGGAUGGACGAAUUGUACAAG
Score of the pre-optimised sequence:
human_optimised_GFP
0.92 of 1.0
Vd_HC_optimizer_GFP
0.78 of 1.0
```

open an R session:

```bash
R
```

```R
GFP_1000_scores <- read.table("~/cluster_mount/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/codon_usage/by_expression/optimised/GFP_1000_scores.tsv", quote="\"", comment.char="")
library("ggplot2")
p <- ggplot(data=GFP_1000_scores, aes(GFP_1000_scores[,1])) + geom_histogram(bins = '50')
p <- p + xlim(0.3, 1.0)
p <- p + xlab("Total codon score") + ylab("Count")
p <- p + xlab("Total codon score") + ylab("Count")

p <- p + geom_vline(xintercept = 0.92, na.rm = FALSE, show.legend = NA)
p
```

### mCherry

```bash
OutDir=analysis/codon_usage/by_expression/optimised
mkdir -p $OutDir
printf ">GFP\nMVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK\n" > $OutDir/mCherry_unoptimised.aa

printf ">Vd_optimised_mCherry\nATGGTCTCCAAGGGCGAGGAGGACAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTCCACATGGAGGGCAGCGTCAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCTTACGAGGGCACCCAGACCGCCAAGCTCAAGGTCACGAAGGGCGGCCCTCTCCCGTTCGCCTGGGACATCCTCTCCCCCCAGTTCATGTACGGCTCCAAGGCCTACGTTAAGCACCCCGCCGACATTCCCGACTACCTCAAGCTCTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTCATGAACTTCGAGGACGGCGGCGTTGTCACCGTCACCCAGGACTCCTCCCTCCAGGACGGCGAGTTCATCTACAAGGTCAAGCTCCGCGGCACCAACTTCCCCTCCGACGGCCCTGTTATGCAGAAGAAGACCATGGGCTGGGAGGCCTCTAGCGAGCGCATGTACCCCGAGGACGGCGCCCTCAAGGGCGAGATCAAGCAGCGCCTCAAGCTCAAGGACGGCGGCCACTACGACGCCGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTCCAGCTCCCCGGCGCCTACAACGTTAACATCAAGCTGGACATCACCTCCCACAACGAGGACTACACGATTGTCGAGCAGTACGAGCGCGCCGAGGGCCGCCACTCTACCGGCGGCATGGACGAGCTCTACAAGTAG\n" > $OutDir/mCherry_unoptimised.cds

printf 'UUU: 0.55; UCU: 0.85; UAU: 0.40; UGU: 0.44; UUC: 1.45; UCC: 1.41; UAC: 1.60; UGC: 1.56; UUA: 0.07; UCA: 0.51; UAA: 1.04; UGA: 1.06; UUG: 0.55; UCG: 1.36; UAG: 0.90; UGG: 1.00; CUU: 0.84; CCU: 0.93; CAU: 0.50; CGU: 0.97; CUC: 2.49; CCC: 1.66; CAC: 1.50; CGC: 2.45; CUA: 0.23; CCA: 0.53; CAA: 0.50; CGA: 0.75; CUG: 1.81; CCG: 0.89; CAG: 1.50; CGG: 0.71; AUU: 0.95; ACU: 0.58; AAU: 0.37; AGU: 0.39; AUC: 1.91; ACC: 1.62; AAC: 1.63; AGC: 1.49; AUA: 0.14; ACA: 0.58; AAA: 0.26; AGA: 0.36; AUG: 1.00; ACG: 1.22; AAG: 1.74; AGG: 0.76; GUU: 0.73; GCU: 0.80; GAU: 0.61; GGU: 0.91; GUC: 2.20; GCC: 1.98; GAC: 1.39; GGC: 2.32; GUA: 0.18; GCA: 0.44; GAA: 0.48; GGA: 0.46; GUG: 0.88; GCG: 0.77; GAG: 1.52; GGG: 0.31' > $OutDir/Vd_codon_weights.cds

ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/codon_optimisation
$ProgDir/score_codons.py --fasta_aa $OutDir/mCherry_unoptimised.aa  --fasta_cds $OutDir/mCherry_unoptimised.cds --codon_table $OutDir/Vd_codon_weights.cds --prefix $OutDir/mCherry
```

```
randomised codons:
AUGGUGUCAAAGGGGGAGGAGGACAAUAUGGCAAUCAUCAAGGAAUUCAUGCGAUUCAAGGUGCACAUGGAGGGGUCUGUAAACGGGCACGAGUUCGAAAUUGAGGGCGAGGGGGAGGGGCGUCCAUACGAGGGAACCCAGACGGCCAAGCUAAAGGUGACGAAAGGGGGACCACUCCCGUUCGCCUGGGACAUAUUGUCACCACAGUUCAUGUAUGGGUCGAAGGCGUACGUCAAGCACCCGGCGGACAUUCCCGAUUAUCUAAAGCUGAGUUUCCCGGAGGGUUUCAAGUGGGAGCGGGUAAUGAACUUCGAGGACGGGGGGGUAGUGACGGUAACGCAGGACAGUUCACUACAGGACGGGGAAUUCAUAUACAAGGUGAAGCUGAGGGGAACGAACUUCCCGUCUGACGGGCCAGUGAUGCAAAAGAAGACGAUGGGAUGGGAGGCGAGUUCAGAGAGGAUGUACCCAGAGGACGGAGCGCUGAAAGGGGAGAUAAAGCAGCGGCUAAAGCUAAAGGAUGGUGGGCAUUACGACGCUGAGGUGAAGACUACGUAUAAGGCGAAGAAACCGGUACAGCUGCCGGGGGCGUACAACGUGAACAUCAAGCUGGACAUAACAAGUCACAACGAGGACUACACCAUCGUCGAGCAGUACGAGAGGGCGGAGGGGCGUCACUCUACAGGCGGCAUGGACGAGCUGUACAAG
0.73 of 1.0
optimum sequence:
AUGGUCAGCAAGGGCGAGGAGGACAACAUGGCCAUCAUCAAGGAGUUCAUGCGCUUCAAGGUCCACAUGGAGGGCAGCGUCAACGGCCACGAGUUCGAGAUCGAGGGCGAGGGCGAGGGCCGCCCCUACGAGGGCACCCAGACCGCCAAGCUCAAGGUCACCAAGGGCGGCCCCCUCCCCUUCGCCUGGGACAUCCUCAGCCCCCAGUUCAUGUACGGCAGCAAGGCCUACGUCAAGCACCCCGCCGACAUCCCCGACUACCUCAAGCUCAGCUUCCCCGAGGGCUUCAAGUGGGAGCGCGUCAUGAACUUCGAGGACGGCGGCGUCGUCACCGUCACCCAGGACAGCAGCCUCCAGGACGGCGAGUUCAUCUACAAGGUCAAGCUCCGCGGCACCAACUUCCCCAGCGACGGCCCCGUCAUGCAGAAGAAGACCAUGGGCUGGGAGGCCAGCAGCGAGCGCAUGUACCCCGAGGACGGCGCCCUCAAGGGCGAGAUCAAGCAGCGCCUCAAGCUCAAGGACGGCGGCCACUACGACGCCGAGGUCAAGACCACCUACAAGGCCAAGAAGCCCGUCCAGCUCCCCGGCGCCUACAACGUCAACAUCAAGCUCGACAUCACCAGCCACAACGAGGACUACACCAUCGUCGAGCAGUACGAGCGCGCCGAGGGCCGCCACAGCACCGGCGGCAUGGACGAGCUCUACAAG
1.0 of 1.0
worst sequence:
AUGGUAAGUAAAGGGGAAGAAGAUAAUAUGGCAAUAAUAAAAGAAUUUAUGAGAUUUAAAGUACAUAUGGAAGGGAGUGUAAAUGGGCAUGAAUUUGAAAUAGAAGGGGAAGGGGAAGGGAGACCAUAUGAAGGGACACAAACAGCAAAAUUAAAAGUAACAAAAGGGGGGCCAUUACCAUUUGCAUGGGAUAUAUUAAGUCCACAAUUUAUGUAUGGGAGUAAAGCAUAUGUAAAACAUCCAGCAGAUAUACCAGAUUAUUUAAAAUUAAGUUUUCCAGAAGGGUUUAAAUGGGAAAGAGUAAUGAAUUUUGAAGAUGGGGGGGUAGUAACAGUAACACAAGAUAGUAGUUUACAAGAUGGGGAAUUUAUAUAUAAAGUAAAAUUAAGAGGGACAAAUUUUCCAAGUGAUGGGCCAGUAAUGCAAAAAAAAACAAUGGGGUGGGAAGCAAGUAGUGAAAGAAUGUAUCCAGAAGAUGGGGCAUUAAAAGGGGAAAUAAAACAAAGAUUAAAAUUAAAAGAUGGGGGGCAUUAUGAUGCAGAAGUAAAAACAACAUAUAAAGCAAAAAAACCAGUACAAUUACCAGGGGCAUAUAAUGUAAAUAUAAAAUUAGAUAUAACAAGUCAUAAUGAAGAUUAUACAAUAGUAGAACAAUAUGAAAGAGCAGAAGGGAGACAUAGUACAGGGGGGAUGGAUGAAUUAUAUAAA
0.3 of 1.0
midpoint sequence:
0.68
AUGGUAUCUAAGGGGGAGGAAGAUAACAUGGCGAUAAUAAAAGAGUUCAUGCGAUUCAAGGUACACAUGGAGGGAUCAGUGAACGGCCAUGAGUUCGAGAUCGAGGGGGAGGGGGAGGGGAGACCAUAUGAGGGGACGCAGACGGCGAAGCUGAAGGUAACUAAGGGGGGGCCACUCCCAUUCGCAUGGGAUAUACUCUCGCCACAGUUUAUGUACGGGUCGAAGGCCUACGUGAAGCACCCAGCUGACAUACCAGAUUACCUCAAGCUGUCAUUCCCUGAGGGUUUUAAGUGGGAACGGGUAAUGAACUUCGAGGACGGGGGGGUAGUGACGGUAACGCAGGAUUCGUCACUACAGGACGGCGAAUUCAUAUAUAAGGUUAAGCUCCGCGGGACAAACUUUCCGAGCGACGGGCCGGUGAUGCAGAAGAAGACGAUGGGUUGGGAAGCAUCUAGUGAGAGAAUGUACCCGGAGGAUGGGGCGCUGAAGGGCGAGAUAAAACAGAGGCUGAAGCUGAAAGACGGGGGGCACUACGAUGCAGAAGUGAAGACAACUUACAAGGCUAAAAAGCCAGUGCAACUACCAGGGGCCUACAACGUGAACAUAAAACUCGACAUCACGUCGCAUAACGAGGACUAUACAAUCGUGGAGCAAUAUGAGAGGGCGGAAGGGAGGCACAGUACUGGUGGAAUGGACGAGCUCUACAAG
Score of the pre-optimised sequence:
Vd_optimised_mCherry
0.98 of 1.0
```

```R
GFP_1000_scores <- read.table("~/cluster_mount/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/codon_usage/by_expression/optimised/mCherry_1000_scores.tsv", quote="\"", comment.char="")
library("ggplot2")
p <- ggplot(data=GFP_1000_scores, aes(GFP_1000_scores[,1])) + geom_histogram(bins = '50')
p <- p + xlim(0.3, 1.0)
p <- p + xlab("Total codon score") + ylab("Count")
p <- p + xlab("Total codon score") + ylab("Count")

p <- p + geom_vline(xintercept = 0.98, na.rm = FALSE, show.legend = NA)
p
```
