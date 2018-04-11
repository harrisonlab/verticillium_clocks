#!/usr/bin/python

'''python
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

'''
Run with commands:

for GeneGff in $(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33_parsed.gff3); do
    Strain=JR2
    Organism=V.dahliae
    Assembly=$(ls public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
    InterPro=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
    SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa)
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_53)
    Dir2=$(ls -d RNA_alignment/featureCounts/experiment_53/WT53)
    DEG_Files=$(ls \
        $Dir1/Frq53_LD_06h.txt \
        $Dir1/Wc153_LD_06h.txt \
        $Dir1/Wt53_Frq53_bl06h.txt \
        $Dir1/Wt53_Frq53_d06h.txt \
        $Dir1/Wt53_bl06h_vs_d06h.txt \
        $Dir1/Wt53_Wc153_bl06h.txt \
        $Dir1/Wt53_Wc153_d06h.txt \
        $Dir2/Wt53_d06h_d12h.txt \
        $Dir2/Wt53_d06h_d18h.txt \
        $Dir2/Wt53_d06h_d24h.txt \
        $Dir2/Wt53_d12h_d18h.txt \
        $Dir2/Wt53_d12h_d24h.txt \
        $Dir2/Wt53_d24h_d18h.txt \
        $Dir2/Wt53_LD_06h.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts_53.txt)
    FPKM=$(ls $Dir1/countData_53.fpkm)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_incl_exp.tsv
done

'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--gene_gff',required=True,type=str,help='Gff file of predicyted gene models')
ap.add_argument('--gene_fasta',required=True,type=str,help='amino acid sequence of predicted proteins')
ap.add_argument('--DEG_files',required=True,nargs='+',type=str,help='space spererated list of files containing DEG information')
ap.add_argument('--raw_counts',required=True,type=str,help='raw count data as output from DESeq')
ap.add_argument('--fpkm',required=True,type=str,help='normalised fpkm count data as output from DESeq')
ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
conf = ap.parse_args()



with open(conf.gene_gff) as f:
    gene_lines = f.readlines()

with open(conf.gene_fasta) as f:
    prot_lines = f.readlines()

DEG_files = conf.DEG_files
DEG_dict = defaultdict(list)
for DEG_file in DEG_files:
    with open(DEG_file) as f:
        filename = DEG_file
        DEG_lines = f.readlines()
        for line in DEG_lines:
            if line.startswith('baseMean'):
                continue
            else:
                split_line = line.split()
                gene_name = split_line[0]
                log_change = split_line[2]
                P_val = split_line[6]
                entryname = "_".join([filename, gene_name])
                DEG_dict[entryname].extend([log_change, P_val])

with open(conf.raw_counts) as f:
    raw_count_lines = f.readlines()

with open(conf.fpkm) as f:
    fpkm_lines = f.readlines()

with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()

with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()


#-----------------------------------------------------
# Load protein sequence data into a dictionary
#-----------------------------------------------------

prot_dict = defaultdict(list)
for line in prot_lines:
    line = line.rstrip()
    if line.startswith('>'):
        header = line.split(' ')[0]
        header = header.replace('>', '')
    else:
        prot_dict[header] += line


#-----------------------------------------------------
#
# Build a dictionary of raw count data
#
#-----------------------------------------------------

raw_read_count_dict = defaultdict(list)

line1 = raw_count_lines.pop(0)
line1 = line1.rstrip("\n")
count_treatment_list = line1.split("\t")
count_treatment_list = list(filter(None, count_treatment_list))

for line in raw_count_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    transcript_id = split_line.pop(0)
    for i, treatment in enumerate(count_treatment_list):
        # i = i-1
        raw_read_count = float(split_line[i])
        dict_key = "_".join([transcript_id, treatment])
        raw_read_count_dict[dict_key].append(raw_read_count)

#-----------------------------------------------------
#
# Build a dictionary of normalised fpkm data
#
#-----------------------------------------------------

fpkm_dict = defaultdict(list)

line1 = fpkm_lines.pop(0)
line1 = line1.rstrip("\n")
fpkm_treatment_list = line1.split("\t")
fpkm_treatment_list = list(filter(None, fpkm_treatment_list))

for line in fpkm_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    transcript_id = split_line.pop(0)
    for i, treatment in enumerate(fpkm_treatment_list[1:]):
        if 'NA' in split_line[i]:
            continue
        fpkm = float(split_line[i])
        dict_key = "_".join([transcript_id, treatment])
        fpkm_dict[dict_key].append(fpkm)

#-----------------------------------------------------
#
# Build a dictionary of interproscan annotations
# Annotations first need to be filtered to remove
# redundancy. This is done by first loading anntoations
# into a set.
#-----------------------------------------------------

interpro_set =  Set([])
interpro_dict = defaultdict(list)

for line in InterPro_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    interpro_columns = []
    index_list = [0, 4, 5, 11, 12]
    for x in index_list:
        if len(split_line) > x:
            interpro_columns.append(split_line[x])
    set_line = ";".join(interpro_columns)
    if set_line not in interpro_set:
        gene_id = interpro_columns[0]
        gene_id = gene_id.replace('.p', '.t')
        interpro_feat = ";".join(interpro_columns[1:])
        interpro_dict[gene_id].append(interpro_feat)
    interpro_set.add(set_line)


#-----------------------------------------------------
#
# Build a dictionary of Swissprot annotations
#-----------------------------------------------------

swissprot_dict = defaultdict(list)

for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    gene_id = gene_id.replace('.p', '.t')
    swissprot_columns = itemgetter(14, 12, 13)(split_line)

    swissprot_dict[gene_id].extend(swissprot_columns)

#-----------------------------------------------------
# Step 3
# Itterate through genes in file, identifying if
# they ahve associated information
#-----------------------------------------------------

# Print header line:
header_line = ['transcript_id']
header_line.extend(['contig', 'start', 'stop', 'strand'])
for treatment in set(count_treatment_list):
    treatment = "raw_count_" + treatment
    header_line.append(treatment)

for treatment in set(fpkm_treatment_list):
    treatment = "fpkm_" + treatment
    header_line.append(treatment)

for DEG_file in DEG_files:
    file_name = DEG_file.split('/')[-1]
    header_line.append("LFC_" + file_name)
    header_line.append("P-val_" + file_name)
header_line.append('CDS_seq')
header_line.append('Swissprot')
header_line.append('Interpro')
print ("\t".join(header_line))

transcript_lines = []

for line in gene_lines:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    split_line = line.split()
    if 'transcript' in split_line[2] or 'mRNA' in split_line[2]:
        transcript_lines.append("\t".join(split_line))

for line in transcript_lines:
    split_line = line.split("\t")
    useful_cols = [split_line[0],  split_line[3], split_line[4], split_line[6]]
    prot_seq = ''
    swissprot_cols = []
    interpro_col = []
    # Identify gene id
    if 'ID' in split_line[8]:
        split_col9 = split_line[8].split(';')
        transcript_id = "".join([ x for x in split_col9 if 'ID' in x ])
        transcript_id = transcript_id.replace('ID=', '').replace('transcript:', '')
    else:
        transcript_id = split_line[8]
    gene_id = transcript_id.split(".")[0]
    DEG_out = []
    for DEG_file in DEG_files:
        entryname = "_".join([DEG_file, gene_id])
        if DEG_dict[entryname]:
            DEG_out.append(DEG_dict[entryname][0])
            DEG_out.append(DEG_dict[entryname][1])
        else:
            DEG_out.append('.')
            DEG_out.append('.')

    # # Add in read count data:
    mean_count_cols = []
    for treatment in set(count_treatment_list):
        dict_key = "_".join([gene_id, treatment])
        expression_values = raw_read_count_dict[dict_key]
        # print expression_values
        mean_count = np.mean(expression_values)
        mean_count = np.round_(mean_count, decimals=0)
        mean_count_cols.append(mean_count.astype(str))
    # print mean_count_cols
    mean_fpkm_cols = []
    for treatment in set(fpkm_treatment_list):
        dict_key = "_".join([gene_id, treatment])
        # print dict_key
        expression_values = fpkm_dict[dict_key]
        # print expression_values
        mean_fpkm = np.mean(expression_values)
        # print mean_fpkm
        mean_fpkm = np.round_(mean_fpkm, decimals=0)
        mean_fpkm_cols.append(mean_fpkm.astype(str))

    # # Add in Swissprot info
    if swissprot_dict[transcript_id]:
        swissprot_cols = swissprot_dict[transcript_id]
    else:
        swissprot_cols = ['.','.','.']
    # Add in interproscan info
    if interpro_dict[transcript_id]:
        interpro_col = "|".join(interpro_dict[transcript_id])
    else:
        interpro_col = '.'

    prot_seq = "".join(prot_dict[transcript_id])
    outline = [transcript_id]
    outline.extend(useful_cols)
    outline.extend(mean_count_cols)
    outline.extend(mean_fpkm_cols)
    outline.extend(DEG_out)
    outline.append(prot_seq)
    outline.extend(swissprot_cols)
    outline.append(interpro_col)


    print "\t".join(outline)
