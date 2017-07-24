#!/usr/bin/python

'''
This script uses text files of upregulated genes to create a count table for which genes are upregulated at which timepoint
'''

from sets import Set
import sys,argparse
from collections import defaultdict
import re
import numpy
import csv

#-----------------------------------------------------
# Step 1
# Import variables, load input files & create set of genes
# If using a different number of files, arguments & appending to list of genes will need to be changed
#-----------------------------------------------------

#These commands use the argparse module to import files specified in the command line
ap = argparse.ArgumentParser()
ap.add_argument('--input_1',required=True,type=str,help='text file of genes at 24hrs')
ap.add_argument('--input_2',required=True,type=str,help='text file of genes at 48hrs')
ap.add_argument('--out_file',required=True,type=str,help='the tsv file where the count table is output to')
conf = ap.parse_args()

#These commands create a list of the genes and loads the text files into a dictionary where

inp1_dict = defaultdict(list)
with open(conf.input_1) as f1:
    inp1_lines = f1.readlines()[1:]
    genes_list = []
    inp1 = []
    for x in inp1_lines:
        if ("NA" in x):
            continue
        else:
            genes_list.append(x.split('\t')[0])
            inp1.append(x.split('\t')[0])
            gene_name = x.split('\t')[0]
            value = float(x.split('\t')[2])
            inp1_dict[gene_name].append(value)

inp2_dict = defaultdict(list)
with open(conf.input_2) as f2:
    inp2_lines = f2.readlines()[1:]
    inp2 = []
    for x in inp2_lines:
        if ("NA" in x):
            continue
        else:
            genes_list.append(x.split('\t')[0])
            inp2.append(x.split('\t')[0])
            gene_name = x.split('\t')[0]
            value = float(x.split('\t')[2])
            inp2_dict[gene_name].append(value)

#This removes all duplicate entries from the list

genes = set(genes_list)

#-----------------------------------------------------
# Step 2
# Load gene names to a numpy array and create new columns
# If doing with a different number of files, change the number in the numpy.reshape() command
#-----------------------------------------------------

#This creates a numpy array, the headings are arbitrary, this does not have to be a timecourse, but the original script was

a = numpy.array(["Gene_Name", "24hr", "48hr"])

#These commands test if a gene is present in the DEG file or not, then if it has an absolute log2 fold change greater than the threshold of 1
#It then adds the name of the gene and the result of the test (1 or 0) to a numpy array

for x in genes:
    to_add = []
    to_add.append(x)
    try:
        b = inp1.index(x)
    except ValueError:
        to_add.append('0')
    else:
        for y in inp1_dict[x]:
            test = abs(y)
            if test > 1:
                to_add.append('1')
            else:
                to_add.append('0')
    try:
        c = inp2.index(x)
    except ValueError:
        to_add.append('0')
    else:
        for y in inp2_dict[x]:
            test = abs(y)
            if test > 1:
                to_add.append('1')
            else:
                to_add.append('0')
    a = numpy.append(a, to_add, axis=0)

#These commands reshape the array to ensure it writes to a file correctly

z = len(genes) + 1
a = numpy.reshape(a, (z, 3))

#These commands write out the results generated above as .tsv specified above

outfile = str(conf.out_file)
with open(outfile, "w") as o:
    writer = csv.writer(o, delimiter='\t')
    writer.writerows(a)
