#!/usr/bin/python

'''
This script generates a codon optimised protein based upon a fasta protein
sequence and a table of relative codon usage.
'''

from sets import Set
import sys,argparse
from collections import defaultdict
import re
import numpy as np
import csv
import random
from Bio import SeqIO

#-----------------------------------------------------
# Step 1
# Import variables, load input files & create set of genes
# If using a different number of files, arguments & appending to list of genes will need to be changed
#-----------------------------------------------------

#These commands use the argparse module to import files specified in the command line
ap = argparse.ArgumentParser()
ap.add_argument('--fasta_aa',required=True,type=str,help='protein sequence for conversion')
ap.add_argument('--fasta_cds',required=True,type=str,help='cds for conversion')
ap.add_argument('--codon_table',required=True,type=str,help='text file containing codon usage table')
ap.add_argument('--prefix',required=True,type=str,help='output directory/filename prefix for output files')
conf = ap.parse_args()


#-----------------------------------------------------
# Step 1
# Import variables, load input files & create set of genes
# If using a different number of files, arguments & appending to list of genes will need to be changed
#-----------------------------------------------------

class AA_weight_obj(object):
    """

    """
    def __init__(self, aa):
        """ """
        self.aa = aa
        self.weightings = defaultdict(float)
        self.weightings_adj = defaultdict(float)
        self.max = float()
        self.optimal = str()
        self.codons = []
        self.sorted_adj_weightings = []
        self.sorted_codons = []
        self.weight_list = []
        self.weight_list_adj = []
    def add_weight(self, codon, weight):
        """ """
        # print codon
        # print weight
        self.weightings[codon] = float(weight)
        # if float(weight) > self.max:
        #     self.max = float(weight)
        #     self.optimal = codon
        self.codons.append(codon)
        self.weight_list.append(weight)
    def random_codon(self):
        """ """
        num_codons = len(self.codons)
        r = float(random.randrange(0,10000, 1))
        # r = float(random.randrange(0,num_codons*100, 1))
        # print (self.aa)
        # print(r)
        r = np.divide(r, 10000)
        # r = np.divide(r, 100)
        # print(" of max ".join([str(r), str(num_codons)]))
        for x,y in zip(self.codons,self.sorted_adj_weightings):
            # print(" - ".join([str(r), str(x), str(y)]))
            selected_codon = x
            if float(y) >= float(r):
                break
            else:
                r = r - float(y)
        return selected_codon
    def get_opt(self):
        """ """
        # sorted_weightings = sorted(self.weight_list)
        # sorted_codons = [x for _,x in sorted(zip(self.weight_list,self.codons))]
        # print sorted_weightings
        # print sorted_codons
        # return sorted_codons[-1]
        return self.sorted_codons[-1]
    def adjust_weight(self):
        """ """
        num_codons = len(self.weight_list)
        # print num_codons
        # print(self.weight_list)
        self.weight_list_adj = [round(np.divide(float(x), num_codons),5) for x in self.weight_list]
        # print self.weight_list_adj
        self.sorted_adj_weightings = sorted(self.weight_list_adj)
        self.sorted_codons = [x for _,x in sorted(zip(self.weight_list_adj,self.codons))]
        for x,y in zip(self.sorted_codons, self.sorted_adj_weightings):
            self.weightings_adj[x] = y
        self.max = self.sorted_adj_weightings[-1]


class CodonTab_obj(object):
    """

    """
    def __init__(self):
        """Return a Expression_obj whose name is *gene_id*"""
        # self.organism = []
        self.weighting_dict = defaultdict(list)
        # self.codon_obj_dict = {}
        self.codon_dict = {
        'UUU':'F','UUC':'F',
        'UUA':'L','UUG':'L','CUU':'L','CUC':'L','CUA':'L','CUG':'L',
        'AUU':'I','AUC':'I','AUA':'I',
        'AUG':'M',
        'GUU':'V', 'GUC':'V','GUA':'V','GUG':'V',
        'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
        'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
        'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
        'UAU':'Y','UAC':'Y',
        'UAA':'X','UAG':'X',
        'CAU':'H','CAC':'H',
        'CAA':'Q','CAG':'Q',
        'AAU':'N','AAC':'N',
        'AAA':'K','AAG':'K',
        'GAU':'D','GAC':'D',
        'GAA':'E','GAG':'E',
        'UGU':'C','UGC':'C',
        'UGA':'X',
        'UGG':'W',
        'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGU':'S','AGC':'S',
        'AGA':'R','AGG':'R',
        'GGU':'G','GGC':'G', 'GGA':'G','GGG':'G'
        }
    def add_table(self, table):
        """"""
        table = table.replace(' ', '')
        table_lines = table.split(';')
        for line in table_lines:
            split_line = line.split(':')
            codon = split_line[0]
            # print codon
            weighting = split_line[1]
            # print weighting
            aa = self.codon_dict[codon]
            if self.weighting_dict[aa] and self.weighting_dict[aa][0]:
                obj = self.weighting_dict[aa][0]
                # print obj.weightings
            else:
                obj = AA_weight_obj(aa)
            obj.add_weight(codon, weighting)
            self.weighting_dict[aa].append(obj)
        for aa in self.weighting_dict.keys():
            self.weighting_dict[aa][0].adjust_weight()




def optimise_rand(prot):
    new_seq = ''
    for aa in prot:
        new_aa = vd_table_obj.weighting_dict[aa][0].random_codon()
        new_seq = new_seq + new_aa
    return(new_seq)

def optimise_best(prot):
    new_seq = ''
    for aa in prot:
        # print aa
        # new_aa = vd_table_obj.weighting_dict[aa][0].get_opt()
        new_aa = vd_table_obj.weighting_dict[aa][0].sorted_codons[-1]
        new_seq = new_seq + new_aa
    return(new_seq)

def optimise_worst(prot):
    new_seq = ''
    for aa in prot:
        # print aa
        new_aa = vd_table_obj.weighting_dict[aa][0].sorted_codons[0]
        new_seq = new_seq + new_aa
    return(new_seq)

def score_seq(seq, table_obj):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    total_score = float(0)
    total_max = float(0)
    for codon in codons:
        aa = table_obj.codon_dict[codon]
        score = table_obj.weighting_dict[aa][0].weightings_adj[codon]
        # score = score - table_obj.weighting_dict[aa][0].weight_list_adj[0]
        max = table_obj.weighting_dict[aa][0].max
        total_score = total_score + score
        total_max = total_max + max
    return [round(np.divide(total_score, total_max), 2), round(np.divide(total_max, total_max), 2)]
    # scores = []
    # for aa in seq.split(''):
    #     scores.append(score_dict[aa])

#-----------------------------------------------------
# Step X
#
#-----------------------------------------------------

seq_records = list(SeqIO.parse(conf.fasta_aa, "fasta"))
cds_records = list(SeqIO.parse(conf.fasta_cds, "fasta"))
prefix = conf.prefix


with open(conf.codon_table) as f:
    table_lines = []
    for line in f.readlines():
        table_lines.append(line.rstrip())

#-----------------------------------------------------
# Step X
#
#-----------------------------------------------------


record = seq_records[0]
# print record
prot = record.seq


# prot = 'MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK'

table = "".join(table_lines)
# table = 'UUU: 0.55; UCU: 0.85; UAU: 0.40; UGU: 0.44; UUC: 1.45; UCC: 1.41; UAC: 1.60; UGC: 1.56; UUA: 0.07; UCA: 0.51; UAA: 1.04; UGA: 1.06; UUG: 0.55; UCG: 1.36; UAG: 0.90; UGG: 1.00; CUU: 0.84; CCU: 0.93; CAU: 0.50; CGU: 0.97; CUC: 2.49; CCC: 1.66; CAC: 1.50; CGC: 2.45; CUA: 0.23; CCA: 0.53; CAA: 0.50; CGA: 0.75; CUG: 1.81; CCG: 0.89; CAG: 1.50; CGG: 0.71; AUU: 0.95; ACU: 0.58; AAU: 0.37; AGU: 0.39; AUC: 1.91; ACC: 1.62; AAC: 1.63; AGC: 1.49; AUA: 0.14; ACA: 0.58; AAA: 0.26; AGA: 0.36; AUG: 1.00; ACG: 1.22; AAG: 1.74; AGG: 0.76; GUU: 0.73; GCU: 0.80; GAU: 0.61; GGU: 0.91; GUC: 2.20; GCC: 1.98; GAC: 1.39; GGC: 2.32; GUA: 0.18; GCA: 0.44; GAA: 0.48; GGA: 0.46; GUG: 0.88; GCG: 0.77; GAG: 1.52; GGG: 0.31'


vd_table_obj = CodonTab_obj()

vd_table_obj.add_table(table)

# for k in vd_table_obj.weighting_dict.keys():
#     print(vd_table_obj.weighting_dict[k][0].weightings)

# print(prot)

#-----------------------------------------------------
# Step X
# Optimise codons - random weightings
#-----------------------------------------------------

print("randomised codons:")
new_cds = optimise_rand(prot)
print(new_cds)
seq_score, max = score_seq(new_cds, vd_table_obj)
print(" of ".join([str(seq_score), str(max)]))



#-----------------------------------------------------
# Step X
# Optimise codons - optimum codons
#-----------------------------------------------------

print("optimum sequence:")
new_cds = optimise_best(prot)
print(new_cds)
seq_score, max = score_seq(new_cds, vd_table_obj)
print(" of ".join([str(seq_score), str(max)]))


#-----------------------------------------------------
# Step X
# Optimise codons - worst codons
#-----------------------------------------------------

print("worst sequence:")
new_cds = optimise_worst(prot)
print(new_cds)
seq_score, max = score_seq(new_cds, vd_table_obj)
print(" of ".join([str(seq_score), str(max)]))


#-----------------------------------------------------
# Step X
# Score 1000 sequences for optimisation scores
#-----------------------------------------------------

score_list = []
cds_list = []
f = open("_".join([prefix, "1000_seqs.fa"]), "w+")
for i in range(0, 10000, 1):
    new_cds = optimise_rand(prot)
    seq_score, max = score_seq(new_cds, vd_table_obj)
    # print seq_score
    cds_list.append(new_cds)
    score_list.append(str(round(seq_score, 2)))
    f.write(">cds_" + str(i) + "_" + str(seq_score))
    f.write(new_cds)
f.close()

f = open("_".join([prefix, "1000_scores.tsv"]), "w+")
f.write("\n".join(score_list))
f.close()

midpoint_score = sorted(score_list)[500]
sorted_cds = [x for _,x in sorted(zip(score_list,cds_list))]
midpoint_cds = sorted_cds[500]
print("midpoint sequence:")
print midpoint_score
print midpoint_cds


#-----------------------------------------------------
# Step X
# Score the pre-optimised sequence
#-----------------------------------------------------

print("Score of the pre-optimised sequence:")
for record in cds_records:
    print record.id
    old_cds = str(record.seq)
    old_cds = old_cds.replace('T', 'U')
    # print old_cds
    seq_score, max = score_seq(old_cds, vd_table_obj)
    print(" of ".join([str(seq_score), str(max)]))

# print(score_list)

# #set matplotlib to use a backend suitable for headless operation
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
#
# plt.hist(score_list, bins='auto')
# out='tmp.png'
# plt.savefig(out, dpi=300, bbox_inches='tight')

# rng = np.random.RandomState(10)  # deterministic random data
# a = np.hstack((rng.normal(size=1000),
#     rng.normal(loc=5, scale=2, size=1000)))
