#!/usr/bin/env python
import sys
import os
import math
import gzip
marker_dict = {}
input_file = 'PRS.var.04172021.38.txt'
FH = open(input_file)
for line in FH:
    line = line.rstrip()
    t = line.split('\t')
    if t[0] != 'CHR':
        loci = t[0]+'-'+t[1]
        marker_dict[loci]={"OR":float(t[2]),"Risk_allele":t[3]}
FH.close()

input_list = sys.argv[1]
sample_name = {}
sample_prs = {}
FH0 = open(input_list)
for line0 in FH0:
    line0 = line0.rstrip()
    input_file = line0
    FH = gzip.open(input_file,"rt")
    for line in FH:
        line = line.rstrip()
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                t = line.split('\t')
                for i in range(9,len(t)):
                    sample_name[i]=t[i]
        else:
            t = line[0:100].split('\t')
            if t[0].startswith('chr'):
                t[0] = t[0][3:]
            loci = t[0]+'-'+t[1]
            if loci in marker_dict:
                t = line.split('\t')
                risk_allele = marker_dict[loci]["Risk_allele"]
                gt_list = {}
                c=0
                gt_list[c]=t[3]
                for alt in t[4].split(','):
                    c=c+1
                    gt_list[c]=alt
                f_ind = {}
                for i,x in enumerate(t[8].split(':')):
                    f_ind[x]=i
                if "GT" in f_ind:
                    for i in range(9,len(t)):
                        x=t[i].split(':')
                        gt = x[f_ind['GT']]
                        if '/' in gt:
                            gt = gt.split('/')
                        if '|' in gt:
                            gt = gt.split('|')
                        allele_sum = sum([gt_list[int(xx)]==risk_allele for xx in gt])
                        if sample_name[i] not in sample_prs:
                            sample_prs[sample_name[i]] = [math.log(marker_dict[loci]["OR"]) * allele_sum]
                        else:
                            sample_prs[sample_name[i]].append(math.log(marker_dict[loci]["OR"]) * allele_sum)
    FH.close()
FH0.close()

for sample in sample_prs:
    print(sample+'\t'+str(sum(sample_prs[sample]))+'\t'+str(len(sample_prs[sample])))
