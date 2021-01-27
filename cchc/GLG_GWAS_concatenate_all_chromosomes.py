#! /data100t1/gapps/anaconda3/bin/python
# Not finish!!!!

# This code concatnate chr1 to 22 data together to make kinship matrix
# The zgrep code in the analysis plan stopped working somehow

import gzip
fn = '/vgipiper04/CCHC/1KG_postimpute_03302020/chr2.dose.vcf.gz'
fh = gzip.open(fn)
line = fh.readline().strip()
while line != '':
    print(line)
    break

