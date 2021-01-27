# 2020/06/11
# This script read in a .fam file and select random subject ID
# The purpose is to subset original (hla) data to do pca
# Then project rest of the samples on to this pca model

import sys
import numpy as np

try:
    args = sys.argv  # Read file name from commandline
    input_fn = args[1]
    output_fn = args[2]
except:
    print('No input file name provided, try again using following format:')
    print('  python random_list.py [input_file_name] [output_file_name]')
    quit()

# /data100t1/home/wanying/hla_analysis/output/MEGA_94_HLA_QCed_filtered_het_hwe_maf.fam
sub_set_percentage = 0.1  # take 10% of individuals from original file
count = 0  # count number of individuals in the original file
individual_ID_lst = []  # empty list to store individual IDs
with open(input_fn, 'r') as fh:
    line = fh.readline().strip()
    while line != '':
        count += 1
        # Get individual ID
        # plink1.9 --keep reads in a list with family IDs in the first column
        # and within-family IDs in the second column
        individual_ID = line.split()[0] + '\t' + line.split()[1]
        individual_ID_lst.append(individual_ID)
        line = fh.readline().strip()

number_of_individuals_in_subset = int(count * sub_set_percentage)
random_line_index_lst = []  # An empty list to store random line index
for i in range(number_of_individuals_in_subset):
    rand_line_index = np.random.randint(count)
    while rand_line_index in random_line_index_lst:    # Don't add repeated line index
        rand_line_index = np.random.randint(count)
    random_line_index_lst.append(rand_line_index)

subset_individual_IDs = np.array(individual_ID_lst)[tuple(random_line_index_lst), ]
np.savetxt(output_fn, subset_individual_IDs, fmt='%s')

print('Total number of individuals:', count)
print('Subset size (% of total):', sub_set_percentage * 100)
print('Output file is:', output_fn)
