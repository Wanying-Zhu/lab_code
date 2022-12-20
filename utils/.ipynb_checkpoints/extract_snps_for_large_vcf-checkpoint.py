# This code takes in a list of variants, with chromosome nubmer and positions,
# then search in the gzipped (.gz) genotype data files and extract these snps
# Results are outputted in a file with user supplied file name
# This code does not use pandas so should be more appropriate for very large dataset (like BioVU)
import gzip
def find_variants(lst_pos, output_fn, input_fn, compression='gzip', input_col_name='POS'):
    '''
    This function takes in positions of a list of snp (no duplication),
    then checks the input genotype file and output found SNP into a output file
    Parameters:
        - lst_pos: positions (int) of snps. Must be list-like for iteration
        - output_fn: path and file name to output found SNPs
        - input_fn: path and file name of input genotype file. Must be a single chromosome
                    Assume all variants are sorted base on position in the input file
        - compression='gzip': compression type. Most genotype files are in .gz format, so use '.gzip' as default
        - input_col_name='POS': column name of position in input genotype file. Usually 'POS', could be 'current_pos'
    Returns:
        - Output found SNPs and genotyeps into output file.
        - Print out any missing SNPs into console
    '''
    lst_pos = list(lst_pos) # Convert lst_pos to a list, so that .sort() and .remove() will work
    lst_pos.sort() # Sort positions to be found
    input_fh = gzip.open(input_fn, 'rt')
    output_fh = open(output_fn, 'w')
    not_found_fh = open(output_fn+'.not_found_pos', 'w') # Write a position into this file if not found in genotype file
    line = input_fh.readline().strip()

    while line[0:2] == '##': # Read trhough header lines
        line = input_fh.readline().strip()

    output_fh.write(line) # Now reach the column header line, write it into output file
    output_fh.write('\n')

    line = input_fh.readline().strip()
    count = 0
    while len(lst_pos) != 0 and line != '':
        current_pos = int(line.split(maxsplit=2)[1]) # get position of current row

        if current_pos<=lst_pos[0]:
            if current_pos==lst_pos[0]: # If found a variant, write into output file
                output_fh.write(line)
                output_fh.write('\n')
                lst_pos.remove(current_pos) # Update list of positions

            line = input_fh.readline().strip() # Put code below here, since if current_pos > lst_pos[0], does not need to read a new line
            count += 1
            if count % 1000 == 0:
                print('.', end='', flush=True)
            if count % 50000 == 0:
                print(count, 'rows processed, now at position', current_pos)

        elif current_pos > lst_pos[0]: # If this variant is not in the genotype file, skip it and look for the next one
            not_found_fh.write(str(lst_pos[0])+'\n')
            lst_pos.remove(lst_pos[0])

    input_fh.close()
    output_fh.close()
    not_found_fh.close()


'''chr_number='chr21'
indiv_args=([16012318, 16675887, 16781136, 25840043, 26743463, 26914271, 26932893, 27038787, 27692049, 28737592, 28902698, 30093693, 31280183, 31696649, 32272654, 32449607, 32971387, 33556049, 33618410, 33862276, 34022513, 34188414, 34342245, 35995469, 36136716, 36317422, 36378531, 36703978, 36804031, 36964368, 37107585, 38293610, 38299554, 38325459, 38388839, 38433581, 38605152, 38635780, 39161847, 39267524, 39354887, 40056747, 42652963, 42656243, 42894021, 43105391, 43455249, 44805750, 44851537, 45240040, 46003595, 46129715, 46129801, 46387318],
            '/data100t1/home/wanying/GIANT_GLGC/BioVU_PRS_and_pheWAS/output/chr21_selected_SNPs.dose.vcf', '/data100t1/share/BioVU/TOPMed_imputed/EUR/merged_imputation_files/chr21_merged.dose.vcf.gz')


print(indiv_args[1]+'\n'+indiv_args[2])

from multiprocessing import Pool
with Pool(30) as p:
    p.starmap(find_variants, [(indiv_args[0],indiv_args[1],indiv_args[2])])'''
