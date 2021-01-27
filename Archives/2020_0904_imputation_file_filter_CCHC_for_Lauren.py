# This code is to filter CCHC and SOL files for Lauren
# Filtering criteria are (apply four types and get four results):
#   1. MAF<0.01, duplicates, multi-allelic, INDEL, info <0.5
#   2. MAF<0.01, duplicates, multi-allelic, INDEL, info <0.8
#   3. MAF<0.01, duplicates, multi-allelic, INDEL, info <0.5, A/T, C/G
#   4. MAF<0.01, duplicates, multi-allelic, INDEL, info <0.8, A/T, C/G

import gzip
import pandas as pd

# ----------------------- Helper functions -----------------------

# This function get SNP list of cchc files
# Return number of SNPs under each filter condition
def cchc_get_snp_list(cchc_input_dir = '/vgipiper04/CCHC/1KG_postimpute_03302020/', chromosome = 'chr21',
                      output_dir = '/data100t1/home/wanying/Lauren/snp_list/', maf_threshold = 0.01):
    cchc_vcf_fn = chromosome + '.info.gz'
    df = pd.read_csv(cchc_input_dir + cchc_vcf_fn, sep='\t', compression='gzip', dtype={'AvgCall':str, 'Rsq':str})
    df = df[['SNP', 'REF(0)', 'ALT(1)', 'MAF', 'Rsq']]
    df['Rsq'] = pd.to_numeric(df['Rsq'], errors='coerce')
    df.dropna(inplace=True)
    df.drop_duplicates(subset='SNP', inplace=True)

    indl_mask_1 = (df['REF(0)'] == 'A') | (df['REF(0)'] == 'T') | (df['REF(0)'] == 'G') | (df['REF(0)'] == 'C')
    indl_mask_2 = (df['ALT(1)'] == 'A') | (df['ALT(1)'] == 'T') | (df['ALT(1)'] == 'G') | (df['ALT(1)'] == 'C')
    atcg_mask = ((df['REF(0)'] == 'A') & (df['ALT(1)'] == 'T')) | ((df['REF(0)'] == 'T') & (df['ALT(1)'] == 'A')) | \
                ((df['REF(0)'] == 'C') & (df['ALT(1)'] == 'G')) | ((df['REF(0)'] == 'G') & (df['ALT(1)'] == 'C'))
    maf_mask = (df['MAF'] >= maf_threshold)
    info_mask_05 = df['Rsq'] >= 0.5
    info_mask_08 = df['Rsq'] >= 0.8

    # Return number of variants
    number_of_snp_filter_1 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_05].shape[0]
    number_of_snp_filter_2 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_08].shape[0]
    number_of_snp_filter_3 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_05 & atcg_mask].shape[0]
    number_of_snp_filter_4 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_08 & atcg_mask].shape[0]

    return number_of_snp_filter_1, number_of_snp_filter_2, number_of_snp_filter_3, number_of_snp_filter_4

"""
    # Two files will be generated:
    #   - Remove duplicates, MAF<0.01, indels
    #   - Remove duplicates, MAF<0.01, indels, A/T C/G
    (df.loc[mask_no_atcg])['SNP'].to_csv(output_dir + 'cchc_' + chromosome + '_maf_duplct_indl.txt', sep='\t',
                                         index=False, header=True)
    (df.loc[mask_with_atcg])['SNP'].to_csv(output_dir + 'cchc_' + chromosome + '_maf_duplct_indl_ATCG.txt', sep='\t',
                                           index=False, header=True)
"""



"""
# Process CCHC variant counts
output_dir = '/data100t1/home/wanying/Lauren'
output_fn = 'cchc_ouput.txt'
fh = open(output_dir+output_fn, 'w')
fh.write('chromosome'+'\t'+'count_filter1'+'\t'+'count_filter2'+'\t'+'count_filter3'+'\t'+'count_filter1'+'\n')
for i in range(1, 23):
    chromosome = 'chr' + str(i)
    num1, num2, num3, num4 = cchc_get_snp_list(chromosome=chromosome)
    output_line = chromosome+'\t'+str(num1)+'\t'+str(num2)+'\t'+str(num3)+'\t'+str(num4)+'\n'
    fh.write(output_line)
    print(output_line, flush=True)"""




# Get allele and minor allele info from SOL_imputed3_chr*.vcf.gz file
# (.metrics.gz does not have info of reference and alternative allele info for some variants)
# Return a dataframe with 3 columns:
#   - SNP id, reference allele, alternative allele
def get_ref_alt_allele(sol_input_dir = '/data100t1/share/SOL/imputed/', chromosome = 'chr21'):
    sol_input_fn = 'SOL_imputed3_' + chromosome + '.vcf.gz'
    fh = gzip.open(sol_input_dir+sol_input_fn, 'rt')
    snp_id = []
    ref_allele = []
    alt_allele = []
    for i in range(5):  # Skip header lines
        line = fh.readline().strip()

    print('Getting ref and alt allele information from .vcf.gz file of', chromosome)   # Keep console busy
    count = 1
    while line != '':
        lst_tmp = line.split()
        snp_id.append(lst_tmp[2])
        ref_allele.append(lst_tmp[3])
        alt_allele.append(lst_tmp[4])
        line = fh.readline().strip()

        # Keep console busy
        count = count + 1
        if count%10000 == 0:
            print()
            print(count, 'SNPs have been processed')
        elif count%5000 == 0:
            print('.\n')
        elif count%200 == 0:
            print('.', end='', flush=True)
    fh.close()
    df = pd.DataFrame()
    df['SNP_id'] = snp_id
    df['REF'] = ref_allele
    df['alt'] = alt_allele
    print('Done')
    return df

# Process SOL variant counts
# Merge with dataframe from get_ref_alt_allele()
#   - Already saved that dataframe to simplify the whole process
#   - allele_df is a dataframe of none-duplicated SNPs with allele information
#     (read from saved dataframe file )
def sol_get_snp_list(allele_df, sol_input_dir = '/data100t1/share/SOL/imputed/', chromosome = 'chr21',
                     output_dir = '/data100t1/home/wanying/Lauren/', maf_threshold = 0.01):
    sol_input_fn = 'SOL_imputed3_' + chromosome + '.metrics.gz'
    df = pd.read_csv(sol_input_dir + sol_input_fn, sep=' ', compression='gzip', dtype={'exp_freq_a1':str})
    df = df[['rs_id','info', 'exp_freq_a1']]
    df['exp_freq_a1'] = pd.to_numeric(df['exp_freq_a1'], errors='coerce')
    df = df.merge(allele_df, left_on='rs_id', right_on='SNP_id', how='inner')
    df.dropna(inplace=True)
    df.drop_duplicates(subset='rs_id', inplace=True)
    df = df[['rs_id', 'REF', 'alt', 'exp_freq_a1', 'info']]

    indl_mask_1 = (df['REF'] == 'A') | (df['REF'] == 'T') | (df['REF'] == 'G') | (df['REF'] == 'C')
    indl_mask_2 = (df['alt'] == 'A') | (df['alt'] == 'T') | (df['alt'] == 'G') | (df['alt'] == 'C')
    atcg_mask = ((df['REF'] == 'A') & (df['alt'] == 'T')) | ((df['REF'] == 'T') & (df['alt'] == 'A')) | \
                ((df['REF'] == 'C') & (df['alt'] == 'G')) | ((df['REF'] == 'G') & (df['alt'] == 'C'))
    maf_mask = (df['exp_freq_a1'] >= maf_threshold) & (df['exp_freq_a1'] <= (1-maf_threshold))
    info_mask_05 = df['info'] >= 0.5
    info_mask_08 = df['info'] >= 0.8

    # Return number of variants
    number_of_snp_filter_1 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_05].shape[0]
    number_of_snp_filter_2 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_08].shape[0]
    number_of_snp_filter_3 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_05 & atcg_mask].shape[0]
    number_of_snp_filter_4 = df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_08 & atcg_mask].shape[0]
    # df.loc[maf_mask & indl_mask_1 & indl_mask_2 & info_mask_05].to_csv(output_dir+'sol_test.txt')
    return number_of_snp_filter_1, number_of_snp_filter_2, number_of_snp_filter_3, number_of_snp_filter_4



# sol_get_snp_list()
# output_dir = '/data100t1/home/wanying/Lauren/snp_list/'

output_dir = '/data100t1/home/wanying/Lauren/snp_list/'
output_fn = 'sol_ouput.txt'
fh = open(output_dir+output_fn, 'w')
fh.write('chromosome'+'\t'+'count_filter1'+'\t'+'count_filter2'+'\t'+'count_filter3'+'\t'+'count_filter1'+'\n')
for i in range(1, 23):
    chromosome = 'chr' + str(i)
    sol_fn = 'sol_' + chromosome + '.txt'
    allele_df = pd.read_csv(output_dir + sol_fn, sep='\t')
    # Drop duplicates
    allele_df.drop_duplicates(subset='SNP_id', inplace=True)
    num1, num2, num3, num4 = sol_get_snp_list(chromosome=chromosome, allele_df=allele_df)
    output_line = chromosome+'\t'+str(num1)+'\t'+str(num2)+'\t'+str(num3)+'\t'+str(num4)+'\n'
    fh.write(output_line)
    print(output_line, flush=True)
