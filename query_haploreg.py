'''
author: Wanying Zhu
date: 2019/10/23
python_version: 3
'''

import urllib
import pandas as pd

# This function queries haploreg and return a dataframe of results
# Position result in version 3 is hg19 (not sure why they do not have both)
# Position result in version 4 is hg38

# Parameters needed:
#  - version: query haploreg 4.1 by default. set it to 3 if need to query haloreg version 3.
#             Or change the url link as needed
#  - input_snps (query): input SNPs in a string, each SNP is separated by ','
#        SNP number per query is set to be 1000, otherwise it takes too long and haploreg may refuse to do it.
#        15000 SNPs did not work for test run.
#  - gwas_id: the dropdown list to choose GWAS study, when no file or query SNP(s) is provided
#  - r2_threshold (ldThresh): r^2 threshold, default is 0.8 in this code
#  - ldPop: 1000G Phase 1 population for LD calculation. Other options includes AFR, AMR and ASN.
#  - epi: Source for epigenomes
#  - cons: Mammalian conservation algorithm. 'siphy'=SiPhy-omega, 'gerp'=GERP, or 'both'
#  - genetypes: Show position relative to
#  - output: set output result type to 'text' for python code to process

# Return: everyting from Haploreg in a dataframe
def query_haploreg_all_results(input_snp,
                               version=4,
                               r2_threshold=0.8,
                               ldPop='EUR',
                               epi='vanilla',
                               cons='siphy',
                               genetypes='gencode'):

    params_library = {'query':input_snp,
                      'gwas_id':0,
                      'ldThresh': r2_threshold,
                      'ldPop': ldPop,
                      'epi': epi,
                      'cons': cons,
                      'genetypes': genetypes,
                      'output':'text'}
   
    # parameters passed to the website, needs to be parsed and in binary
    params = urllib.parse.urlencode(params_library).encode("utf-8")
    
    # url of HaploReg3
    if version == 3:
        print('Query haploreg version 3')
        url = 'https://pubs.broadinstitute.org/mammals/haploreg/haploreg_v3.php'
    # url of HaploReg4.1
    elif version == 4:
        print('Query haploreg version 4.1')
        url = 'https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php'
    else:
        print('Wrong version of haploreg.\nNo result returned')
        return('')
    
    # Query with parameters
    query = urllib.request.urlopen(url, params)
    
    content = query.read().decode("utf-8") # This is a massive string...
    content_list = content.rstrip().split('\n')  # First element contains the line of column title
    
    if version == 3:
        # 1st line in result of version 3 need to be removed
        content_list = content_list[1:]   
    
    content_df = pd.DataFrame(columns=content_list[0].split('\t'))
    
    for i in range(1, len(content_list)):
        content_df.loc[i] = content_list[i].split('\t')

    return(content_df)