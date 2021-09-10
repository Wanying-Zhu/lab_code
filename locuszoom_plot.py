# Create a python 2.7 environment, run this code in that environment with:
#   python3 code_name.py
import os
from multiprocessing import Pool

input_fn = '/vgipiper05/hannah/stuttering/population_analysis/replication/AddHealth_sugen_b37_summary_rsids.txt'
marker_SNP_fn = '/vgipiper05/hannah/stuttering/population_analysis/replication/AddHealth_sentinelvars.txt'
locuszoom = '/data100t1/home/wanying/lab_code/LocusZoom/locuszoom/bin/locuszoom'
# input_fn = '/vgipiper05/hannah/stuttering/population_analysis/replication/STUT_saige_summary_rsids.txt'
# marker_SNP_fn = '/vgipiper05/hannah/stuttering/population_analysis/replication/Clinical_sentinelvars.txt'

# Get a list of marker names
lst_markers = []
with open(marker_SNP_fn) as fh:
    # LD R2 will be calcualted based on this marker SNP
    marker_variant = fh.readline().strip()
    while marker_variant != '':
        lst_markers.append(marker_variant)
        marker_variant = fh.readline().strip()


# marker_variant: LD R2 will be calcualted based on this marker SNP
def call_locuszoom(marker_variant):
    pop = 'EUR'
    build = 'hg19'
    source = '1000G_Nov2014'
    flank = '500kb'
    marker_col = 'rsid' # Column name of marker
    pvalcol = 'p' # Column nmae of p values
    command = locuszoom + ' --metal ' + input_fn + \
              ' --refsnp \"' + marker_variant + '"' +\
              ' --pop ' + pop +\
              ' --build '+ build +\
              ' --source ' + source +\
              ' --flank ' + flank +\
              ' --markercol ' + marker_col +\
              ' --pvalcol ' + pvalcol
    os.system('echo '+command)
    os.system(command)
        
# Use paralel programming
# multiprocessing.cpu_count() to check number of cores in the current system
with Pool(50) as p: # Max number of parallel processes on our server is 64
    p.map(call_locuszoom, lst_markers)
