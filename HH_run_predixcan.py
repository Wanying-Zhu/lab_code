import os
import sys
import multiprocessing as mp

dbfile=open('mashr_eqtl_db.txt','r')
cpu = 11

dblist = []
for line in dbfile:
    line = line.strip()
    dblist.append(line)


def predix(db):
    DB = db.split('/')[-1]
    file_format = sys.argv[1]
#   print('working on {} ...'.format(DB))
#   os.system('~/tools/PrediXcan/Software/predict_gene_expression.py --dosages /data100t1/share/BioVU/GReX_predixcan/dosage/annotate/ --dosages_prefix chr --weights '+db+' --output ./'+DB+'.out')
    os.system('python /data100t1/gapps/MetaXcan/software/Predict.py --model_db_path '+db+' --model_db_snp_key varID --on_the_fly_mapping METADATA '+file_format+' --skip_palindromic --vcf_genotypes /vgipiper05/hannah/stuttering/population_analysis/merged/potential_ctrls/selected_controls_taketwo_09092020/postTOPMed/chr*.dose.clean.vcf.gz --vcf_mode imputed --prediction_output '+DB+' --prediction_summary_output '+DB+'_summary.txt --verbosity 9 --throw')
#   os.system('gzip '+DB+'.out')
    print('finish {}'.format(DB))

pool = mp.Pool(cpu)

for db in dblist:
    pool.apply_async(predix, args=(db,))

pool.close()
pool.join()

