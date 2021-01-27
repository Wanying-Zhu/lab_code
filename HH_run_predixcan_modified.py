import os
import sys
# import multiprocessing as mp

# jti_eqtl_db.txt contains file names of all file names of eqtls in different tissues
# dbfile=open('/data100t1/home/wanying/cchc/doc/select_individuals_for_RNAseq/jti_eqtl_db.txt','r')
# dbfile=open('/data100t1/home/wanying/cchc/doc/select_individuals_for_RNAseq/mashr_eqtl_db.txt','r')

# cpu = 11

# dblist = []
# for line in dbfile:
#     line = line.strip()
#     dblist.append(line)

# db is the model to use
# chr_number is chromosome number
def predix(db):
    db_tissue = db.split('/')[-1]
    # Use lift over is the vcf files are not in the same build as models, add the below line
    # Actual location of the chain file might be different
    # + ' --liftover /data100t1/home/wanying/lab_code/predixcan_test_run/data/hg19ToHg38.over.chain.gz' \
    command = 'python /data100t1/gapps/MetaXcan/software/Predict.py' \
              + ' --model_db_path ' + db \
              + ' --model_db_snp_key varID' \
              + ' --vcf_genotypes /vgipiper04/CCHC/TOPMed_postimpute_042020/chr*.dose.vcf.gz' \
              + ' --vcf_mode genotyped' \
              + ' --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38"'\
              + ' --prediction_output ./predixcan_mashr_output/' + db_tissue[:-3]+'.txt' \
              + ' --prediction_summary_output ./predixcan_mashr_output/' + db_tissue[:-3]+'_summary.txt' \
              + ' --verbosity 9 --throw'
    print(command)
    os.system(command)
#   os.system('gzip '+DB+'.out')
#     print('finish {}'.format(DB))

db = '/data100t1/share/predixcan_models/gtex_v8/mashr/eqtl/mashr/mashr_Whole_Blood.db'

# pool = mp.Pool(cpu)
#
# for db in dblist:
#     # Use this to run all tissues
#     if db == '/data100t1/share/predixcan_models/gtex_v8/mashr/eqtl/mashr/mashr_Whole_Blood.db':
#         predix(db)
#     pool.apply_async(predix, args=(db,))
#     predix(db)
#
# pool.close()
# pool.join()


"""
# Actual shell command with MASHR model
python /data100t1/gapps/MetaXcan/software/Predict.py \
--model_db_path /data100t1/share/predixcan_models/gtex_v8/mashr/eqtl/mashr/mashr_Whole_Blood.db \
-model_db_snp_key varID \
--vcf_genotypes /vgipiper04/CCHC/TOPMed_postimpute_042020/chr*.dose.vcf.gz \
--vcf_mode genotyped \
--on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
--prediction_output ./20210104_predixcan_mashr_output_ALL/mashr_Whole_Blood.txt \
--prediction_summary_output ./20210104_predixcan_mashr_output_ALL/mashr_Whole_Blood_summary.txt \
--verbosity 9 \
--throw


# To use JTI model
python /data100t1/gapps/MetaXcan/software/Predict.py \
--model_db_path /data100t1/share/predixcan_models/gtex_v8/JTI/JTI_Whole_Blood.db \
--vcf_genotypes /vgipiper04/CCHC/TOPMed_postimpute_042020/chr1.dose.vcf.gz \
--vcf_mode genotyped \
--prediction_output ./predixcan_jti_output/jti_Whole_Blood_chr1.txt \
--prediction_summary_output ./predixcan_jti_output/jti_Whole_Blood_chr1_summary.txt \
--verbosity 9 \
--throw
"""
