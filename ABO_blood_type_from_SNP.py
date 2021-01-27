# This code predict ABO blood type based on genotype .ped file
# The rules are: (/data100t1/home/wanying/shared_data_files/ABO_blood_type_usable_version.txt)
#       variant_alleles                     Blood_type
#       rs8176746(C:C),rs574347(T:T)        A
#       rs8176746(C:C),rs574347(T:C)	    O
#       rs8176746(A:C),rs574347(T:N)	    B

out_dir = '/data100t1/home/wanying/shared_data_files/'
out_fn ='ABO_out.txt'
out_fh = open(out_dir+out_fn, 'w')
with open(out_dir+'ABO_SNPs.ped', 'r') as fh:
    line = fh.readline().strip()
    while line != '':
        lst_tmp = line.split()
        individual_ID = lst_tmp[0]
        rs8176746 =lst_tmp[6 ] +lst_tmp[7]
        rs574347 =lst_tmp[8 ] +lst_tmp[9]
        blood_type ='NA'
        if rs8176746 =='AC' or rs8176746 =='CA':
            if rs574347 =='TC' or rs574347 =='CT': blood_type ='B'
        elif rs8176746 =='CC':
            if rs574347 == 'TT': blood_type ='A'
            elif rs574347 =='TC' or rs574347 == 'CT': blood_type ='O'
        out_fh.write(individual_ID + '\t' + blood_type + '\n')
        line = fh.readline().strip()
