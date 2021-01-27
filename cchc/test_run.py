# import pandas as pd

# fn = 'wanying_cchc.csv'
fn = 'wanying_cchc.txt'
c_line = 0
c_field = 0
col_lst = []
line5_lst = []
line2_lst = []
line3_lst = []
line22_lst = [] # use it as line 38
with open(fn, 'r') as fh:
    line = fh.readline().strip()
    while line != '':
        c_line += 1
        # tmp =  line.split(',')
        tmp = line.split('*')
        for i in tmp:
            c_field += 1
            if c_line == 1:
                col_lst.append(i)
            elif c_line == 5:
                line5_lst.append(i)
            elif c_line == 38:
                line22_lst.append(i)
            # elif c_line == 3:
            #     line3_lst.append(i)
            # elif c_line == 22:
            #     line22_lst.append(i)

        #
        # if c_line == 5:
        #     print('\n================')
        #     print(c_field)
        c_field = 0
        line = fh.readline().strip()

# print(len(col_lst), len(line_lst))
for i in range(len(col_lst)):
    if i <479:
        # print('col #', i, ':', col_lst[i], '|', line2_lst[i], '|',
        #       line3_lst[i], '|', line5_lst[i], '|', line22_lst[i])
        print('col #', i, ':', col_lst[i], '|', line5_lst[i], '|', line22_lst[i])

print('\n================')
print('\nLast item in line 5:',line5_lst[-1])
print('\nLast 5 item in line 38:',line22_lst[-5:])
# print('\nLen of column:',len(col_lst))
# print('\nLen of line 2:',len(line2_lst))
# print('\nLen of line 5:',len(line3_lst))
# print('\nLen of line 5:',len(line5_lst))
# print('\nLen of line 22:',len(line22_lst))





