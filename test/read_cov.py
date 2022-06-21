import pandas as pd
cov = pd.read_table('coverage_test.txt', header=0)
cov_info={}
##########cov_info is a dict storing coverage information###
#存储覆盖率信息
cov = cov.values[0: , [0,2]]
for i in range(cov.shape[0]):
    cov_info[cov[i , 0].split(' ')[0]] = cov[i , 1]
print(cov_info)