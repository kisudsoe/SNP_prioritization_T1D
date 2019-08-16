#!/usr/bin/python
# This file is for downloading LDlink data from risk SNP list.

# System parameter
import requests
import time
import sys
import os
import csv
import pandas as pd

## Command Arg Parameters ##
# Usage: Python ldlink_dn.py [SNP_list_file_path]
path = sys.argv[1]
print('Load SNP file: %s' % path)

t0 = time.time()
dir = 'db/ldlink'
#f = 'data/gwas_5e-08_129.tsv'
os.mkdir(dir)

#########################
## Function start here ##
#########################
snpdf = pd.read_csv(path,sep='\t')
snpid = list(set(snpdf['rsid'])) # unique IDs
print('SNP rsid number = %d' % len(snpid))

i = 1
for id in snpid[0:len(snpid)]: # To remove empty file generation
    print('%d/%d = %s' % (i,len(snpid),id))
    pop  = '%2B'.join(('CEU','TSI','FIN','GBR','IBS'))
    r2_d = 'd'
    token= '&token=669e9dc0b428'
    url_base = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?'
    url_options = ''.join(('var=',id,'&pop=',pop,'&r2_d=',r2_d,token))
    url  = ''.join((url_base,url_options))
    
    response = requests.get(url=url)
    print('  status_code = %d' % response.status_code) # Returned '200' means normal status.
    db = response.text
    db2 = db.splitlines()
    print('  line number = %d' % len(db2))
    path_w = ''.join((dir,'/',id,'.tsv'))
    with open(path_w,'w') as f:
        for row in db:
            f.write(row)
    print('  file saved = %s\n' % path_w)
    if(divmod(i,10)[1]==0):
        t1 = time.time()
        m,s = divmod(t1-t0,60)
        h,m = divmod(m,60)
        print(time.strftime('> Job time= %02d:%02d:%02d\n' % (h,m,s)))
    if(i==len(snpid)): print(">> Download process completed.")
    i+=1
t1 = time.time()
m,s = divmod(t1-t0,60)
h,m = divmod(m,60)
print(time.strftime('>> Job time= %02d:%02d:%02d\n' % (h,m,s)))
##################
## Function end ##
##################