#!/usr/bin/python
import requests
import time
import sys
from datetime import date

t0 = time.time()
dir = 'db/ldlink'
#f = 'SNP_5e-08_129.tsv'
f = sys.argv[1]
print('Load SNP file: %s' %f)
snp = open(f,'r')
snpid = snp.read().split('\n')
i = 1
for id in snpid[0:len(snpid)-1]: # To remove empty file generation
    print('%d/%d = %s' % (i,len(snpid)-1,id))
    pop  = '%2B'.join(('CEU','TSI','FIN','GBR','IBS'))
    r2_d = 'd'
    url1 = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?'
    url2 = ''.join(('var=',id,'&pop=',pop,'&r2_d=',r2_d))
    url  = ''.join((url1,url2))
    
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
    
<<<<<<< HEAD:T1D_ldlink_dn.py
    if(divmod(i,10)[1]==0):
        t1 = time.time()
        m,s = divmod(t1-t0,60)
        h,m = divmod(m,60)
        print(time.strftime('> Job time= %02d:%02d:%02d\n' % (h,m,s)))
=======
>>>>>>> parent of 0c366a1... write: db/roadmap - core functions:T1D_ldlink.py
    if(i==len(snpid)): print("\n>> Download process completed.")
    i+=1
t1 = time.time()
m,s = divmod(t1-t0,60)
h,m = divmod(m,60)
print(time.strftime('>> Job time= %02d:%02d:%02d\n' % (h,m,s)))