#!/usr/bin/env python
from params import *

largs = [('i',workf+'/simyields','','input_folder'),
         ('name','','','file_prefix'),
         ('nf',300,'','number_of_files'),
         ('r',True,'','remove_files'),
         ('m',30,'','number_of_merging')]

[indir,name,nf,remove,merge] = ArgParser(largs).argl

error = []
for i in range(nf/merge):
  x = os.system('hadd -f '+indir+'/'+name+'_m_'+str(i)+'.root '+' '.join([indir+'/'+name+'_'+str(j)+'.root' for j in range(merge*i,merge*(i+1))]))
  error.append(x)

if remove and error==[0 for i in range(nf/merge)]:
  for i in range(nf/merge):
    if error[i]!=0: continue
    for j in range(merge):
      os.remove(indir+'/'+name+'_'+str(merge*i+j)+'.root')

  for i in range(nf/merge): os.rename(indir+'/'+name+'_m_'+str(i)+'.root',indir+'/'+name+'_'+str(i)+'.root')
      
      
