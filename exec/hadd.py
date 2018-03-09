#!/usr/bin/env python

import os,re 
from arg_parser import *

largs = [('i','.','','input_dir'),
         ('o','','','output_dir'),
         ('s','','','suffix'),
         ('n','','','name'),
         ('r',-1,'several','range'),
         ('f',False,'','force')]

[indir,outdir,suffix,name,rg,force] = ArgParser(largs).argl

if not outdir: outdir = indir

def extract_number(s):
  i0 = s.index(name)+len(name)
  i2 = s.index(suffix+'.root',i0)
  i1 = s.index('_',i0)
  print s,i0,i1,i2
  if i1!=i2: return int(s[i1+1:i2])
  
l = [x for x in os.listdir(indir) if re.match(name+'_[0-9]*'+suffix+'\.root',x)]

i0 = 0
for i in range(len(l)):
  not_ok = False
  run = extract_number(l[i-i0])
  print run
  if not run:  not_ok = True
  elif len(rg)==2 and (run<rg[0] or run>rg[1]): not_ok = True
  elif len(rg)>2 and run not in rg: not_ok = True
  if not_ok:
    del l[i-i0]
    i0+=1
    continue

print l

os.system('hadd {0}{1}/{2}{3}.root {4}'.format('-f '*force,outdir,name,suffix,' '.join([indir+'/'+x for x in l])))

s = raw_input('delete input files (y/n):')
if s in ['y','yes','Y','Yes']:
  for x in l: os.remove(indir+'/'+x)
