#!/usr/bin/env python
from params import *
from struct import *

types = {'c':str,'b':int,'?':bool,'h':int,'i':int,'l':int,'q':int,'f':float,'d':float,'s':str,'p':str}

largs = [('f','','','format'),
         ('i','','','input file'),
         ('o','','','output file'),
         ('p',0,'several','float precision'),
         ('b',False,'','convert to binary')]
         
aparser = ArgParser(largs)
[fformat,infile,outfile,precision,binary] = aparser.argl


count = ''
fformat2 = []
for x in fformat:
  if 48<=ord(x)<=57: count+=x
  else:
    if x=='p' or x=='s': fformat2.append(types[x])
    else: fformat2.extend([types[x]]*int(count))
    count = ''

if binary:
  l = readfile(infile,fformat2)
  fout = open(outfile,'wb')
  for x in l:
    y = pack(fformat,*x)
    fout.write(y)
  fout.close()

else:
  size = calcsize(fformat)
  f = open(infile,'rb')
  l = []
  while f.read(size): 
    f.seek(-size,1)
    l.append(list(unpack(fformat,f.read(size))))
  writelist(outfile,l,ljust=min(10,max(precision)+2),ro=precision)
