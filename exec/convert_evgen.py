#!/usr/bin/env python
from misc import *
from struct import *

types = {'c':str,'b':int,'?':bool,'h':int,'i':int,'l':int,'q':int,'f':float,'d':float,'s':str,'p':str}

fformat = sys.argv[1]
infile = sys.argv[2]

count = ''
fformat2 = []
for x in fformat:
  if 48<=ord(x)<=57: count+=x
  else:
    if x=='p' or x=='s': fformat2.append(types[x])
    else: fformat2.extend([types[x]]*int(count))
    count = ''

if infile[-1]!='2':
  l = readfile(infile,fformat2)
  fout = open(infile+'2','wb')
  for x in l:
    y = pack(fformat,*x)
    fout.write(y)
  fout.close()

else:
  size = calcsize(fformat)
  f = open(infile)
  l = []
  while f.read(size): 
    f.seek(-size,1)
    l.append(list(unpack(fformat,f.read(size))))
  writelist(infile[:-1],l,ljust=10)
