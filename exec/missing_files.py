#!/usr/bin/env python

from misc import *
from arg_parser import *

largs = [('p','.*','','pattern'),
         ('n',10,'','number_of_files'),
         ('e','.','','extension'),
         ('i','.','','input_folder'),
         ('x','','','exclude_pattern'),
         ('s',0.,'','percentage_expected_size'),
         ('t',0,'','expected_time')]

[pattern,n,extension,indir,exclude,size,tim] = ArgParser(largs).argl

print 'match pattern',pattern, 'in', indir
if exclude!='': print 'exclude pattern',exclude

l = [x for x in os.listdir(indir) if re.match(pattern,x) and (exclude=='' or (not re.match(exclude,x)))] 

if size!=0.:
  filesize = [os.stat(indir+'/'+x).st_size for x in l]
  mean = median(filesize)
  std = (sum((x-mean)**2 for x in filesize)/len(filesize))**0.5
  if mean>2**30: mean2,std2,unit = mean/2**20, int(std)/2**20, 'Mo'
  elif mean>2**20: mean2,std2,unit = round(float(mean)/2**20,1), round(std/2**20,1), 'Mo'
  elif mean>2**10: mean2,std2,unit = mean/2**10, int(std)/2**10, 'ko'
  else: mean2,std2,unit = mean, int(std), 'o'
  print 'average size of',mean2,unit,'with standard deviation of',std2,unit
  l = [x for i,x in enumerate(l) if abs(filesize[i]-mean)/float(mean)<size] 

if tim!=0:
  filetime = [os.stat(indir+'/'+x).st_mtime for x in l]
  mean = median(filetime)
  std = (sum((x-mean)**2 for x in filetime)/len(filetime))**0.5
  print 'created in average around',time.ctime(mean),'with a standard deviation of',parse_time(std)
  l = [x for i,x in enumerate(l) if filetime[i]-mean>-tim*3600] 

if len(l)==n:
  print 'all files are present'
  sys.exit()

l2 = []
for x in l:
  j = x.rindex(extension)
  i = x[:j].rindex('_')
  l2.append(int(x[i+1:j]))
  
l3 = list(set(range(n)) - set(l2))
l3.sort()
print ' '.join([str(x) for x in l3])

sys.exit()
