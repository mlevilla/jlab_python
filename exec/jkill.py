#!/usr/bin/env python

import os,sys,re
from misc import *
from arg_parser import *
from subprocess import *

largs = [('n',0,'several','delete_between_id1_id2'),
         ('s','','','delete_match')]
[n,s] = ArgParser(largs).argl

if n==[0] and s=='':
  os.system('jkill 0')
  sys.exit()

output = check_output(['jobstat','-u','mlevilla'])
l = [[y for y in x.split(' ') if y!=''] for x in output.split('\n')]
if n[0]<=0: n[0] = int(l[1][0])
if len(n)==1: n.append(int(l[-2][0]))

ids = []
for x in l[1:-1]:
  if re.match(s,x[5]) and n[0]<=int(x[0])<=n[1]: ids.append(x[0])

if len(ids)!=0: os.system('jkill '+' '.join(ids))

