#!/usr/bin/env python

from misc import *

def add_arg(d,k,v):
  if isinstance(d,dict):
    if not isistance(v,list):
      if k in d: d[k]+=v
      else: d[k]=v
    else:
      d = [d.copy() for x in v]
      argin = k in d[0]
      for x,y in zip(d,v): 
        if argin: x[k]+=y
        else: x[k]=y
  elif isinstance(d,list): d = [add_arg(x,k,v) for x in d]
  else: 
    sys.stdout.write('type not recognized')
    return []
  return d

def replace_args(d,x,y):
  if isinstance(d,dict):
    lk,lv = [],[]
    for k,v in d.items(): 
      if x in v:
        lk.append(k)
        lv.append([v.replace(x,e) for e in y])
        del d[k]
    d = [d.copy() for v in y]
    for k,v in zip(lk,lv):
      argin = k in d[0]
      for i,v2 in enumerate(v):
        if argin: d[i][k]+=v2
        else: d[i][k]=v2
    return d
  elif isinstance(d,list):
    d = [replace_args(z,x,y) for z in d]
  else:
    sys.stdout.write('type not recognized')
    return []
  return d

args = {'project':'prad','track':'analysis','memory':'1 GB'}

for i,x in enumerate(sys.argv):
  if x=='-list':
    args = replace_args(args,sys.argv[i+1],sys.argv[i+2].split(','))
  elif x[0]=='-': 
    if ',' in sys.argv[i+1]: v = sys.argv[i+1].split(',')
    else: v = sys.argv[i+1]
    args = add_arg(args,x[1:],v)
 
#print args 
jsub(**args)
