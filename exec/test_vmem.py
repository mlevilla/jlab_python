#!/usr/bin/env python

from params import *
from subprocess import *
from misc_root import *

def check_mem(message):
  output = check_output(['top','-u','mlevilla','-n','1'])
  l = [x.split(' ') for x in output.split('\n')]
  virt,res,shr = 0,0,0
  for i,x in enumerate(l):
    while '' in l[i]: l[i].remove('')
    if 'mlevilla' not in l[i]: continue
    j = l[i].index('mlevilla')
    virt+=int(l[i][j+3])
    res+=int(l[i][j+4])
    shr+=int(l[i][j+5])
  print message,virt,res,shr
    
check_mem('start')
# f = TFile(work+'/pradsim-1/output/moller_elastic_11.root')
f = TFile(work+'/trees/island/tree_island_1288.root')
fout = TFile('test.root','recreate')
check_mem('file opening')
# t = f.T
t = f.event
n = t.GetEntries()
h = listofhisto('h',range(100),1000,-10,10)

check_mem('histo made')

for i,_ in enumerate(t):
  if i%(n/10)==0: check_mem('event '+str(i))
  h[randint(0,90)].Fill(randint(-10,10))

writelisto(h,fout)
check_mem('end')
fout.Close()
f.Close()
