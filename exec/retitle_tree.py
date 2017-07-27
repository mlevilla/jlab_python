#!/usr/bin/env python

from params import *
from ROOT import TFile,TTree

largs = [('i',work+'/tree/island','','input_directory'),
         ('o',work,'','output_directory'),
         ('m','island','','method'),
         ('b',False,'','batch_mode')]

[indir,outdir,method,batch] = [catch_arg(x,y,z) for x,y,z,_ in largs]

if batch:
  jsub(track='analysis',project='prad',command=os.path.abspath('retitle_tree.py'),options='-i {0} -o {1} -m {2}'.format(indir,outdir,method),jobname='retitle_tree_'+method,memory='4 GB')
  sys.exit()

l = os.listdir(indir)
l.sort()
for x in l:
  if not (('tree' in x) and ('.root' in x)): continue 
  print x
  idx1 = x.rindex('_')
  idx2 = x.index('.root')
  run = x[idx1+1:idx2]
  f = TFile(indir+'/'+x)
  fout = TFile(outdir+'/tree_'+method+'_'+run+'.root','recreate')
  t = f.event
  t2 = t.CloneTree()
  t2.SetTitle(method+'_'+run)
  fout.cd()
  t2.Write()
  fout.Close()
  f.Close()

sys.exit()
