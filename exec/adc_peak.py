#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('o',workf,'','output_folder'),
         ('i',replayfm,'','input_folder'),
         ('b',False,'','batch_mode'),
         ('s',0,'','progress_bar')]
print_help(largs)
[lrun,outdir,indir,batch,show] = catch_args(largs)

if batch:
  runs = get_run_between(lrun[0],lrun[-1])
  args = del_args(largs,['r','b','s'])
  input_files = add_args_list(args,'r',runs,workf+'/jobs/adc_peak')
  jsub(project='prad',track='analysis',jobname='adcpeak_{0}_{1}'.format(lrun[0],lrun[-1]),command=execf+'/adc_peak.py args.in',input_files=input_files,input_data='args.in',memory='1024 MB',os='centos7')
  sys.exit()

from readdst import *
from misc_root import *

load_names()

t = chao(run=lrun[0],folder=indir)
t.goto(20000)
t.lload = ['','adc']

fout = TFile('{0}/adc_peak_{1}.root'.format(outdir,lrun[0]),'recreate')
hadc = listofhisto('hadc',module_names.values(),8192,0,8192)

for ev in progress(t,modulo=10000,show=show):
  if t.isepics: continue
  if event_is_bad(ev.iev,t.run): continue
  if ev.trigger not in [1,2]: continue
  for i,v in zip(ev.adc_id,ev.adc_val): 
    if 0<=i<1728: hadc[i].Fill(v)
  
writelisto(hadc,fout)
fout.Close()
