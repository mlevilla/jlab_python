#!/usr/bin/env python

from params import *

largs = [('p','','','pass_number'),
         ('m','fprod','','clustering_method'),
         ('n','','','extension_name'),
         ('g','','','gains_input_file'),
         ('w',False,'','write_output_gains'),
         ('d','all','','lead_glass or crystal'),
         ('pe',['',''],'','period of production')]
print_help(largs)
[ipass,method,name,gainfile_in,write_gain,det,period] = [catch_arg(x,y,z) for x,y,z,c in largs]

from misc_root import *
from selection import *

load_names()
load_positions()
indir = work+'/calib_results/{0}/d{1}/'.format(method,name)
f = TFile(indir+'histo_'+period[0]+'_'+period[1]+'.root')
hr_ep = {x:f.Get(y).Get('hr_ep_'+y) for x,y in module_names.items()}
hr_moller = {x:f.Get(y).Get('hr_moller_'+y) for x,y in module_names.items()}
he_ep = {x:f.Get(y).Get('he_ep_'+y) for x,y in module_names.items()}
he_moller = {x:f.Get(y).Get('he_moller_'+y) for x,y in module_names.items()}
ebeam = 1100.44 if int(period[0])<1350 else 2142.0

# all modules

lratio_ep, lsigma_ep, lratio_moller, lsigma_moller, lalpha, ldalpha = [],[],[],[],[],[]
lmean_ep, lmean_moller = [],[]
linteg_ep, linteg_moller, lexp_ep, lexp_moller, = [],[],[],[]
lerror = []

for i,x in module_names.items():
  if ((det=='pwo' and x[0]=='G') or (det=='lg' and x[0]=='W')): continue
  hr_ep[i].GetXaxis().SetRangeUser(0.5,1.5)
  hr_moller[i].GetXaxis().SetRangeUser(0.2,1.5)
  if hr_ep[i].Integral()<10 or hr_moller[i].Integral()<10:
    print x,'no enough stat'
    lmean_ep.append(-1); lsigma_ep.append(0); 
    lmean_moller.append(-1); lsigma_moller.append(0);
    lratio_ep.append(-1); lratio_moller.append(-1)
    lalpha.append(-1)
    linteg_ep.append(0); linteg_moller.append(0)
    lexp_ep.append(0); lexp_moller.append(0)
    lerror.append([x,'stat'])
    continue

  theta = module_theta[i]*deg
  exp_ep_energy = Ef_el(theta,ebeam/1000.)*1000 
  exp_ep_energy -= energy_loss(theta,exp_ep_energy)
  exp_moller_energy = 1000*moller_energy(theta,ebeam/1000.)
  exp_moller_energy -= energy_loss(theta,exp_moller_energy)

  if i in get_regions('inner2'):    
    r_ep = hr_ep[i].Fit('gaus','rwS0Q','',0.96,1.2)
  else:
    r_ep = hr_ep[i].Fit('gaus','rwS0Q','',0.8,1.2)
  r_moller = hr_moller[i].Fit('gaus','rS0Q','',0.5,1.2)
  
  integ_ep = hr_ep[i].Integral()
  integ_moller = hr_moller[i].Integral()
  ratio_ep = r_ep.GetParams()[1]
  ratio_moller = r_moller.GetParams()[1]
  dratio_ep = r_ep.GetErrors()[1]
  dratio_moller = r_moller.GetErrors()[1]
  sratio_ep = r_ep.GetParams()[2]
  sratio_moller = r_moller.GetParams()[2]

  if i>1000: e_ep = he_ep[i].Fit('gaus','rwS0Q','',exp_ep_energy-50,exp_ep_energy+50)
  else: e_ep = he_ep[i].Fit('gaus','rwS0Q','',exp_ep_energy-100,exp_ep_energy+100)
  e_moller = he_moller[i].Fit('gaus','rwS0Q','',exp_moller_energy-50,exp_moller_energy+50)
  
  mean_ep = e_ep.GetParams()[1]
  mean_moller = e_moller.GetParams()[1]

  error = []
  if abs(ratio_ep-1)>0.05: error.append('mean')
  if abs(sratio_ep)>0.15: error.append('sigma')
  if error!=[]: lerror.append([x,' '.join(error)])

  lratio_ep.append(round(ratio_ep,5))
  lratio_moller.append(round(ratio_moller,5))
  lsigma_ep.append(round(sratio_ep,5))
  lsigma_moller.append(round(sratio_moller,5))
  lalpha.append(round(1000*(mean_moller/exp_moller_energy-1.)/(mean_moller-exp_ep_energy),5))
  linteg_ep.append(int(integ_ep))
  linteg_moller.append(int(integ_moller))
  lexp_ep.append(round(exp_ep_energy,2))
  lexp_moller.append(round(exp_moller_energy,2))
  lmean_ep.append(round(mean_ep,2))
  lmean_moller.append(round(mean_moller,2))
  
names = [x for x in module_names.values() if (det=='all') or (det=='lg' and x[0]=='G') or (det=='pwo' and x[0]=='W')]
writelist(indir+'ep_'+period[0]+'_'+period[1]+'_'+det+name+'.txt',[[u,v,w,x,y,z] for u,v,w,x,y,z in zip(names,lratio_ep,lsigma_ep,linteg_ep,lmean_ep,lexp_ep)],ljust=9)
writelist(indir+'moller_'+period[0]+'_'+period[1]+'_'+det+name+'.txt',[[u,v,w,x,y,z] for u,v,w,x,y,z in zip(names,lratio_moller,lsigma_moller,linteg_moller,lmean_moller,lexp_moller)],ljust=9)
writelist(indir+'alpha_'+period[0]+'_'+period[1]+'_'+det+name+'.txt',[[u,v,w,x,y,z] for u,v,w,x,y,z in zip(names,lalpha,lmean_ep,lmean_moller,lexp_ep,lexp_moller)],ljust=9)

writelist(indir+'error_'+period[0]+'_'+period[1]+'_'+det+name+'.txt',lerror,ljust=9)

if gainfile_in!='':
  gains0 = readfile(gainfile_in,[str]+[float]*6)
  names0 = [x[0] for x in gains0]
  gains1 = [gains0[names0.index(x)] for x in module_names.values()]
  gains2 = [[w[0],round(w[1]/x,6),y,z,w[4],w[5],w[6]] if abs(x-1)<0.05 else w for w,x,y,z in zip(gains1,lratio_ep,lexp_ep,lalpha)]
  writelist(indir+'cal_'+period[0]+'_'+period[1]+'.txt',gains2,ljust=10)
