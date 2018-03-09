#!/usr/bin/env python
from params import *
from moddensity_utils import *
from ROOT import TF2

largs = [('i','','','input_file'),
         ('o','.','','output_folder'),
         ('n','','','suffix')]

aparser = ArgParser(largs)
[infile,outdir,suffix] = aparser.argl

from misc_root import *
f = TFile(infile) 
h = getlisto(f.dist)
hE = getlisto(f.Edist)
correct = True if f.Edist_correct else False
if correct: hEc = getlisto(f.Edist_correct)
nbin = h[0].GetNbinsX()

if suffix=='':
  i1 = infile.index('ecorrect2D_')
  i2 = infile.index('.root')
  suffix = infile[i1+10:i2]

fout = TFile(outdir+'/ecorrect2D_fit'+suffix+'.root','recreate')

region_names = ['pwo','lg_top','lg_right','lg_bottom','lg_left','lg']
h_regions = listofhisto('h',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-1.,1.,nbin,-1.,-1.,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;')
hE_regions = listofhisto('hE',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-1.,1.,nbin,-1.,-1.,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;E\'/E_{theo}')

if correct: hEc_regions = listofhisto('hEc',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-1.,1.,nbin,-1.,-1,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;E\'/E_{theo}')

ffit = TF2('ffit','[0]*(1+[1]*x^2+[2]*y^2+[3]*x^2*y^2+[4]*x^4+[5]*y^4)',-0.5,0.5,-0.5,0.5)
ffit.SetParameters(1,0,0,0,0,0)
ffit.SetParLimits(0,0.9,1.5)

lparam = [[] for i in range(706)]
for i in range(len(h)):
  name = h[i].GetName()
  i1 = name.index('_')+1
  i2 = name.rindex('_')
  index = int(name[i1:i2])
  ep = int(name[i2+1:])
  if index<16: sector=1
  elif index<32: sector=2
  elif index<48: sector=3
  elif index<64: sector=4
  else: sector=0
  h_regions[sector*2+ep].Add(h[i])
  hE_regions[sector*2+ep].Add(hE[i]) 
  if correct: hEc_regions[sector*2+ep].Add(hEc[i])
  if sector>0:
    h_regions[10+ep].Add(h[i])
    hE_regions[10+ep].Add(hE[i]) 
    if correct: hEc_regions[10+ep].Add(hEc[i])

lparam = [[] for i in range(706)]
for i in range(len(h)):
  hE[i].Divide(h[i])
  if correct: hEc[i].Divide(h[i])
  if h[i].GetEntries()>200:
    hE[i].Fit(ffit,'q')
    lparam[i/2+353*(i%2)] = [str(i/2)+['_ep','_ee'][i%2]]+[hE[i].GetFunction('ffit').GetParameter(k) for k in range(6)]
  else:
    lparam[i/2+353*(i%2)] = [str(i/2)+['_ep','_ee'][i%2]]+[1,0,0,0,0,0]

lparam_mod_ep, lparam_mod_ee = [],[]
for x,y in module_names.items():
  idx = module2index_asym(x)
  lparam_mod_ep.append([y]+lparam[idx][1:])
  lparam_mod_ee.append([y]+lparam[idx+353][1:])

lparam_regions = [[] for i in range(12)]  
# regions
ffit.ReleaseParameter(0)
ffit.SetParameters(1,0,0,0,0,0)
for i in range(12):
  hE_regions[i].Divide(h_regions[i])
  if correct: hEc_regions[i].Divide(h_regions[i])
  hE_regions[i].Fit(ffit,'q')
  lparam_regions[i/2+6*(i%2)] =  [region_names[i/2]+['_ep','_ee'][i%2]]+[hE_regions[i].GetFunction('ffit').GetParameter(k) for k in range(6)]

writelisto(hE,fout,'Edist',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.9,1.1])
writelisto(hE_regions,fout,'regions',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.95,1.05])
if correct:
  writelisto(hEc,fout,'Edist_correct',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.9,1.1])
  writelisto(hEc_regions,fout,'regions_correct',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.95,1.05])
writelist(outdir+'/ecorrect2D_params'+suffix+'.txt',lparam,ljust=10,ro=5)
writelist(outdir+'/ecorrect2D_params'+suffix+'_ep.txt',lparam_mod_ep,ljust=10,ro=5)
writelist(outdir+'/ecorrect2D_params'+suffix+'_ee.txt',lparam_mod_ee,ljust=10,ro=5)
writelist(outdir+'/ecorrect2D_params_regions'+suffix+'.txt',lparam_regions,ljust=[15]+[10]*6,ro=5)
fout.Close()
f.Close()
