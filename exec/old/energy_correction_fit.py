#!/usr/bin/env python
from params import *
from moddensity_utils import *

largs = [('i','','','input_file'),
         ('o','.','','output_folder'),
         ('n','','','suffix'),
         ('c',False,'','correction')]

aparser = ArgParser(largs)
[infile,outdir,suffix,correct] = aparser.argl

from misc_root import *
f = TFile(infile) 
h = getlisto(f.dist)
hE = getlisto(f.Edist)
if correct: hEc = getlisto(f.Edist_correct)
nbin = h[0].GetNbinsX()

if suffix=='':
  i1 = infile.index('ecorrect_')
  i2 = infile.index('.root')
  suffix = infile[i1+8:i2]

fout = TFile(outdir+'/ecorrect_fit'+suffix+'.root','recreate')

h_regions = listofhisto('h_regions',range(20),nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hE_regions = listofhisto('hE_regions',range(20),nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;E\'/E_{theo}')

if correct: hEc_regions = listofhisto('hEc_regions',range(20),nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;E\'/E_{theo}')

ffit = TF1('ffit','[0]*(1+[1]*x^2+[2]*x^4)',-0.5,0.5)
ffit.SetParameters(1.,0.04,0.3)
ffit.SetParLimits(0,0.9,1.5)
ffit.SetParLimits(2,-1,10.)

lparam = [[] for i in range(706)]
for i in range(len(h)):
  name = h[i].GetName()
  i1 = name.index('_')+1
  i2 = name.rindex('_')
  index = int(name[i1:i2])
  axis = 0 if 'hx' in name else 1
  ep = int(name[i2+1:])
  if index<16: sector=1
  elif index<32: sector=2
  elif index<48: sector=3
  elif index<64: sector=4
  else: sector=0
  h_regions[10*axis+sector*2+ep].Add(h[i])
  hE_regions[10*axis+sector*2+ep].Add(hE[i]) 
  if correct: hEc_regions[10*axis+sector*2+ep].Add(hEc[i])

lparam = [[] for i in range(706)]
for i in range(len(h)):
  hE[i].Divide(h[i])
  if correct: hEc[i].Divide(h[i])
  if h[i].GetEntries()>200:
    hE[i].Fit(ffit,'q','',-0.5,0.5)

for i in range(len(h)/2):
  if hE[i].GetFunction('ffit') and hE[i+706].GetFunction('ffit'):
    meanx = hE[i].GetFunction('ffit').GetParameter(0)
    meany = hE[i+706].GetFunction('ffit').GetParameter(0)
    chi2x = hE[i].GetFunction('ffit').GetChisquare()
    chi2y = hE[i+706].GetFunction('ffit').GetChisquare()
    mean = (1/chi2x*meanx+1/chi2y*meany)/(1/chi2x+1/chi2y)
    ffit.FixParameter(0,mean)
    hE[i].Fit(ffit,'q','',-0.5,0.5)
    hE[i+706].Fit(ffit,'q','',-0.5,0.5)
    lparam[i/2+353*(i%2)] =  [str(i/2)+['_ep','_ee'][i%2]]+[mean]+[hE[i].GetFunction('ffit').GetParameter(k) for k in range(1,3)]+ [hE[i+706].GetFunction('ffit').GetParameter(k) for k in range(1,3)]
  else:
    lparam[i/2+353*(i%2)] = [str(i/2)+['_ep','_ee'][i%2]]+[1,0,0,0,0]

lparam_mod_ep, lparam_mod_ee = [],[]
for x,y in module_names.items():
  idx = module2index_asym(x)
  lparam_mod_ep.append([y]+lparam[idx][1:])
  lparam_mod_ee.append([y]+lparam[idx+353][1:])

# regions
ffit.ReleaseParameter(0)
ffit.SetParameters(1.,0.04,0.3)
for i in range(20):
  hE_regions[i].Divide(h_regions[i])
  if correct: hEc_regions[i].Divide(h_regions[i])
  hE_regions[i].Fit(ffit,'q','',-0.5,0.5)
  
region_names = ['pwo','lg_top','lg_right','lg_bottom','lg_left']
lparam_regions = [[] for i in range(10)]  
for i in range(10):
  meanx = hE_regions[i].GetFunction('ffit').GetParameter(0)
  meany = hE_regions[i+10].GetFunction('ffit').GetParameter(0)
  chi2x = hE_regions[i].GetFunction('ffit').GetChisquare()
  chi2y = hE_regions[i+10].GetFunction('ffit').GetChisquare()
  mean = (1/chi2x*meanx+1/chi2y*meany)/(1/chi2x+1/chi2y)
  ffit.FixParameter(0,mean)
  hE_regions[i].Fit(ffit,'q','',-0.5,0.5)
  hE_regions[i+10].Fit(ffit,'q','',-0.5,0.5)
  lparam_regions[i/2+5*(i%2)] =  [region_names[i/2]+['_ep','_ee'][i%2]]+[mean]+[hE_regions[i].GetFunction('ffit').GetParameter(k) for k in range(1,3)]+ [hE_regions[i+10].GetFunction('ffit').GetParameter(k) for k in range(1,3)]

writelisto(hE,fout,'Edist',xrg=[-0.6,0.6],yrg=[0.9,1.1])
writelisto(hE_regions,fout,'regions',xrg=[-0.6,0.6],yrg=[0.95,1.05])
if correct:
  writelisto(hEc,fout,'Edist_correct',xrg=[-0.6,0.6],yrg=[0.9,1.1])
  writelisto(hEc_regions,fout,'regions_correct',xrg=[-0.6,0.6],yrg=[0.95,1.05])
writelist(outdir+'/ecorrect_params'+suffix+'.txt',lparam,ljust=10,ro=5)
writelist(outdir+'/ecorrect_params'+suffix+'_ep.txt',lparam_mod_ep,ljust=10,ro=5)
writelist(outdir+'/ecorrect_params'+suffix+'_ee.txt',lparam_mod_ee,ljust=10,ro=5)
writelist(outdir+'/ecorrect_params_regions'+suffix+'.txt',lparam_regions,ljust=[15]+[10]*5,ro=5)
fout.Close()
f.Close()
