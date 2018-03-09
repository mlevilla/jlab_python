#!/usr/bin/env python
from params import *
from moddensity_utils import *

largs = [('i','','','input_file'),
         ('o','.','','output_folder'),
         ('n','','','suffix')]

aparser = ArgParser(largs)
[infile,outdir,suffix] = aparser.argl

from misc_root import *
f = TFile(infile) 
h = getlisto(f.dist)
hE = getlisto(f.Edist)
nbin = h[0].GetNbinsX()

if suffix=='':
  i1 = infile.index('ecorrect_')
  i2 = infile.index('.root')
  suffix = infile[i1+8:i2]

fout = TFile(outdir+'/ecorrect_fit'+suffix+'.root','recreate')

ffit = TF1('ffit','[0]*(1+[1]*x^2+[2]*x^4)',-0.5,0.5)
ffit.SetParameters(1.,0.04,0.3)
ffit.SetParLimits(0,0.9,1.5)
ffit.SetParLimits(2,-1,10.)

lparam = [[] for i in range(3456)]
for i in range(len(h)):
  h[i].Rebin(5)
  hE[i].Rebin(5)
  hE[i].Divide(h[i])
  if h[i].GetEntries()>200:
    hE[i].Fit(ffit,'q','',-0.5,0.5)

for i in range(len(h)/2):
  if hE[i].GetFunction('ffit') and hE[i+3456].GetFunction('ffit'):
    meanx = hE[i].GetFunction('ffit').GetParameter(0)
    meany = hE[i+3456].GetFunction('ffit').GetParameter(0)
    chi2x = hE[i].GetFunction('ffit').GetChisquare()
    chi2y = hE[i+3456].GetFunction('ffit').GetChisquare()
    mean = (1/chi2x*meanx+1/chi2y*meany)/(1/chi2x+1/chi2y)
    ffit.FixParameter(0,mean)
    hE[i].Fit(ffit,'q','',-0.5,0.5)
    hE[i+3456].Fit(ffit,'q','',-0.5,0.5)
    lparam[i/2+1728*(i%2)] =  [hE[i].GetName()[4:-2]+['_ep','_ee'][i%2]]+[mean]+[hE[i].GetFunction('ffit').GetParameter(k) for k in range(1,3)]+ [hE[i+3456].GetFunction('ffit').GetParameter(k) for k in range(1,3)]
  else:
    lparam[i/2+1728*(i%2)] = [hE[i].GetName()[4:-2]+['_ep','_ee'][i%2]]+[1,0,0,0,0]

writelisto(hE,fout,'Edist',yrg=[0.9,1.1])
writelist(outdir+'/ecorrect_params_mod'+suffix+'.txt',lparam,ljust=10,ro=5)
fout.Close()
f.Close()
