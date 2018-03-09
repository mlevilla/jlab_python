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
h = getlisto(f.raw_correct) if 'correct' in infile else getlisto(f.raw)
nbin = h[0].GetNbinsX()

if suffix=='':
  i1 = infile.index('density_')
  i2 = infile.index('.root')
  suffix = infile[i1+7:i2]

fout = TFile(outdir+'/density_correction'+suffix+'.root','recreate')

l,std_dev = [],[]
for i in range(len(h)):
  mean = sum(h[i][k] for k in range(nbin/4+1,3*nbin/4+1))/float(nbin/2)
  if mean!=0: h[i].Scale(1./mean)
  std_dev.append([i/5,i%5,round((sum((h[i][k]-1)**2 for k in range(nbin/4+1,3*nbin/4+1))/float(nbin/2))**0.5,5)])
  ffit = TF1('ffit','2*[0]*x^2*(-[3] + x^2)*([1]*x^2 + [2] + x^4) + 2*[0]*x^2*(x^2 - 0.25)*([1]*x^2 + [2] + x^4) + [0]*x*(-[3] + x^2)*(x^2 - 0.25)*(2*[1]*x + 4*x^3) + [0]*(-[3] + x^2)*(x^2 - 0.25)*([1]*x^2 + [2] + x^4)+[4]',-0.5,0.5)
  h[i].Fit(ffit,'lq','',0.01,0.49)
  h[i].Fit(ffit,'lq','',0.01,0.49)
  h[i].Fit(ffit,'q','',0.01,0.49)
  l.append([h[i].GetName()+'_0']+['%.2e'%(ffit.GetParameter(j)) for j in range(5)])
  h[i].Fit(ffit,'lq','',-0.49,-0.01)
  h[i].Fit(ffit,'lq','',-0.49,-0.01)
  h[i].Fit(ffit,'q','',-0.49,-0.01)
  l.append([h[i].GetName()+'_1']+['%.2e'%(ffit.GetParameter(j)) for j in range(5)])

names = ['lg_out','lg_in','pwo_out','pwo_center','pwo_in']
hregions = listofhisto('h',[['x','y'],regions,range(10)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
count = [0 for i in range(5)]
for i in range(len(h)):
  hname = h[i].GetName()
  i1 = 0 if 'hx' in hname else 1
  i3 = int(hname[-1])
  j = int(hname[hname.index('_')+1:hname.rindex('_')])
  i2 = regions_id(j)
  hregions[i1][i2][i3].Add(h[i])
  count[i2]+=1

lregions,std_dev_regions = [],[]
for x in unnest(hregions):
  mean = sum(x[k] for k in range(nbin/4+1,3*nbin/4+1))/float(nbin/2)
  if mean!=0: x.Scale(1/mean)
  std_dev_regions.append([x.GetName(),round((sum((x[k]-1)**2 for k in range(nbin/4+1,3*nbin/4+1))/float(nbin/2))**0.5,5)])
  ffit = TF1('ffit','2*[0]*x^2*(-[3] + x^2)*([1]*x^2 + [2] + x^4) + 2*[0]*x^2*(x^2 - 0.25)*([1]*x^2 + [2] + x^4) + [0]*x*(-[3] + x^2)*(x^2 - 0.25)*(2*[1]*x + 4*x^3) + [0]*(-[3] + x^2)*(x^2 - 0.25)*([1]*x^2 + [2] + x^4)+[4]',-0.5,0.5)
  x.Fit(ffit,'lq','',0.01,0.49)
  x.Fit(ffit,'lq','',0.01,0.49)
  x.Fit(ffit,'q','',0.01,0.49)
  lregions.append([x.GetName()]+['%.2e'%(ffit.GetParameter(j)) for j in range(5)])
  x.Fit(ffit,'lq','',-0.49,-0.01)
  x.Fit(ffit,'lq','',-0.49,-0.01)
  x.Fit(ffit,'q','',-0.49,-0.01)
  lregions.append([x.GetName()]+['%.2e'%(ffit.GetParameter(j)) for j in range(5)])

writelisto(h,fout,'final',style=[1,0,0],yrg=[0,1.5])
writelisto(hregions,fout,'regions',style=[1,0,0],yrg=[0,1.5])
writelist(outdir+'/density_params'+suffix+'.txt',l,ljust=10)
writelist(outdir+'/density_params_regions'+suffix+'.txt',lregions,ljust=[17]+[10 for i in range(5)])
writelist(outdir+'/density_stddev'+suffix+'.txt',std_dev,ljust=10)
writelist(outdir+'/density_stddev_regions'+suffix+'.txt',std_dev_regions,ljust=[17,10])
fout.Close()
f.Close()
