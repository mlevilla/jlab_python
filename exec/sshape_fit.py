#!/usr/bin/env python
from params import *
from moddensity_utils import *
from ROOT import TF2

largs = [('i','','','input_file'),
         ('o','.','','output_folder'),
         ('n','','','suffix'),
         ('e',2,'','beam_energy')]

aparser = ArgParser(largs)
[infile,outdir,suffix,gev] = aparser.argl

from misc_root import *
f = TFile(infile) 
h = getlisto(f.dist)
hE = getlisto(f.Edist)
correct = 'Edist_correct' in [x.GetName() for x in f.GetListOfKeys()]
if correct: hEc = getlisto(f.Edist_correct)
nbin = h[0].GetNbinsX()

if suffix=='':
  i1 = infile.index('sshape')+len('sshape')
  i2 = infile.index('.root')
  suffix = infile[i1:i2]

fout = TFile(outdir+'/sshape_fit'+suffix+'.root','recreate')

region_names = ['pwo','lg_top','lg_right','lg_bottom','lg_left','lg']
h_regions = listofhisto('h',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-0.5,0.5,nbin,-0.5,0.5,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;')
hE_regions = listofhisto('hE',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-0.5,0.5,nbin,-0.5,0.5,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;E\'/E_{theo}')

if correct: hEc_regions = listofhisto('hEc',[x+'_'+y for x in region_names for y in ['ep','ee']],nbin,-0.5,0.5,nbin,-0.5,0.5,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;E\'/E_{theo}')

ffit = TF2('ffit','[0]*(1+[1]*x^2+[2]*y^2+[3]*x^2*y^2+[4]*x^4+[5]*y^4+[6]*x+[7]*y)',-0.5,0.5,-0.5,0.5)
ffit.SetParameters(1,0,0,0,0,0,0,0)
ffit.SetParLimits(0,0.9,1.5)

lindex = readfile(conff+'/'+str(gev)+'GeV_ep_groupList.txt',[int],cols=[2],start=2)+readfile(conff+'/'+str(gev)+'GeV_ee_groupList.txt',[int],cols=[2],start=2)
nindex = len(lindex)
lep = readfile(conff+'/groupindex_'+str(gev)+'GeV_ep.txt',[int]*2,out=dict)
lee = readfile(conff+'/groupindex_'+str(gev)+'GeV_ee.txt',[int]*2,out=dict)
nindex_ep = max(lep.values())+1

print nindex
for i in range(len(h)):
  name = h[i].GetName()
  i1 = name.rindex('_')+1
  index = int(name[i1:])
  ep = int('ee' in name)
  sector = which_sector(lindex[i])
  h_regions[sector*2+ep].Add(h[i])
  hE_regions[sector*2+ep].Add(hE[i]) 
  if correct: hEc_regions[sector*2+ep].Add(hEc[i])
  if sector>0:
    h_regions[10+ep].Add(h[i])
    hE_regions[10+ep].Add(hE[i]) 
    if correct: hEc_regions[10+ep].Add(hEc[i])

lparam = [[] for i in range(nindex)]
for i in range(len(h)):
  name = h[i].GetName()
  i1 = name.rindex('_')+1
  index = name[i1:]
  ep = int('ep' not in name)
  hE[i].Divide(h[i])
  if correct: hEc[i].Divide(h[i])
  if i in [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 344, 345, 360, 361, 564, 565, 566, 567, 587, 588, 607, 608, 627, 628, 647, 648, 668, 669, 670, 671, 760, 761, 795, 796]:
    ffit.ReleaseParameter(6)
    ffit.ReleaseParameter(7)
  else:
    ffit.FixParameter(6,0)
    ffit.FixParameter(7,0)
  if h[i].GetEntries()>200: hE[i].Fit(ffit,'q')
  if hE[i].GetFunction('ffit'):
    lparam[i] = [name[name.index('e'):]]+[hE[i].GetFunction('ffit').GetParameter(k) for k in range(8)]
  else:
    lparam[i] = [name[name.index('e'):]]+[1,0,0,0,0,0,0,0]

if gev==1:
  ep_inner = [208, 209, 210, 211, 221, 222, 237, 238, 247, 248, 263, 264, 274, 275, 276, 277]
  h2D_ep = listofhisto('h2D_ep',ep_inner,nbin,-0.5,0.5,nbin,-0.5,0.5,title=';(x_{hycal}-x_{center})/d;(y_{hycal}-y_{center})/d;')
  h1 = []
  s = TSpectrum(2)
  f2gaus = TF1('2gaus','gaus(0)+gaus(3)+[6]',0.9,1.1)
  for k,x in enumerate(ep_inner):
    h3 = f.inner.Get('h3D_ep_'+str(x))
    for i in range(1,nbin+1):
      for j in range(1,nbin+1):
        h1.append(h3.ProjectionZ('h1_'+str(x)+'_ep_pz_'+str(i)+'_'+str(j),i,i,j,j))
        s.Search(h1[-1])
        if s.GetNPeaks()<2 or s.GetPositionY()[1]<0.25*s.GetPositionY()[0] or abs(s.GetPositionX()[0]-s.GetPositionX()[1])<0.02:
          l = fit_gaus(h1[-1],sigma=1)
          h2D_ep[k].SetBinContent(i,j,l[0])
        else:
          index = 0 if s.GetPositionX()[0]>s.GetPositionX()[1] else 1
          xp = [s.GetPositionX()[index],s.GetPositionX()[1-index]]
          yp = [s.GetPositionY()[index],s.GetPositionY()[1-index]]
          f2gaus.SetParameter(0,yp[0])
          f2gaus.SetParLimits(0,0.9*yp[0],1.1*yp[0])
          f2gaus.SetParameter(1,xp[0])
          f2gaus.SetParLimits(1,xp[0]-0.02,xp[0]+0.02)
          f2gaus.SetParameter(2,0.025)
          f2gaus.SetParLimits(2,0.02,0.05)
          f2gaus.SetParameter(3,yp[1])
          f2gaus.SetParLimits(3,0.9*yp[1],1.1*yp[1])
          f2gaus.SetParameter(4,xp[1])
          f2gaus.SetParLimits(4,xp[1]-0.02,xp[1]+0.02)
          f2gaus.SetParameter(5,0.025)
          f2gaus.SetParLimits(5,0.02,0.05)
          h1[-1].Fit(f2gaus)
          h2D_ep[k].SetBinContent(i,j,f2gaus.GetParameter(1))
        
  ffit.ReleaseParameter(6)
  ffit.ReleaseParameter(7)
  for i in range(len(h2D_ep)):
    name = h2D_ep[i].GetName()
    i1 = name.rindex('_')+1
    index = name[i1:]
    h2D_ep[i].Fit(ffit,'q')
    lparam[int(index)] = [name[name.index('e'):]]+[h2D_ep[i].GetFunction('ffit').GetParameter(k) for k in range(8)]

lparam_mod_ep, lparam_mod_ee = [],[]
for x,y in module_names.items():
  if lep[x]!=-1: lparam_mod_ep.append([y]+lparam[lep[x]][1:])
  else: lparam_mod_ep.append([y,1,0,0,0,0,0,0,0])
  if lee[x]!=-1: lparam_mod_ee.append([y]+lparam[lee[x]+nindex_ep][1:])
  else: lparam_mod_ee.append([y,1,0,0,0,0,0,0,0])
lparam_regions = [[] for i in range(12)]  
# regions
ffit.SetParameters(1,0,0,0,0,0,0,0)
ffit.FixParameter(6,0)
ffit.FixParameter(7,0)
for i in range(12):
  hE_regions[i].Divide(h_regions[i])
  if correct: hEc_regions[i].Divide(h_regions[i])
  if h_regions[i].GetEntries()>200: hE_regions[i].Fit(ffit,'q')
  if hE_regions[i].GetFunction('ffit'):
    lparam_regions[i/2+6*(i%2)] =  [region_names[i/2]+['_ep','_ee'][i%2]]+[hE_regions[i].GetFunction('ffit').GetParameter(k) for k in range(8)]
  else:
    lparam_regions[i/2+6*(i%2)] = [1,0,0,0,0,0,0,0]

writelisto(hE,fout,'Edist',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.9,1.1])
writelisto(hE_regions,fout,'regions',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.95,1.05])
if gev==1:
  writelisto(h2D_ep,fout,'inner',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.98,1.02])
  writelisto(h1,fout,'inner_proj',xrg=[0.9,1.1])
if correct:
  writelisto(hEc,fout,'Edist_correct',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.9,1.1])
  writelisto(hEc_regions,fout,'regions_correct',xrg=[-0.5,0.5],yrg=[-0.5,0.5],zrg=[0.95,1.05])
writelist(outdir+'/sshape_params'+suffix+'.txt',lparam,ljust=10,ro=5)
writelist(outdir+'/sshape_params'+suffix+'_ep.txt',lparam_mod_ep,ljust=10,ro=5)
writelist(outdir+'/sshape_params'+suffix+'_ee.txt',lparam_mod_ee,ljust=10,ro=5)
writelist(outdir+'/sshape_params_regions'+suffix+'.txt',lparam_regions,ljust=[15]+[10]*8,ro=5)
fout.Close()
f.Close()
