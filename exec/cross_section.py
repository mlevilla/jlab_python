#!/usr/bin/env python
from params import *

largs = [('n','','several','suffix'),
         ('r',[-1,-1],'','list_of_run'),
         ('i',workf,'','input_dir'),
         ('int',[0.8,2.],'','integration range'),
         ('spacer',False,'','spacer_efficiency')]

[suffix,lr,indir,intrange,spacer] = ArgParser(largs).argl

from misc_root import *
from selection import *

print 'initialisation'

fout = TFile(indir+'/cross_section_{0}_{1}{2}.root'.format(lr[0],lr[1],suffix[-1]),'recreate')

load_run_params(lr[0])
Ebeam = beam_energy[0]/1000.
alpham = (Ebeam-m_e)/(Ebeam+m_e)
max_energy = 1500. if Ebeam<2. else 2500.
ftmp = TFile(indir+'/yields_'+str(lr[0])+suffix[0]+'.root')
thbin = get_binning(ftmp.yields.Get('hth_ep'))
ftmp.Close()
q2bin = [ep_q2(th*degrad,Ebeam) if th!=0 else 0. for th in thbin]
nbin = len(thbin)-1

idx_intrange = [0,0]
for i in range(nbin):
  if thbin[i-1]<intrange[0]<=thbin[i]: idx_intrange[0] = i+1
  if thbin[i]<=intrange[1]<thbin[i+1]: idx_intrange[1] = i+1
intrange = [thbin[idx_intrange[0]-1],thbin[idx_intrange[1]-1]]

# target and beam
load_live_charge()
load_thickness() 

# periods
periods = get_bg_runs(lr[0],lr[1])

# barycenter histograms
hbarth = listofhisto('hbarth',['ep','ee1','ee2','ee3'],nbin,xbin=thbin,title=';#theta (deg);#theta (deg)')
hbarq2 = listofhisto('hbarq2',[],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});Q^{2} (GeV^{2})')
# selection histograms
helas1 = listofhisto('helas1',[['prod','empty'],['ep','ee1','ee2','ee3'],range(2)],1000,-1,1,title=';E/E_{theo}-1;')
helas2 = listofhisto('helas2',[['prod','empty'],['ee2','ee3'],range(5)],1000,-1,1,title=';(#sum E)/E_{beam}-1;')
hdphi = listofhisto('hdphi',[['prod','empty'],['ee2','ee3'],range(5)],1000,-30,30,title=';#Delta#phi (deg);')
hzvertex = listofhisto('hz',[['prod','empty'],['ee2','ee3'],range(5)],1000,-1000,1000,title=';z_{vertex} (mm);')
h2cut1 = listofhisto('h2cut1',[['prod','empty'],['ep','ee1'],range(2)],200,0.,10.,150,0.,max_energy,title=';#theta (deg);E (MeV)')
h2cut2 = listofhisto('h2cut2',[['prod','empty'],['ee2','ee3'],range(5)],200,0.,10.,150,0.,max_energy,title=';#theta (deg);E (MeV)')
hxy = listofhisto('hxy',[['prod','empty'],['ep','ee1','ee2','ee3'],['hycal','gem']],1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hselection = [helas1,helas2,hdphi,hzvertex,h2cut1,h2cut2,hxy]

# kinematic histograms
names = ['prod','empty','diff','ratio','cross','eff']
snames = ['cross','eff']
ynames = ['N','n','N-c#cdot n','n/(N-c#cdot n)','#frac{d#sigma}{d#Omega} (b/sr)','#frac{d#sigma}{d#Omega} (b/sr)']
enames = ['ep','ee1','ee2','ee3']

hth_yield = listofhisto('hth_yield',[names,enames,range(len(periods))],nbin,xbin=thbin,title=[';#theta (deg);'+x for x in ynames])
hrth_yield = listofhisto('hrth_yield',[enames[1:],range(len(periods))],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hth = listofhisto('hth',[names,enames],nbin,xbin=thbin,title=[';#theta (deg);'+x for x in ynames])
hrth = listofhisto('hrth',[snames,enames[1:]],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')
hrth2 = listofhisto('hrth2',[snames,enames[1:]],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hth_int = listofhisto('hth_int',[names,enames],1,xbin=intrange,title=[';#theta (deg);'+x for x in ynames])
hrth_int = listofhisto('hrth_int',[snames,enames[1:]],1,xbin=intrange,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hq2 = listofhisto('hq2_ep',snames,nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma}{d#Omega} (b/sr)')
hrq2 = listofhisto('hrq2',[snames,enames[1:]],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma_{ep}}{d#sigma_{ee}}')
hrq22 = listofhisto('hrq22',[snames,enames[1:]],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma_{ep}}{d#sigma_{ee}}')

# pulls
pth = listofhisto('pth',[names,enames],100,-10,10,title=';pulls;')
prth = listofhisto('prth',[enames[1:]],100,-10,10,title=';pulls;')
  
# live charge and thickness
lcharge = [[sum([live_charge[run]*1e4 for run in x]) for x in y] for y in periods]
ratio_lcharge = [x[0]/x[1] for x in lcharge]
sum_lcharge = sum(x[0] for x in lcharge)
thick = [[[sum([live_charge[run]*1e4*y[run] for run in x])/lcharge[i][j] if lcharge[i][j]!=0 else 0. for y in [thickness,dthickness]] for j,x in enumerate(z)] for i,z in enumerate(periods)]
average_thickness = sum(x[0][0] for x in thick)/len(thick)

# gem efficiency
st_energy = '1GeV' if Ebeam<2. else '2GeV'
st_spacer = 'remove' if spacer else 'with'
[gem_eff,gem_deff] = readfile(workf+'/gem_efficiencies/{0}_ep_{1}_spacer.txt'.format(st_energy,st_spacer),[float,float],cols=[4,5],tr=True)
hth_gem_eff = listofhisto('hth_gem_eff',[],nbin,xbin=thbin,title=';#theta (deg); efficiency')
for i in range(nbin):
  hth_gem_eff[i+1] = gem_eff[i]
  hth_gem_eff.SetBinError(i+1,gem_deff[i])

# input files
fs = [[[TFile(indir+'/yields_'+str(run)+suffix[0]+'.root') for run in x] for x in y] for y in periods]

# barycenter
print 'getting histograms and yields'
nrun = sum(len(p[0]) for p in periods)
for i in range(len(periods)):
  for m in range(len(periods[i][0])):
    if periods[i][0][m]>lr[1] or periods[i][0][m]<lr[0]: continue
    u = fs[i][0][m].Get('barycenter')
    for h in hbarth+[hbarq2]:
      h.Add(u.Get(h.GetName()))
for h in hbarth+[hbarq2]: h.Scale(1./nrun)
    
# selection
for i in range(2):
  for j in range(len(periods)):
    for m in range(len(periods[j][i])):
      if periods[j][i][m]>lr[1] or periods[j][i][m]<lr[0]: continue
      u = fs[j][i][m].Get('selection')
      u2 = fs[j][i][m].Get('occupancy')
      for h1 in hselection:
        for h2 in h1[i]:
          for h3 in h2:
            name_temp = h3.GetName().replace('_empty','').replace('_prod','')
            if 'hxy' in h3.GetName(): h3.Add(u2.Get(name_temp))
            else: h3.Add(u.Get(name_temp))

# yields
for i in range(2):
  for j in range(len(periods)):
    for k in range(len(periods[j][i])):
      if periods[j][i][k]>lr[1] or periods[j][i][k]<lr[0]: continue
      u = fs[j][i][k].Get('yields')
      for m,suf in enumerate(['ep','ee1','ee2','ee3']):
        hth_yield[i][m][j].Add(u.Get('hth_'+suf))

# background subtraction and summation
print 'background subtraction'
for i in range(4):
  for j in range(len(periods)):
    hth_yield[2][i][j].Add(hth_yield[0][i][j],hth_yield[1][i][j],1.,-ratio_lcharge[j])
    hth_yield[3][i][j].Divide(hth_yield[1][i][j],hth_yield[2][i][j])
    hth_yield[4][i][j].Add(hth_yield[2][i][j])
    hth_yield[4][i][j].Scale(1./lcharge[j][0]*1e28/thick[j][0][0]/6.242e9)
    hth_yield[5][i][j].Divide(hth_yield[4][i][j],hth_gem_eff)
    for k in range(nbin):
      hth_yield[4][i][j][k+1] = hth_yield[4][i][j][k+1]/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad))
      hth_yield[4][i][j].SetBinError(k+1,hth_yield[4][i][j].GetBinError(k+1)/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad)))
      hth_yield[5][i][j][k+1] = hth_yield[5][i][j][k+1]/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad))
      hth_yield[5][i][j].SetBinError(k+1,hth_yield[5][i][j].GetBinError(k+1)/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad)))
    if i!=0:
      hrth_yield[i-1][j].Divide(hth_yield[4][0][j],hth_yield[4][i][j])
    hth[0][i].Add(hth_yield[0][i][j])
    hth[1][i].Add(hth_yield[1][i][j],ratio_lcharge[j])
    hth[2][i].Add(hth_yield[2][i][j])
  
  hth[3][i].Divide(hth[1][i],hth[0][i])
  # correlation using several times same empty run
  for j in range(len(periods)-1):
    s = set(periods[j][1])-set(periods[j+1][1])
    for r in s:
      idx_r = periods[j][1].index(r)
      h = fs[j][1][idx_r].Get('yields').Get('hth_'+['ep','ee1','ee2','ee3'][i])
      for k in range(nbin):
        hth[2][i].SetBinError(k+1,(hth[2][i].GetBinError(k+1)**2+2*ratio_lcharge[j]*ratio_lcharge[j+1]*h[k+1])**0.5)

# integrated yields
fout.cd()
for i in range(4):
  for j in range(3):
    hth_int[j][i] = hth[j][i].Rebin(1,hth[j][i].GetName().replace('hth','hth_int'),array('d',intrange))
                 
# cross sections
print 'cross section calculation'
for i in range(4):
  hth[4][i].Add(hth[2][i])
  hth[4][i].Scale(1./sum_lcharge*1e28/average_thickness/6.242e9)
  for k in range(nbin):
    hth[4][i][k+1] = hth[4][i][k+1]/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad))
    hth[4][i].SetBinError(k+1,hth[4][i].GetBinError(k+1)/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad)))
  # efficiency correction
  hth[5][i].Divide(hth[4][i],hth_gem_eff)
  # integrated cross section
  for j in range(2):
    deff_temp = 0.
    for k in range(*idx_intrange):
      if hth[4+j][i][k]!=0:
        hth_int[4+j][i][1] += hth[4+j][i][k]*sin(hbarth[i][k]*degrad)*(thbin[k]-thbin[k-1])*degrad
        deff_temp += (hth[4+j][i].GetBinError(k)*sin(hbarth[i][k]*degrad)*(thbin[k]-thbin[k-1])*degrad)**(-2)
    if deff_temp!=0: hth_int[4+j][i].SetBinError(1,deff_temp**(-0.5))
    hth_int[4+j][i].Scale(1./(cos(intrange[0]*degrad)-cos(intrange[1]*degrad)))
    if i==0: continue
    # ratio ep/moller
    hrth[j][i-1].Divide(hth[4+j][0],hth[4+j][i])
    if hth_int[4+j][i][1]!=0: hrth2[j][i-1].Add(hth[4+j][0],1./hth_int[4+j][i][1])
    hrth_int[j][i-1].Divide(hth_int[4+j][0],hth_int[4+j][i])

# q2 cross-sections
for j in range(2):
  for k in range(nbin):
    hq2[j][k+1] = hth[4+j][0][k+1]
    hq2[j].SetBinError(k+1,hth[4+j][0].GetBinError(k+1))
    for i in range(1,4):
      hrq2[j][i-1][k+1] = hrth[j][i-1][k+1]
      hrq2[j][i-1].SetBinError(k+1,hrth[j][i-1].GetBinError(k+1))
      hrq22[j][i-1][k+1] = hrth2[j][i-1][k+1]
      hrq22[j][i-1].SetBinError(k+1,hrth2[j][i-1].GetBinError(k+1))

# pulls
for i in range(4):
  for j in range(6):
    for m in range(len(periods)):
      for k in range(nbin):
        # if hth_yield[j][i][m].GetBinError(k+1)**2-hth[j][i].GetBinError(k+1)**2>0:
          # pth[j][i].Fill((hth_yield[j][i][m][k+1]-hth[j][i][k+1])/(hth_yield[j][i][m].GetBinError(k+1)**2-hth[j][i].GetBinError(k+1)**2)**0.5)
        if hth_yield[j][i][m].GetBinError(k+1)!=0:
          pth[j][i].Fill((hth_yield[j][i][m][k+1]-hth[j][i][k+1])/hth_yield[j][i][m].GetBinError(k+1))
        if j==0 and i!=0 and hrth_yield[i-1][m].GetBinError(k+1)!=0:
          # if hrth_yield[i-1][m].GetBinError(k+1)**2-hrth[i-1].GetBinError(k+1)**2>0:
          #   prth[i-1].Fill((hrth_yield[i-1][m][k+1]-hrth[i-1][k+1])/(hrth_yield[i-1][m].GetBinError(k+1)**2-hrth[i-1].GetBinError(k+1)**2)**0.5)
          prth[i-1].Fill((hrth_yield[i-1][m][k+1]-hrth[0][i-1][k+1])/hrth_yield[i-1][m].GetBinError(k+1))

print 'graph extraction'

# strings for titles
thn = ';#theta (deg);'
q2n = ';Q^{2} (GeV^{2});'
xsn = '#frac{d#sigma}{d#Omega} (b/sr)'
rxsn = '#frac{d#sigma_{ep}}{d#sigma_{ee}}'

# theoretical graphs
ftheo = TFile(dataf+'/graph_theo0.root')
gomtheo = [[ftheo.Get('gom_{0}_{1}_{2}'.format(int(Ebeam),x,y)) for x in ['ep','moller']] for y in ['elas','brem']]
gththeo = [[ftheo.Get('gth_{0}_{1}_{2}'.format(int(Ebeam),x,y)) for x in ['ep','moller']] for y in ['elas','brem']]
gq2theo = [tgraph([ep_q2(th*degrad,Ebeam) for th in x[0].GetX()],x[0],name='gq2_theo_'+['elas','brem'][i],title=q2n+xsn) for i,x in enumerate(gomtheo)]
grththeo = [tgraph(x[0],op=['/',x[1]],name='grththeo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]
grq2theo = [tgraph(gq2theo[i],x[0],op=['/',x[1]],name='grq2theo_'+['elas','brem'][i],title=q2n+rxsn) for i,x in enumerate(gomtheo)]
theo_int = [[y.Integral(int((intrange[0]-0.5)/0.01+1),int((intrange[1]-0.5)/0.01+1))/(cos(intrange[0]*degrad)-cos(intrange[1]*degrad))/2./pi for y in x] for x in gththeo]
grth2theo = [tgraph(x[0],x[0],op=['/',[theo_int[i][1] for _ in range(x[0].GetN())]],name='grth2theo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]
grq22theo = [tgraph(gq2theo[i],x[0],op=['/',[theo_int[i][1] for _ in range(x[0].GetN())]],name='grq22theo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]

# experimental graphs
gthexp = [[tgraph(hbarth[i],hth[4+j][i],name='gthexp_'+['','eff_'][j]+x,title=thn+xsn) for i,x in enumerate(enames)] for j in range(2)]
gq2exp = [tgraph(hbarq2,hq2[j],name='gq2exp_'+['','eff_'][j]+'ep',title=q2n+xsn) for j in range(2)]

# experimental ratio graphs
grthexp = [[[tgraph(hbarth[0],h[j][i],name='grth'+['','2'][k]+'exp_'+['','eff_'][j]+x,title=thn+rxsn) for i,x in enumerate(enames[1:])] for j in range(2)] for k,h in enumerate([hrth,hrth2])] 
grq2exp = [[[tgraph(hbarq2,h[j][i],name='grq2'+['','2'][k]+'exp_'+['','eff_'][j]+x,title=q2n+rxsn) for i,x in enumerate(enames[1:])] for j in range(2)] for k,h in enumerate([hrq2,hrq22])] 

# remove points with low stats / weird behaviour
for gth,gq2 in zip(unnest(grthexp[0]),unnest(grq2exp[0])):
  for i in range(gth.GetN()-1):
    if abs(gth.GetY()[i+1]-gth.GetY()[i])/(gth.GetX()[i+1]-gth.GetX()[i])>2:
      gth.SetPoint(i+1,gth.GetX()[i+1],0.)
      gth.SetPointError(i+1,gth.GetEX()[i+1],0.)
      gq2.SetPoint(i+1,gq2.GetX()[i+1],0.)
      gq2.SetPointError(i+1,gq2.GetEX()[i+1],0.)

# comparison experimental/theory
gthcomp = [[tgraph(hbarth[0],h[i],op=['/',[grththeo,grth2theo][j][0]],name='gth'+['','2'][j]+'comp_'+x+'_elas',title=thn+xsn+'_{data}/'+xsn+'_{theo}') for i,x in enumerate(enames[1:])] for j,h in enumerate([hrth[0],hrth2[1]])]
gq2comp = [[tgraph(hbarq2,h[i],op=['/',[grq2theo,grq22theo][j][0]],name='gq2'+['','2'][j]+'comp_'+x+'_elas',title=thn+xsn+'_{data}/'+xsn+'_{theo}') for i,x in enumerate(enames[1:])] for j,y in enumerate([hrq2[0],hrq22[1]])]

                      
# writing
print 'writing  histogramms and graphs'
writelisto(hselection,fout,['selection'])
writelisto([hbarth,hbarq2],fout,['barycenter'])
writelisto([hth,hth_int],fout,['theta'])
writelisto([hrth,hrth2,hrth_int,gthexp,grthexp,gomtheo,grththeo],fout,['theta'],yrg=[0.,1.1])
writelisto(gthcomp,fout,['theta'],yrg=[0.5,2.])
writelisto(hq2,fout,['q2'])
writelisto([hq2,hrq2,gq2exp,grq2exp,gq2theo,grq2theo],fout,['q2'],yrg=[0,1.1])
writelisto(gq2comp,fout,['q2'],yrg=[0.5,2.])
writelisto([pth,prth],fout,['pulls'])
  
fout.Close()
sys.exit()
