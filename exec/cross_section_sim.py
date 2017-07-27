#!/usr/bin/env python
from params import *

largs = [('n','','several','suffix'),
         ('i',workf,'','input_dir'),
         ('int',[0.8,2.],'','integration range'),
         ('nf',1,'','number_of_sim_files'),
         ('ebeam',2,'','beam_energy_for_simulation'),
         ('lum',[1.,1.],'','luminosity (ub^-1)')]
print_help(largs)
[suffix,indir,intrange,nf,ebeam,luminosity] = [catch_arg(x,y,z) for x,y,z,c in largs]

from misc_root import *
from esepp import ep_q2, m_e

print 'initialisation'

fout = TFile(indir+'/cross_section_sim{0}.root'.format(suffix[-1]),'recreate')

if ebeam==2: run = 1443
else: run = 1288
load_run_params(run)
Ebeam = beam_energy[0]/1000.
alpham = (Ebeam-m_e)/(Ebeam+m_e)
max_energy = 1500. if Ebeam<2. else 2500.
ftmp = TFile(indir+'/simyields'+suffix[0]+'_0.root') if nf!=1 else TFile(indir+'/simyields'+suffix[0]+'.root')
thbin = get_binning(ftmp.yields.Get('hth_ep'))
ftmp.Close()
q2bin = [ep_q2(th*degrad,Ebeam) if th!=0 else 0. for th in thbin] 
nbin = len(thbin)-1

idx_intrange = [0,0]
for i in range(nbin):
  if thbin[i-1]<intrange[0]<=thbin[i]: idx_intrange[0] = i+1
  if thbin[i]<=intrange[1]<thbin[i+1]: idx_intrange[1] = i+1
intrange = [thbin[idx_intrange[0]-1],thbin[idx_intrange[1]-1]]

# periods
periods = range(nf)

# barycenter histograms
hbarth = listofhisto('hbarth',['ep','ee1','ee2','ee3'],nbin,xbin=thbin,title=';#theta (deg);#theta (deg)')
hbarq2 = listofhisto('hbarq2',[],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});Q^{2} (GeV^{2})')
# selection histograms
helas1 = listofhisto('helas1',[['ep','ee1','ee2','ee3'],range(2)],1000,-1,1,title=';E/E_{theo}-1;')
helas2 = listofhisto('helas2',[['ee2','ee3'],range(5)],1000,-1,1,title=';(#sum E)/E_{beam}-1;')
hdphi = listofhisto('hdphi',[['ee2','ee3'],range(5)],1000,-30,30,title=';#Delta#phi (deg);')
hzvertex = listofhisto('hz',[['ee2','ee3'],range(5)],1000,-1000,1000,title=';z_{vertex} (mm);')
h2cut1 = listofhisto('h2cut1',[['ep','ee1'],range(2)],200,0.,10.,150,0.,max_energy,title=';#theta (deg);E (MeV)')
h2cut2 = listofhisto('h2cut2',[['ee2','ee3'],range(5)],200,0.,10.,150,0.,max_energy,title=';#theta (deg);E (MeV)')
hxy = listofhisto('hxy',[['ep','ee1','ee2','ee3'],['hycal','gem']],1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hselection = [helas1,helas2,hdphi,hzvertex,h2cut1,h2cut2,hxy]

# kinematic histograms
names = ['prod','empty','diff','ratio','cross']
ynames = ['N','n','N-c#cdot n','n/(N-c#cdot n)','#frac{d#sigma}{d#Omega} (b/sr)']
enames = ['ep','ee1','ee2','ee3']

hth_yield = listofhisto('hth_yield',[names,enames,range(len(periods))],nbin,xbin=thbin,title=[';#theta (deg);'+x for x in ynames])
hrth_yield = listofhisto('hrth_yield',[enames[1:],range(len(periods))],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hth = listofhisto('hth',[names,enames],nbin,xbin=thbin,title=[';#theta (deg);'+x for x in ynames])
hrth = listofhisto('hrth',enames[1:],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')
hrth2 = listofhisto('hrth2',enames[1:],nbin,xbin=thbin,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hth_int = listofhisto('hth_int',[names,enames],1,xbin=intrange,title=[';#theta (deg);'+x for x in ynames])
hrth_int = listofhisto('hrth_int',enames[1:],1,xbin=intrange,title=';#theta (deg);#frac{d#sigma_{ep}}{d#sigma_{ee}}')

hq2 = listofhisto('hq2_ep',[],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma}{d#Omega} (b/sr)')
hrq2 = listofhisto('hrq2',enames[1:],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma_{ep}}{d#sigma_{ee}}')
hrq22 = listofhisto('hrq22',enames[1:],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});#frac{d#sigma_{ep}}{d#sigma_{ee}}')

# pulls
# pth = listofhisto('pth',[names,enames,range(nbin)],100,-10,10,title=';pulls;')
# prth = listofhisto('prth',[enames[1:],range(nbin)],100,-10,10,title=';pulls;')
  
# luminosity = [3.3693e8,7.2168e7,7.2168e7,7.2168e7]
luminosity.extend([luminosity[-1],luminosity[-1]])
luminosity = [x*1e6 for x in luminosity]

# input files
fs = [TFile(indir+'/simyields'+suffix[0]+'_'+str(i)+'.root') if nf!=1 else TFile(indir+'/simyields'+suffix[0]+'.root') for i in range(nf)]

# barycenter
print 'getting histograms and yields'
for i in range(nf):
    u = fs[i].Get('barycenter')
    for h in hbarth+[hbarq2]:
      h.Add(u.Get(h.GetName()))
for h in hbarth+[hbarq2]: h.Scale(1./nf)
    
# selection
for i in range(nf):
  u = fs[i].Get('selection')
  u2 = fs[i].Get('occupancy')
  for h1 in hselection:
    for h2 in h1:
      for h3 in h2:
        if 'hxy' in h3.GetName(): h3.Add(u2.Get(h3.GetName()))
        else: h3.Add(u.Get(h3.GetName()))

# yields
for i in range(nf):
  u = fs[i].Get('yields')
  for m,suf in enumerate(['ep','ee1','ee2','ee3']):
    hth_yield[0][m][i].Add(u.Get('hth_'+suf))

# background subtraction and summation
print 'background subtraction'
for i in range(4):
  for j in range(nf):
    hth_yield[2][i][j].Add(hth_yield[0][i][j])
    hth_yield[4][i][j].Add(hth_yield[2][i][j])
    hth_yield[4][i][j].Scale(1./luminosity[i]/nf)
    if i!=0:
      hrth_yield[i-1][j].Divide(hth_yield[2][0][j],hth_yield[2][i][j])
    hth[0][i].Add(hth_yield[0][i][j])
    hth[2][i].Add(hth_yield[2][i][j])

# integrated yields
fout.cd()
for i in range(4):
  for j in range(3):
    hth_int[j][i] = hth[j][i].Rebin(1,hth[j][i].GetName().replace('hth','hth_int'),array('d',intrange))

# cross-section
print 'cross section calculation'
for i in range(4):
  hth[4][i].Add(hth[2][i])    
  hth[4][i].Scale(1./luminosity[i]/nf)
  for k in range(nbin):
    hth[4][i][k+1] = hth[4][i][k+1]/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad))
    hth[4][i].SetBinError(k+1,hth[4][i].GetBinError(k+1)/2./pi/(cos(thbin[k]*degrad)-cos(thbin[k+1]*degrad)))
  # integrated cross section
  deff_temp = 0.
  for k in range(*idx_intrange):
    if hth[4][i][k]!=0:
      hth_int[4][i][1] += hth[4][i][k]*sin(hbarth[i][k]*degrad)*(thbin[k]-thbin[k-1])*degrad
      deff_temp += (hth[4][i].GetBinError(k)*sin(hbarth[i][k]*degrad)*(thbin[k]-thbin[k-1])*degrad)**(-2)
  if deff_temp!=0: hth_int[4][i].SetBinError(1,deff_temp**(-0.5))
  hth_int[4][i].Scale(1./(cos(intrange[0]*degrad)-cos(intrange[1]*degrad)))
  if i==0: continue
  # ratio ep/moller
  hrth[i-1].Divide(hth[4][0],hth[4][i])
  if hth_int[4][i][1]!=0: hrth2[i-1].Add(hth[4][0],1./hth_int[4][i][1])
  hrth_int[i-1].Divide(hth_int[4][0],hth_int[4][i])

# q2 cross-sections
for k in range(nbin):
  hq2[k+1] = hth[4][0][k+1]
  hq2.SetBinError(k+1,hth[4][0].GetBinError(k+1))
  for i in range(1,4):
    hrq2[i-1][k+1] = hrth[i-1][k+1]
    hrq2[i-1].SetBinError(k+1,hrth[i-1].GetBinError(k+1))
    hrq22[i-1][k+1] = hrth2[i-1][k+1]
    hrq22[i-1].SetBinError(k+1,hrth2[i-1].GetBinError(k+1))

# pulls
# for i in range(4):
#   for j in range(4):
#     for m in range(nf):
#       for k in range(nbin):
#         # if hth_yield[j][i][m].GetBinError(k+1)**2-hth[j][i].GetBinError(k+1)**2>0:
#         #   pth[j][i][k].Fill((hth_yield[j][i][m][k+1]-hth[j][i][k+1])/(hth_yield[j][i][m].GetBinError(k+1)**2-hth[j][i].GetBinError(k+1)**2)**0.5)
#         if hth_yield[j][i][m].GetBinError(k+1)!=0:
#           pth[j][i].Fill((hth_yield[j][i][m][k+1]-hth[j][i][k+1])/hth_yield[j][i][m].GetBinError(k+1))
#         if j==0 and i!=0:
#           if hrth_yield[i-1][m].GetBinError(k+1)**2-hrth[i-1].GetBinError(k+1)**2>0:
#             prth[i-1][k].Fill((hrth_yield[i-1][m][k+1]-hrth[i-1][k+1])/(hrth_yield[i-1][m].GetBinError(k+1)**2-hrth[i-1].GetBinError(k+1)**2)**0.5)

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
grththeo = [tgraph(x[0],x[0],op=['/',x[1]],name='grththeo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]
grq2theo = [tgraph(gq2theo[i],x[0],op=['/',x[1]],name='grq2theo_'+['elas','brem'][i],title=q2n+rxsn) for i,x in enumerate(gomtheo)]
theo_int = [[y.Integral(int((intrange[0]-0.5)/0.01+1),int((intrange[1]-0.5)/0.01+1))/(cos(intrange[0]*degrad)-cos(intrange[1]*degrad))/2./pi for y in x] for x in gththeo]
grth2theo = [tgraph(x[0],x[0],op=['/',[theo_int[i][1] for _ in range(x[0].GetN())]],name='grth2theo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]
grq22theo = [tgraph(gq2theo[i],x[0],op=['/',[theo_int[i][1] for _ in range(x[0].GetN())]],name='grq22theo_'+['elas','brem'][i],title=thn+rxsn) for i,x in enumerate(gomtheo)]

# experimental graphs
gthexp = [tgraph(hbarth[i],hth[4][i],name='gthsim_'+x,title=thn+xsn) for i,x in enumerate(enames)]
gq2exp = tgraph(hbarq2,hq2,name='gq2sim_ep',title=q2n+xsn)

# experimental ratio graphs
grthexp = [[tgraph(hbarth[0],h[i],name='grth'+['','2'][j]+'sim_'+x,title=thn+rxsn) for i,x in enumerate(enames[1:])] for j,h in enumerate([hrth,hrth2])]
grq2exp = [[tgraph(hbarq2,hrq22[i],name='grq2'+['','2'][j]+'sim_'+x,title=q2n+rxsn) for i,x in enumerate(enames[1:])] for j,h in enumerate([hrq2,hrq22])]

# comparison experimental/theory
gthcomp = [[tgraph(hbarth[0],h[i],op=['/',[grththeo,grth2theo][j][0]],name='gth'+['','2'][j]+'comp_'+x+'_elas',title=thn+xsn+'_{data}/'+xsn+'_{theo}') for i,x in enumerate(enames[1:])] for j,h in enumerate([hrth,hrth2])]
gq2comp = [[tgraph(hbarq2,h[i],op=['/',[grq2theo,grq22theo][j][0]],name='gq2'+['','2'][j]+'comp_'+x+'_elas',title=thn+xsn+'_{data}/'+xsn+'_{theo}') for i,x in enumerate(enames[1:])] for j,y in enumerate([hrq2,hrq22])]
                      
# writing
print 'writing  histogramms and graphs'
writelisto(hselection,fout,['selection'])
writelisto([hbarth,hbarq2],fout,['barycenter'])
writelisto([hth,hth_int,hrth,hrth2,hrth_int,gthexp,grthexp,gomtheo,grththeo,gthcomp],fout,['theta'])
writelisto([hq2,hrq2,gq2exp,grq2exp,gq2theo,grq2theo,gq2comp],fout,['q2'])
# writelisto([pth,prth],fout,['pulls'])
  
fout.Close()
sys.exit()
