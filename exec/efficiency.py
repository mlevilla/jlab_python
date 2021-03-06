#!/usr/bin/env python

from params import *

largs = [('i',workf,'','input_directory'),
         ('n','','several','suffix')]
[indir,suffix] = ArgParser(largs).argl

from misc_root import *
gStyle.SetOptStat(1)
f = TFile(indir+'/histo'+suffix[0]+'.root')
fout = TFile(indir+'/eff'+suffix[-1]+'.root','recreate')

load_names()

heg = {k:[f.Get(x).Get('heg_'+x+'_'+str(j)) for j in [1,2,5]] if f.Get(x) else [] for k,x in module_names.items()}

l = {k:[[y.Integral(201+50*i,250+50*i) for i in range(17)] for y in z] for k,z in heg.items()}
stat = [{x:[int(y[j][k]) for k in range(17)] if y!=[] else [0 for k in range(17)] for x,y in l.items()} for j in range(3)]

stat12 = {x:[a+b for a,b in zip(stat[0][x],stat[1][x])] for x in module_names.keys()}
stat125 = {x:[a+b+c for a,b,c in zip(stat[0][x],stat[1][x],stat[2][x])] for x in module_names.keys()}

for i in range(3): writelist(indir+'/stat_tr'+['1','2','5'][i]+suffix[-1]+'.txt',[[module_names[x]]+y for x,y in stat[i].items()])

##############################
#### efficiency by module ####
##############################

eff = {x:[float(a)/b if b!=0 else 0 for a,b in zip(stat12[x],stat125[x])] for x in module_names}
deff = {x:[((a+1.)*(a+2.)/(b+2.)/(b+3.)-(a+1.)**2/(b+2.)**2)**0.5 for a,b in zip(stat12[x],stat125[x])] for x in module_names}

writelist(indir+'/eff_modules'+suffix[-1]+'.txt',[[module_names[x]]+y for x,y in eff.items()])
writelist(indir+'/deff_modules'+suffix[-1]+'.txt',[[module_names[x]]+y for x,y in deff.items()])

# mean
mean_eff = {x:sum([y[i] for i in range(5,15) if y[i]>0])/len([y[i] for i in range(5,15) if y[i]>0]) if len([y[i] for i in range(5,15) if y[i]>0])!=0 else 0. for x,y in eff.items()}
sum_stat125 = {x:float(sum([y[i] for i in range(5,15)])) for x,y in stat125.items()}
sum_stat12 = {x:float(sum([y[i] for i in range(5,15)])) for x,y in stat12.items()}
err_stat = {x:((sum_stat12[x]+1)*(sum_stat12[x]+2)/(sum_stat125[x]+2)/(sum_stat125[x]+3)-(sum_stat12[x]+1)**2/(sum_stat125[x]+2)**2)**0.5 for x in module_names}

gmean_eff = tgraph(mean_eff.keys(),mean_eff.values(),name='gmean_eff')
gmean_eff.Write()

#######################
#### normalization ####
#######################

zones = [range(1,25)+range(31,55)+range(61,85)+range(91,115),[x for x in range(1001) if ((x not in (range(1,25)+range(31,55)+range(61,85)+range(91,115))) and (x in module_names))],range(1001,1478),list(set(range(1478,1749))-{1561,1562,1595,1596}),range(1749,1919),range(1919,2157)]

h_zones = listofhisto('h_zones',range(len(zones)),1000,0.9,1.1)
for i in module_names.keys():
  for j in range(len(zones)):
    if i in zones[j] and err_stat[i]!=0: 
      _ = h_zones[j].Fill(mean_eff[i],1/err_stat[i]**2)


zone_means = [0. for i in range(len(zones))]
zone_sigmas = [0. for i in range(len(zones))]
for i in range(len(zones)):
  r = h_zones[i].Fit('gaus','rwsq')
  zone_means[i] = r.GetParams()[1]
  zone_sigmas[i] = r.GetParams()[2]

writelisto(h_zones,fout)

choice_mean = (zone_means[1]+zone_means[-1])/2.
choice_sigma = (zone_sigmas[1]+zone_sigmas[-1])/2.
mean_eff1, eff1, deff1 = {},{},{}
for i in module_names.keys():
  for j in range(len(zones)):
    if i in zones[j]:
      mean_eff1[i] = choice_mean+choice_sigma*(mean_eff[i]-zone_means[j])/zone_sigmas[j]
      eff1[i] = [eff[i][k]*mean_eff1[i]/mean_eff[i] if mean_eff[i]!=0 else 0 for k in range(17)]
      deff1[i] = [deff[i][k]*mean_eff1[i]/mean_eff[i] if mean_eff[i]!=0 else 0 for k in range(17)]
      break

# correction for modules with no stat
correct_av_lg = [10,101,102,129,130,131,132,133,134,135,136,137,159,393,394,571,783,784,793,794,795,813,820,832,900]
not_working_lg = [900]
correct_av_pwo = [1393,1419,1420,1470,1471,1472,1527,1528,1560,1563,1629,1630,1637,1752,1835,1891]
not_working_pwo = [1835,1891]
for i in correct_av_lg+correct_av_pwo:
  l = []
  step = 34 if i>900 else 30
  for j in [i-step-1,i-step,i-step+1,i-1,i+1,i+step-1,i+step,i+step+1]:
    if j in correct_av_lg+correct_av_pwo: continue
    if j<0 or (i<=900 and j>900) or (i>900 and j<900): continue
    if j not in module_names: continue
    l.append(j)
  if i in not_working_lg+not_working_pwo: 
    mean_eff1[i] = 0.
    eff1[i] = [0 for k in range(17)]
    deff1[i] = [0 for k in range(17)]
  else: 
    mean_eff1[i] = sum([mean_eff1[j] for j in l])/len(l)
    eff1[i] = [sum([eff1[j][k] for j in l])/len(l) for k in range(17)]
    deff1[i] = [sum([eff1[j][k] for j in l])/len(l) for k in range(17)]

writelist(indir+'/eff_modules'+suffix[-1]+'_corr.txt',[[module_names[x]]+y for x,y in eff1.items()])
writelist(indir+'/deff_modules'+suffix[-1]+'_corr.txt',[[module_names[x]]+y for x,y in deff1.items()])
writelist(indir+'/mean_eff_modules'+suffix[-1]+'_corr.txt',[[module_names[x],y] for x,y in mean_eff1.items()])

gmean_eff1 = tgraph(mean_eff1.keys(),mean_eff1.values(),name='gmean_eff1')
gmean_eff1.Write()


##############################
#### efficiency by region ####
##############################

region_names = ['center_pwo','trans','center_lg']
regions = [get_regions(x) for x in region_names]

h = listofhisto('heff',[range(3),range(17)],2000,0.7,1.1)
for i in module_names.keys():
  for j in range(3):
    if i in regions[j]: 
      for k in range(17):
        if deff1[i][k]!=0: _ = h[j][k].Fill(eff1[i][k],1/deff1[i][k]**2)

region_means = [[0. for k in range(17)] for i in range(3)]
region_sigmas = [[0. for k in range(17)] for i in range(3)]
region_ent = [[0. for k in range(17)] for i in range(3)]
region_range = [[[] for k in range(17)] for i in range(3)]
c = TCanvas()
for i in range(3):
  for k in range(17):
    h[i][k].GetXaxis().SetRangeUser(0.97,1.01)
    h[i][k].Draw('hist')
    c.Update()
    raw_input()
    if h[i][k].Integral()==0: continue
    r = h[i][k].Fit('gaus','rsl','',0.97,1.01)
    h[i][k].GetFunction('gaus').Draw('same')
    c.Update()
    raw_input()
    if h[i][k].Integral()==0: continue
    region_means[i][k] = r.GetParams()[1]
    region_sigmas[i][k] = r.GetParams()[2]
    region_ent[i][k] = h[i][k].GetEffectiveEntries()
    h[i][k].UseCurrentStyle()
    region_range[i][k] = [h[i][k].GetBinLowEdge(h[i][k].GetXaxis().GetFirst()),h[i][k].GetXaxis().GetBinUpEdge(h[i][k].GetXaxis().GetLast())]


region_sigmas2 = [[region_sigmas[i][k]/(region_ent[i][k]-1)**0.5 for k in range(17)] for i in range(3)]

writelisto(h,fout,['regions'])
Egbin = [225+50*i for i in range(17)]
gregions = [tgraph(Egbin,region_means[i],dy=region_sigmas2[i],name='g'+region_names[i],title=';E_{#gamma} (MeV);#epsilon_{trigger}') for i in range(3)]

ffit = TF1('ffit','[0]*(1-exp(-[1]*x/1000-[2]))',200,1200)
ffit.SetParameters(0.995,0.1,0.1)
ffit.SetParLimits(0,0.994,1.0)
ffit.SetParLimits(1,0.01,30)
ffit.SetParLimits(2,0,20)
for x in gregions: x.Fit(ffit)

writelisto(gregions,fout,['regions'])

writelist(indir+'/eff_region'+suffix[-1]+'.txt',[[x]+y for x,y in zip(region_names,region_means)])
writelist(indir+'/deff_region'+suffix[-1]+'.txt',[[x]+y for x,y in zip(region_names,region_sigmas)])

mean_eff1_param0 = []
regions2 = [get_regions('center_pwo')+get_regions('inner2'),get_regions('trans'),get_regions('center_lg')+get_regions('outer2')]
for x,y in module_names.items():
  for j in range(3):
    if x in regions2[j]: 
      mean_eff1_param0.append([y,mean_eff1[x],gregions[j].GetFunction('ffit').GetParameter(1),gregions[j].GetFunction('ffit').GetParameter(2)])

writelist(indir+'/hycal_trgeff_regions.txt',mean_eff1_param0,ljust=10,ro=4)

#################################
#### efficiency by tdc group ####
#################################

load_tdc_map()
htdc = listofhisto('heff_tdc',[tdc_groups,range(17)],2000,0.7,1.1)
for i in module_names.keys():
  j = tdc_groups.index(tdc_map[i])
  for k in range(17):
    if deff1[i][k]!=0: _ = htdc[j][k].Fill(eff1[i][k],1/deff1[i][k]**2)

tdc_means = [[0. for k in range(17)] for i in range(52)]
tdc_sigmas = [[0. for k in range(17)] for i in range(52)]
tdc_ent = [[0. for k in range(17)] for i in range(52)]
c = TCanvas()
for i in range(52):
  for k in range(17):
    htdc[i][k].GetXaxis().SetRangeUser(0.97,1.01)
    htdc[i][k].Draw('hist')
    if htdc[i][k].Integral()==0: continue
    r = htdc[i][k].Fit('gaus','rsl','',0.97,1.01)
    htdc[i][k].GetFunction('gaus').Draw('same')
    if htdc[i][k].Integral()==0: continue
    tdc_means[i][k] = r.GetParams()[1]
    tdc_sigmas[i][k] = r.GetParams()[2]
    tdc_ent[i][k] = htdc[i][k].GetEffectiveEntries()
    htdc[i][k].UseCurrentStyle()


tdc_sigmas2 = [[tdc_sigmas[i][k]/(tdc_ent[i][k]-1)**0.5 if tdc_ent[i][k]>1 else 0. for k in range(17)] for i in range(52)]

#writelisto(htdc,fout,['regions'])
gtdc = [tgraph(Egbin,tdc_means[i],dy=tdc_sigmas2[i],name='g'+tdc_groups[i],title=';E_{#gamma} (MeV);#epsilon_{trigger}') for i in range(52)]

for i,x in enumerate(gtdc): 
  if tdc_groups[i] in ['W19','W20','W21','W22','W23','W24']:
    x.GetXaxis().SetRangeUser(200,650)
  x.Fit(ffit)

writelisto(gtdc,fout,['tdc'])

writelist(indir+'/eff_tdc'+suffix[-1]+'.txt',[[x]+y for x,y in zip(tdc_groups,tdc_means)])
writelist(indir+'/deff_tdc'+suffix[-1]+'.txt',[[x]+y for x,y in zip(tdc_groups,tdc_sigmas)])

mean_eff1_param1 = []
for x,y in module_names.items():
  j = tdc_groups.index(tdc_map[x])
  mean_eff1_param1.append([y,mean_eff1[x],gtdc[j].GetFunction('ffit').GetParameter(1),gtdc[j].GetFunction('ffit').GetParameter(2)])

writelist(indir+'/hycal_trgeff_tdc.txt',mean_eff1_param1,ljust=10,ro=4)

fout.Close()
