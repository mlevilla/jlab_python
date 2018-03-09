# efficiency vs Energy plots 

methods = [['calib','13'],['fsnake','0']]
folders = ['/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/'.format(m[0],m[1]) for m in methods]
cuts = ['0','t40','t40_xy','t40_xy_s4','t40_xy_s3','t40_xy_s2','t40_xy_s3_n2','t40_xy_s3_n3']
cuts = ['t40_xy_s3']
f = [[TFile(fo+'/d_'+c+'/eff.root') for c in cuts] for fo in folders]
regions = ['lg_center','pwo_center']
regions = ['edge','lg_center','trans_pwo','pwo_center']
g = [[[y.Get(re).Get('g0_'+re) for re in regions] for y in x] for x in f]


# normalized efficicency map
load_names()
method,ipass = 'fsnake',0
folder = '/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/d_t40_xy_s3/'.format(method,ipass)
folder = '/work/hallb/prad/mlevilla/histo_calib/test2_3/'
stat1 = readfile(folder+'stat_tr1_3.txt',[str]+[float]*17,cols=range(18),out=dict)
stat2 = readfile(folder+'stat_tr2_3.txt',[str]+[float]*17,cols=range(18),out=dict)
stat5 = readfile(folder+'stat_tr5_3.txt',[str]+[float]*17,cols=range(18),out=dict)
stat125 = {x:[a+b+c for a,b,c in zip(stat1[y],stat2[y],stat5[y])] for x,y in module_names.items()}
stat12 = {x:[a+b for a,b in zip(stat1[y],stat2[y])] for x,y in module_names.items()}
eff = {x:[a/b if b!=0 else 0 for a,b in zip(stat12[x],stat125[x])] for x in module_names}
deff = {x:[((a+1)*(a+2)/(b+2)/(b+3)-(a+1)**2/(b+2)**2)**0.5 for a,b in zip(stat12[x],stat125[x])] for x in module_names}
mean_eff = {x:sum([y[i] for i in range(5,15) if y[i]>0])/len([y[i] for i in range(5,15) if y[i]>0]) if len([y[i] for i in range(5,15) if y[i]>0])!=0 else 0. for x,y in eff.items()}
sum_stat125 = {x:sum([y[i] for i in range(5,15)]) for x,y in stat125.items()}
sum_stat12 = {x:sum([y[i] for i in range(5,15)]) for x,y in stat12.items()}
err_stat = {x:((sum_stat12[x]+1)*(sum_stat12[x]+2)/(sum_stat125[x]+2)/(sum_stat125[x]+3)-(sum_stat12[x]+1)**2/(sum_stat125[x]+2)**2)**0.5 for x in module_names}

zones = [range(1,25)+range(31,55)+range(61,85)+range(91,115),[x for x in range(1001) if ((x not in (range(1,25)+range(31,55)+range(61,85)+range(91,115))) and (x in module_names))],range(1001,1478),list(set(range(1478,1749))-{1561,1562,1595,1596}),range(1749,1919),range(1919,2157)]

h = [TH1F('h'+str(i),'',1000,0.9,1.1) for i in range(6)]
for i in module_names.keys():
  for j in range(6):
    if i in zones[j]: _ = h[j].Fill(mean_eff[i],1/err_stat[i]**2)

zone_means = [0. for i in range(6)]
zone_sigmas = [0. for i in range(6)]
for i in range(6):
  h[i].Draw()
  s = raw_input()
  r = h[i].Fit('gaus','rws')
  zone_means[i] = r.GetParams()[1]
  zone_sigmas[i] = r.GetParams()[2]

choice_mean = 0.995
choice_sigma = 1e-3
mean_eff1, eff1, deff1 = {},{},{}
for i in module_names.keys():
  for j in range(6):
    if i in zones[j]:
      mean_eff1[i] = choice_mean+choice_sigma*(mean_eff[i]-zone_means[j])/zone_sigmas[j]
      eff1[i] = [eff[i][k]*mean_eff1[i]/mean_eff[i] if mean_eff[i]!=0 else 0 for k in range(17)]
      deff1[i] = [deff[i][k]*mean_eff1[i]/mean_eff[i] if mean_eff[i]!=0 else 0 for k in range(17)]
      break

correct_av_lg = [10,101,102,129,130,131,132,133,134,135,136,159,393,394,783,784,793,794,795,813,820,832,900]
not_working_lg = [900]
correct_av_pwo = [1393,1419,1420,1470,1471,1472,1527,1528,1560,1563,1629,1630,1637,1752,1824,1835,1891]
# correct_av_pwo = [1256,1393,1419,1420,1470,1471,1472,1527,1528,1560,1563,1597,1610,1629,1630,1752,1824,1835,1891]
not_working_pwo = [1824,1835,1891]
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

g = tgraph(module_names.keys(),mean_eff1.values(),dy=err_stat.values())
g.Draw('ap')
v = hcviewer('eff',False,0.95,1.0)
v.fill([primex_to_prad(mean_eff1)])
v.draw()

# efficiency versus theta
load_positions()
hth = TH1F('hth_eff',';theta (deg);efficiency',100,0.,10.)
hthn = TH1F('hth_n',';theta (deg);n',100,0.,10.)
hths = TH1F('hth_s',';theta (deg);sigma',100,0.,10.)
for i in module_names.keys():
  _ = hth.Fill(module_theta[i],mean_eff1[i])
  _ = hthn.Fill(module_theta[i])
  _ = hths.Fill(module_theta[i],1/err_stat[i]**2)
  
hth.Divide(hthn)
for i in range(100):
  if hths[i+1]!=0: hth.SetBinError(i+1,1/hths[i+1]**0.5)

# by energy

zone_stat12 = [[sum(stat12[x][i] for x in y) for i in range(17)] for y in zones]
zone_stat125 = [[sum(stat125[x][i] for x in y) for i in range(17)] for y in zones]
zone_eff = [[a/b if b!=0 else 0 for a,b in zip(x,y)] for x,y in zip(zone_stat12,zone_stat125)]
zone_deff = [[((a+1)*(a+2)/(b+2)/(b+3)-(a+1)**2/(b+2)**2)**0.5 for a,b in zip(x,y)] for x,y in zip(zone_stat12,zone_stat125)]
g2 = [tgraph([0.225+0.050*i for i in range(17)],zone_eff[j],dy=zone_deff[j]) for j in [1,-1]]

h = listofhisto('heff',[range(6),range(17)],2000,0.7,1.1)
for i in module_names.keys():
  for j in range(6):
    if i in zones[j]: 
      for k in range(17): _ = h[j][k].Fill(eff[i][k],1/deff[i][k]**2)


zone_means = [[0. for k in range(17)] for i in range(6)]
zone_sigmas = [[0. for k in range(17)] for i in range(6)]
zone_ent = [[0. for k in range(17)] for i in range(6)]
for i in range(6):
  for k in range(17):
    h[i][k].GetXaxis().SetRangeUser(0.97,1.01)
    h[i][k].Draw()
    if h[i][k].Integral()==0: continue
    s = raw_input()
    r = h[i][k].Fit('gaus','rsl','',0.97,1.01)
    if h[i][k].Integral()==0: continue
    zone_means[i][k] = r.GetParams()[1]
    zone_sigmas[i][k] = r.GetParams()[2]
    zone_ent[i][k] = h[i][k].GetEffectiveEntries()

g = [tgraph([0.225+0.050*i for i in range(17)],zone_means[j],dy=[zone_sigmas[j][k]/(zone_ent[j][k]-1)**0.5 if zone_ent[j][k]>1 else 0.3 for k in range(17)]) for j in [1,-1]]
ffit = TF1('ffit','[0]*(1-exp([1]*x+[2]))',0.2,1.1)
ffit.SetParLimits(0,0.9,1.)
ffit.SetParLimits(1,-50,0)
ffit.SetParLimits(2,-10,10) 



choice_mean = [[0.993,0.993,0.993,0.9936,0.9943,0.9945,0.9946 for k in range(17)] for i in range(2)]
choice_sigma = [[0 for k in range(17)] for i in range(2)]
mean_eff2 = {}
for i in module_names.keys():
  belong = False
  for j in range(6):
    if i in zones[j]:
      belong = True
      mean_eff2[i] = [choice_mean+choice_sigma*(x[i]-zone_means[j])/zone_sigmas[j] for k in range(17)]
      break
  if not belong: mean_eff1[i]=0


# by regions
regions = [get_regions('center_pwo'),get_regions('trans_pwo')+get_regions('trans_lg'),get_regions('center_lg')]
h = listofhisto('heff',[range(3),range(17)],2000,0.7,1.1)
for i in module_names.keys():
  for j in range(3):
    if i in regions[j]: 
      for k in range(17):
        if deff1[i][k]!=0: _ = h[j][k].Fill(eff1[i][k],1/deff1[i][k]**2)

zone_means = [[0. for k in range(17)] for i in range(3)]
zone_sigmas = [[0. for k in range(17)] for i in range(3)]
zone_ent = [[0. for k in range(17)] for i in range(3)]
for i in range(3):
  for k in range(17):
    h[i][k].GetXaxis().SetRangeUser(0.97,1.01)
    h[i][k].Draw()
    if h[i][k].Integral()==0: continue
    s = raw_input()
    r = h[i][k].Fit('gaus','rsl','',0.97,1.01)
    if h[i][k].Integral()==0: continue
    zone_means[i][k] = r.GetParams()[1]
    zone_sigmas[i][k] = r.GetParams()[2]
    zone_ent[i][k] = h[i][k].GetEffectiveEntries()

g = [tgraph([0.225+0.050*i+0.005*j for i in range(17)],zone_means[j],dy=[zone_sigmas[j][k]/(zone_ent[j][k]-1)**0.5 if zone_ent[j][k]>1 else 0.3 for k in range(17)]) for j in range(3)]
