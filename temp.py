cp f = TFile('gains.root')

for x in module_names:
  if x in ['G900','G107','G16','W835']: continue
  g = [f.Get(x).Get('gped_'+x),f.Get(x).Get('gled_'+x),f.Get(x).Get('ggains_'+x+'_LMS1')]
  for i in range(1,g[0].GetN()):
    for j in range(3):
      if abs(g[j].GetY()[i-1]-g[j].GetY()[i])>g[j].GetErrorY(i-1)+g[j].GetErrorY(i):
        print 'problem with:',g[j].GetName(), g[j].GetX()[i]
        g[j].Draw('apl')
        s = raw_input()

for i,run in enumerate(calib_runs['all']):
  t = chao('../replay/prad_'+str(run).zfill(6)+'.dst')
  for j,x in enumerate(module_names):
    if abs(t.ped_sigma[j]-f.Get(x).Get('gped_'+x).GetErrorY(i))/t.ped_sigma[j]>0.1: print 'sigma not compatible for run',run,'and module',x
    if abs(t.ped_mean[j]-f.Get(x).Get('gped_'+x).GetY()[i])>t.ped_sigma[j]: print 'mean not compatible for run',run,'and module',x,':',t.ped_mean[j],f.Get(x).Get('gped_'+x).GetY()[i], t.ped_sigma[j]
  t.close()
    

# test profile of cluster
from readdst import *
from viewer import *
t = chao('../replay/prad_000960.dst')
t.goto(1000002)
t.load('adc')
e = hcevent(t.e,2)
e.clusters


# get bad event numbers from epics
from array import array
runs = calib_runs['all']
l = []
for run in runs:
  print run
  f = TFile('/work/hallb/prad/mlevilla/stability/gains_LMS_'+str(run)+'.root')
  g = f.Get('epics').Get('ghallb_IPM2C24A_CUR')
  n = g.GetN()
  x = [g.GetX()[i] for i in range(n)]
  y = [g.GetY()[i] for i in range(n)]
  i = 1
  lrun = []
  while i<n:
    down = (y[i-1]-y[i])/y[i-1]
    if down>0.2:
      start = i-1
      if i==n-1:
        lrun.append([run,int(x[start])+1,-1])
        break
      i+=1
      per = (y[i]-y[start])/y[start]
      while (per<-0.2 or per>0.1) and i<n-1:
        i+=1
        per = (y[i]-y[start])/y[start]
      if i!=n-1: lrun.append([run,int(x[start])+1,int(x[i])-1])
      else: lrun.append([run,int(x[start])+1,-1])
    i+=1
  l.extend(lrun)
  g2 = [[TGraph(2,array('f',[a[1],a[1]]),array('f',[0.,10.])),TGraph(2,array('f',[a[2],a[2]]),array('f',[0.,10.]))] for a in lrun]
  g.Draw('alp')
  for a in g2: 
    for b in a:
      b.SetLineColor(kRed)
      b.Draw('l')
  _ = raw_input()
      

# length of runs
from readdst import *
runs = calib_runs['all']
l = []
for run in runs:
  print run
  t = chao(run=run,step=100000)
  for ev in progress(t,modulo=10000): pass
  l.append([run,t.ancre])
  print t.ancre

# combine calib_gains files
runs = calib_runs['all']
l = [readfile('calib_gains_'+str(run)+'.txt',[float,int],cols=[1,2]) for run in runs]
l2 = [[] for i in range(1728)]
for x in l:
  for i,y in enumerate(x):
    l2[i].append(y)

a = [module_names[i] for i in range(len(l2)) if l2[i]==[]]
b = [module_names[i] for i in range(len(l2)) if len(l2[i])>1]
c = [sum([y[0]*y[1] for y in x])/sum([y[1] for y in x]) if sum([y[1] for y in x])!=0 else 0.25 for x in l2]
d = [sum([(y[0]-z)*y[1] for y in x])/sum([y[1] for y in x]) if sum([y[1] for y in x])!=0 else 0. for x,z in zip(l2,c)]
e = [sum([y[1] for y in x]) for x in l2]

writelist('gains2.txt',[module_names,c],tr=True)

# group by 2x2 modules
nx_pwo,ny_pwo=2,2
nx_lg,ny_lg=2,2

l = [[30*i+1+j for i in range(a,a+nx_lg) for j in range(b,b+ny_lg)] for a in range(0,30,nx_lg) for b in range(0,30,ny_lg)]+[[1000+34*i+1+j for i in range(a,a+nx_pwo) for j in range(b,b+ny_pwo)] for a in range(0,34,nx_pwo) for b in range(0,34,ny_pwo)]

lid = [[module_names.index(['G','W'][int(x>1000)]+str(x-1000*int(x>1000))) for x in y if ['G','W'][int(x>1000)]+str(x-1000*int(x>1000)) in module_names] for y in l]
while [] in lid: del lid[lid.index([])]


# efficiency
from efficiency_calc import efficiency_calc

# lid = [[i] for i in range(1728)]
runs = sum_dict_lists(runs_list,'calib')
fout = TFile('test.root','recreate')
htrig = [[TH1F('htrig_'+'_'.join([module_names[y] for y in x])+'_'+str(j),'',500,0.,1500) for j in range(6)] for x in lid]
for run in runs:
  f = TFile('efficiency_calib/island_10/efficiency_calib_island_'+str(run)+'.root')
  for i in range(len(lid)):
    for j in range(len(lid[i])):
      for k in range(6): _ = htrig[i][k].Add(f.Get(module_names[lid[i][j]]).Get('htrigger_'+module_names[lid[i][j]]+'_'+str(k)))

h = [efficiency_calc([htrig[i][1],htrig[i][2]],[htrig[i][1],htrig[i][2],htrig[i][5]],name='heff_'+'_'.join([module_names[y] for y in x])) for i,x in enumerate(lid)]
  

# create file for pedestal and lms gains
f = TFile(os.environ['WORK']+'/stability/gains_LMS.root')
lped_mean,lped_sigma,lgain = [],[],[[] for i in range(3)]
for x in module_names:
  gped = f.Get(x).Get('gped_'+x)
  gped_sigma = f.Get(x).Get('gped_sigma_'+x)
  ggain = [f.Get(x).Get('ggains_'+x+'_LMS'+str(i+1)) for i in range(3)]
  n = gped.GetN()
  lped_mean.append([x]+[round(gped.GetY()[i],2) for i in range(n)])
  lped_sigma.append([x]+[round(gped_sigma.GetY()[i],2) for i in range(n)])
  for j in range(3): lgain[j].append([x]+[round(ggain[j].GetY()[i],2) for i in range(n)])

lheader = ['run']+[int(gped.GetX()[i]) for i in range(n)]
lped_mean = [lheader]+lped_mean
lped_sigma = [lheader]+lped_sigma
writelist('pedestal_mean2.txt',lped_mean)
writelist('pedestal_sigma2.txt',lped_sigma)
for i in range(3): 
  lgain[i] = [lheader]+lgain[i]
  writelist('lms_gains2_'+str(i+1)+'.txt',lgain[i])
  

# fake modules surrounding hycal
module_names_more = ['W561','W562','W595','W596']+['O'+str(i) for i in range(1,125)]
xpos_more = [-1.0384,1.0384,-1.0384,1.0384]+[xpos[module_names.index('G1')]-3.815]+[xpos[module_names.index('G'+str(i+1))] for i in range(30)]+[xpos[module_names.index('G30')]+3.815]+[xpos[module_names.index('G'+str(30*(i+1)))]+3.815 for i in range(30)]+[xpos[module_names.index('G900')]+3.815]+[xpos[module_names.index('G'+str(900-i))] for i in range(30)]+[xpos[module_names.index('G871')]-3.815]+[xpos[module_names.index('G'+str(871-30*i))]-3.815 for i in range(30)]
ypos_more = [1.0375,1.0375,-1.0375,-1.0375]+[ypos[module_names.index('G1')]+3.815]+[ypos[module_names.index('G'+str(i+1))]+3.815 for i in range(30)]+[ypos[module_names.index('G30')]+3.815]+[ypos[module_names.index('G'+str(30*(i+1)))] for i in range(30)]+[ypos[module_names.index('G900')]-3.815]+[ypos[module_names.index('G'+str(900-i))]-3.815 for i in range(30)]+[ypos[module_names.index('G871')]-3.815]+[ypos[module_names.index('G'+str(871-30*i))] for i in range(30)]


# iterate over hycal calib tree
filename = os.environ['WORK']+'/python/calib/calib/island_14/tree_calib14.root'
f = TFile(filename)
t = f.Get('hycal')

hratio = [TH1F('hratio_'+module_names[i],';E_{cluter}/E_{#gamma};',500,0.,2.5) for i in range(1728)]
htrigger = [[TH1F('htrigger_'+module_names[i]+'_'+str(j),'',50,0.,1500) for i in range(1728)] for j in range(6)]
for ev in progress(t,modulo=10000,show=3,n=t.GetEntries()):
  if ev.n_cl!=1: continue
  #if ev.chi2_cl[0]/(ev.nhit_cl[0]+ev.nneigh_cl[0]-2)>10 or ev.status_cl[0]==0: continue
  _ = hratio[ev.mid[0]].Fill(ev.E_cl[0]/ev.E_g)
  _ = htrigger[ev.trigger][ev.mid[0]].Fill(ev.E_g)


  
# iterate over hycal prod tree
filename = os.environ['WORK']+'/python/calib/prod/island_0/tree_calib0.root'
f = TFile(filename)
t = f.Get('hycal')

hratio = [[TH1F('hratio_'+module_names[i]+'_'+str(j),';E_{cluter}/E_{#gamma};',500,0.,2.5) for i in range(1728)] for j in range(2)]
hx = TH1F('hx',';x_{center};',2000,-2.,2.)
hx = TH1F('hy',';y_{center};',2000,-2.,2.)
for ev in progress(t,modulo=10000,show=3,n=t.GetEntries()):
  i = 0
  while (i<ev.n_cl):
    if ev.ep_cl[i]==1:
      if ev.chi2_cl[i]/(ev.nhit_cl[i]+ev.nneigh_cl[i]-2)>10: 
        i+=1
        continue
      _ = hratio[ev.mid[i]].Fill(ev.E_cl[i]/ev.ebeam)
      i+=1
    else:
      if ev.chi2_cl[i]/(ev.nhit_cl[i]+ev.nneigh_cl[i]-2)>10 or ev.chi2_cl[i+1]/(ev.nhit_cl[i+1]+ev.nneigh_cl[i+1]-2)>10: 
        i+=2
        continue
      d,d1,d2,x,y = distance_center_cls(x=[ev.x_cl[i],ev.x_cl[i+1]],y=[ev.y_cl[i],ev.y_cl[i+1]])
      _ = hratio[0][ev.mid[i]].Fill(ev.E_cl[i]/(ev.ebeam*d2/(d1+d2)))
      _ = hratio[0][ev.mid[i+1]].Fill(ev.E_cl[i]/(ev.ebeam*d1/(d1+d2)))
      _ = hx.Fill(x)
      _ = hy.Fill(y)
      i+=2

# group by 2x2 modules
nx_pwo,ny_pwo=2,2
nx_lg,ny_lg=2,2

l = [[30*i+1+j for i in range(a,a+nx_lg) for j in range(b,b+ny_lg)] for a in range(0,30,nx_lg) for b in range(0,30,ny_lg)]+[[1000+34*i+1+j for i in range(a,a+nx_pwo) for j in range(b,b+ny_pwo)] for a in range(0,34,nx_pwo) for b in range(0,34,ny_pwo)]

lid = [[module_names.index(['G','W'][int(x>1000)]+str(x-1000*int(x>1000))) for x in y if ['G','W'][int(x>1000)]+str(x-1000*int(x>1000)) in module_names] for y in l]
while [] in lid: del lid[lid.index([])]


# efficiency
from efficiency_calc import efficiency_calc

# lid = [[i] for i in range(1728)]
runs = sum_dict_lists(runs_list,'calib')
htrig = [[TH1F('htrig_'+'_'.join([module_names[y] for y in x])+'_'+str(j),'',50,0.,1500) for j in range(6)] for x in lid]
for i in range(len(lid)):
  for j in range(len(lid[i])):
    for k in range(6): _ = htrig[i][k].Add(htrigger[k][lid[i][j]])

h = [efficiency_calc([htrig[i][1],htrig[i][2]],[htrig[i][1],htrig[i][2],htrig[i][5]],name='heff_'+'_'.join([module_names[y] for y in x])) for i,x in enumerate(lid)]


# create calibration file for PRadEventViewer
ipass = 11
indir = os.environ['WORK']+'/calib_results/calib/pass'+str(ipass)
gains = readfile(indir+'/results.txt',[float],cols=[1],start=1)
gains = [x if 0<x<1 else 0.25 for x in gains]
alpha = readfile(indir+'/results.txt',[float],cols=[2],start=1)
f = TFile(work_folder+'/stability/gains_LMS.root')
run_start,run_end = 1000,1600
for x in periods_calib+periods_prod:
  run_min,run_max = x[0],x[-1]
  if run_max<run_start or run_min>run_end: continue
  print run_min,run_max
  lms = []
  for z in module_names:
    ggain = [f.Get(z).Get('ggains_'+z+'_LMS'+str(y)) for y in range(1,4)]
    j = [ggain[0].GetX()[i] for i in range(ggain[0].GetN())].index(x[0])
    l = [ggain[i].GetY()[j] for i in range(3)]
    # n = ggain[0].GetN()
    # l = [0.,0.,0.]
    # m = 0
    # for i in range(n):
    #   if run_min<=ggain[0].GetX()[i]<=run_max:
    #     m+=1
    #     for j in range(3): l[j]+=ggain[j].GetY()[i]
    # for j in range(3): l[j]/=m
    lms.append([round(a,2) for a in l])
  lf = [[x,'%.4f'%y,'%.4f'%b]+z for x,y,b,z in zip(module_names,gains,alpha,lms)]
  writelist(os.environ['WORK']+'/calib_files/calib/pass'+str(ipass+1)+'/calib_'+str(run_min)+'-'+str(run_max)+'.txt',lf,ljust=10)

ipass = 0
indir = os.environ['WORK']+'/calib_results/prod/pass'+str(ipass)
gains = [readfile(os.environ['WORK']+'/calib_files/prod/pass{0}/calib_{1}-{2}.txt'.format(ipass,x[0],x[-1]),[float],cols=[1]) for x in periods_prod]
peaks = [readfile(os.environ['WORK']+'/calib_results/prod/pass{0}/results2_{1}-{2}.txt'.format(ipass,x[0],x[-1]),[float],cols=[2]) for x in periods_prod]
ebeam = [1100.44,2142.00] 
theta = [atan2((xpos[i]**2+ypos[i]**2)**0.5,582.) for i in range(1728)]
alpha = [[(ebeam[i]-m_e)/(ebeam[i]+m_e)*cos(x)**2 for x in theta] for i in range(2)]
e_moller = [[m_e*(1+y)/(1-y) for y in x] for x in alpha]
for i,x in enumerate(periods_prod):
  for j in range(1728):
    if peaks[i][j] in [0.,-1.]: continue
    elif peaks[i][j]<50.: gains[i][j] = gains[i][j]/(peaks[i][j]/e_moller[int(i>=4)][j])**0.5
    else: gains[i][j] = gains[i][j]/(peaks[i][j]/e_moller[int(i>=4)][j])

f = TFile(work_folder+'/stability/gains_LMS.root')
for k,x in enumerate(periods_prod):
  run_min,run_max = x[0],x[-1]
  print run_min,run_max
  lms = []
  for x in module_names:
    ggain = [f.Get(x).Get('ggains_'+x+'_LMS'+str(y)) for y in range(1,4)]
    n = ggain[0].GetN()
    l = [0.,0.,0.]
    m = 0
    for i in range(n):
      if run_min<=ggain[0].GetX()[i]<=run_max:
        m+=1
        for j in range(3): l[j]+=ggain[j].GetY()[i]
    for j in range(3): l[j]/=m
    lms.append([round(x,2) for x in l])
  lf = [[x,'%.4f'%y,'%.4f'%b]+z for x,y,b,z in zip(module_names,gains[k],[0. for i in range(1728)],lms)]
  writelist(os.environ['WORK']+'/calib_files/prod/pass'+str(ipass+1)+'/calib_'+str(run_min)+'-'+str(run_max)+'.txt',lf,ljust=10)

# create pedestal and lms files for PRadEventViewer
runs = sum_dict_lists(runs_list)
ped_mean = readfile(home_folder+'/pedestal_mean2.txt',[float]*len(runs),cols=range(1,len(runs)+1),header=[int]*len(runs),tr=True,out=dict)
ped_sigma = readfile(home_folder+'/pedestal_sigma2.txt',[float]*len(runs),cols=range(1,len(runs)+1),header=[int]*len(runs),tr=True,out=dict)
lms_gains = [readfile(home_folder+'/lms_gains2_'+str(i+1)+'.txt',[float]*len(runs),cols=range(1,len(runs)+1),header=[int]*len(runs),tr=True,out=dict) for i in range(3)]
template = readfile(os.environ['WORK']+'/PRadEventViewer/config/pedestal.dat',[int,int,int,float,float])
for x in runs:
  if not (980<x<1107): continue
  print x 
  lped = [a[:3]+[b,c] for a,b,c in zip(template,ped_mean[x],ped_sigma[x])]+template[1728:1735]
  writelist(os.environ['WORK']+'/pedestal_files/pedestal_'+str(x)+'.dat',lped,ljust=7)
  lgain = [[a,b,c,d] for a,b,c,d in zip(module_names,lms_gains[0][x],lms_gains[1][x],lms_gains[2][x])]
  writelist(os.environ['WORK']+'/lms_files/lms_gains_'+str(x)+'.txt',lgain,ljust=8)

# create pedestal and lms file for primana
for x in runs:
  if not (980<x<1107): continue
  print x 
  writelist(os.environ['WORK']+'/pedestal_files_primana/mean_'+str(x)+'.txt',prad_to_primex(ped_mean[x]))
  writelist(os.environ['WORK']+'/pedestal_files_primana/sigma_'+str(x)+'.txt',prad_to_primex(ped_sigma[x]))
  for i in range(3):
    writelist(os.environ['WORK']+'/lms_files_primana/lms_gains'+str(i+1)+'_'+str(x)+'.txt',prad_to_primex(lms_gains[i][x]))


# plot pedestal, lms and lms_gains in pdf

f = TFile(os.environ['WORK']+'/python/stability/gains_LMS.root')
c = TCanvas('c','',1000,1200)
c.Divide(1,4)
pref = 'gled_'
suff = ''
outfile = 'lms_graphs.pdf' 
j = 1
for i in range(1,901):
  name = 'G'+str(i)
  if name not in module_names: continue
  g = f.Get(name).Get(pref+name+suff)
  g.SetMarkerStyle(22)
  g.SetMarkerColor(kRed)
  g.SetMarkerSize(0.8)
  _ = c.cd(j)
  g.Draw('ap')
  if j<4 and i!=900: j+=1
  else: 
    j = 1
    if i==4: c.Print(outfile+'(','pdf')
    else: c.Print(outfile,'pdf')

for i in range(1,1157):
  name = 'W'+str(i)
  if name not in module_names: continue
  g = f.Get(name).Get(pref+name+suff)
  g.SetMarkerStyle(22)
  g.SetMarkerColor(kRed)
  g.SetMarkerSize(0.8)
  _ = c.cd(j)
  g.Draw('ap')
  if j<4 and i!=1156: j+=1
  else: 
    j = 1
    if i==1156: c.Print(outfile+')','pdf')
    else: c.Print(outfile,'pdf')


# get ypos-y_cl for each calibration run
ipass = 3
runs = sum_dict_lists(runs_list,'calib')
ldiffx = []
for i,x in enumerate(runs):
  f = TFile(os.environ['WORK']+'/calib_results/calib/pass'+str(ipass)+'/calib_'+str(x)+'.root')
  t = f.Get('hycal')
  t.Draw('xpos-x_cl>>htemp_'+str(x)+'(10000,-700,-500)','E_cl>300 && (trigger==1 || trigger==2)')
  htemp = gROOT.FindObject('htemp_'+str(x))
  m = htemp.GetBinCenter(htemp.GetMaximumBin())
  ldiffx.append(m)


# produce final alpha and gain files
outdir = os.environ['WORK']+'/../config_calibration/maxime'
ipass = 9
alpha = readfile(os.environ['WORK']+'/calib_results/calib/pass'+str(ipass)+'/results.txt',[float],cols=[2],start=1)
gains = readfile(os.environ['WORK']+'/calib_results/calib/pass'+str(ipass)+'/results.txt',[float],cols=[1],start=1)
writelist(outdir+'/non_linearity_snake.txt',[[x,round(y,3)] for x,y in zip(module_names_primex,prad_to_primex(alpha,0))],ljust=5)
writelist(outdir+'/gain_factors_snake.txt',[[x,round(y,3)] for x,y in zip(module_names_primex,prad_to_primex(gains))],ljust=5)

# efficiency and stuff
g = []
nbin = 17
# colors = [kRed-9,kRed-3,kRed,kRed+2,kBlue-9,kBlue-3,kBlue,kBlue+2,kGreen-7,kGreen+1, kGreen+3]
colors = [kMagenta,kBlue,kGreen+1,kRed]
lnames = ['sigma10','sigma5','sigma3', 'sigma2'] #,'sigma_e25','sigma_e50','sigma_e100','sigma_e200','sigma_h2','sigma_h3','sigma_h4']
for x in lnames: 
  folder = os.environ['WORK']+'/eff_results/eff_'+x
  l = [readfile(folder+'/stats&eff2.txt',[float]*nbin,cols=range(1+i*nbin,1+(i+1)*nbin)) for i in range(3)]
  
  l2 = [[[l[i][j][k] for i in range(3)] for k in range(nbin)] for j in range(1728)]
  eff1 = [[y[0]/(y[0]+y[2]) if y[0]+y[2]!=0. else 0. for y in x] for x in l2]
  eff2 = [[y[1]/(y[1]+y[2]) if y[1]+y[2]!=0. else 0. for y in x] for x in l2]
  eff12 = [[(y[0]+y[1])/(y[0]+y[1]+y[2]) if y[0]+y[1]+y[2]!=0. else 0. for y in x] for x in l2]
  eff1_err = [[1./(y[0]+y[2])*(y[0]*(1-y[0]/(y[0]+y[2])))**0.5 if y[0]+y[2]!=0. else 0. for y in x] for x in l2]
  eff2_err = [[1./(y[1]+y[2])*(y[1]*(1-y[1]/(y[1]+y[2])))**0.5 if y[1]+y[2]!=0. else 0. for y in x] for x in l2]
  eff12_err = [[1./(y[0]+y[1]+y[2])*((y[0]+y[1])*(1-(y[0]+y[1])/(y[0]+y[1]+y[2])))**0.5 if y[0]+y[1]+y[2]!=0. else 0. for y in x] for x in l2]
  
  l3 = [[[sum([l[i][j][k] for j in range(1728) if module_names[j][0]==['G','W'][a]]) for i in range(3)] for k in range(nbin)] for a in range(2)]
  eff1_edep = [[y[0]/(y[0]+y[2]) if y[0]+y[2]!=0. else 0. for y in x] for x in l3]
  eff2_edep = [[y[1]/(y[1]+y[2]) if y[1]+y[2]!=0. else 0. for y in x] for x in l3]
  eff12_edep = [[(y[0]+y[1])/(y[0]+y[1]+y[2]) if y[0]+y[1]+y[2]!=0. else 0. for y in x] for x in l3]
  eff1_err_edep = [[1./(y[0]+y[2])*(y[0]*(1-y[0]/(y[0]+y[2])))**0.5 if y[0]+y[2]!=0. else 0. for y in x] for x in l3]
  eff2_err_edep = [[1./(y[1]+y[2])*(y[1]*(1-y[1]/(y[1]+y[2])))**0.5 if y[1]+y[2]!=0. else 0. for y in x] for x in l3]
  eff12_err_edep = [[1./(y[0]+y[1]+y[2])*((y[0]+y[1])*(1-(y[0]+y[1])/(y[0]+y[1]+y[2])))**0.5 if y[0]+y[1]+y[2]!=0. else 0. for y in x] for x in l3]
  
  g.append([TGraphErrors(nbin,array('f',[225.+50*i for i in range(nbin)]),array('f',eff12_edep[j]),array('f',[0. for i in range(nbin)]), array('f',eff12_err_edep[j])) for j in range(2)])
  for i in range(2):
    g[-1][i].SetTitle('')
    g[-1][i].GetXaxis().SetTitle('E_{#gamma} (MeV)')
    g[-1][i].GetXaxis().SetTitleSize(0.07)
    g[-1][i].GetXaxis().SetTitleOffset(0.9)
    g[-1][i].GetXaxis().SetRangeUser(320,1080)
    g[-1][i].GetXaxis().SetNdivisions(505)
    g[-1][i].GetXaxis().SetLabelSize(0.05)
    g[-1][i].GetYaxis().SetTitle('Trigger Efficiency')
    g[-1][i].GetYaxis().SetTitleSize(0.07)
    g[-1][i].GetYaxis().SetTitleOffset(1.0)
    g[-1][i].GetYaxis().SetRangeUser(0.98,0.999)
    g[-1][i].GetYaxis().SetNdivisions(505)
    g[-1][i].GetYaxis().SetLabelSize(0.05)
    g[-1][i].SetMarkerStyle([20,23,22,21][len(g)-1])
    g[-1][i].SetMarkerSize(1.2)
    g[-1][i].SetMarkerColor(colors[len(g)-1])
    g[-1][i].SetLineColor(colors[len(g)-1])
    g[-1][i].SetLineWidth(2)
    g[-1][i].SetFillColor(colors[len(g)-1])
    g[-1][i].SetFillStyle(0)


g[0][1].Draw('ap')
for i in range(1,len(g)): g[i][1].Draw('p')

g0 = TGraph(2,array('f',[450,450]),array('f',[0,2]))
g0.SetLineWidth(3)
g0.SetLineStyle(2)
g0.Draw('l')

t = TText(.0,.3,"Preliminary")
t.SetTextAlign(22)
t.SetTextFont(43)
t.SetTextSize(40)
t.SetTextAngle(45)
t.Draw()

# n_{hits} > 
# #chi^{2} <
# |E/E_{#gamma} - <E/E_{#gamma}>| < #sigma

t = TLatex(0.1,0.3,"Preliminary")
t.SetTextAngle(45)
t.Draw()

from viewer import *
v = hcviewer('eff',false,0.96,1.01)
v.fill([[x[3] for x in eff12]])
v.draw()

# deal with particular cases
ipass = 0
l = [readfile(os.environ['WORK']+'/calib_results/prod/pass{0}/errors_{1}-{2}.txt'.format(ipass,x[0],x[-1]),[int]) for x in periods_prod]
l2 = [readfile(os.environ['WORK']+'/calib_results/prod/pass{0}/results_{1}-{2}.txt'.format(ipass,x[0],x[-1]),[str,float,float,float]) for x in periods_prod]
f = TFile(os.environ['WORK']+'/calib_results/prod/pass{0}/results.root'.format(ipass))
s = TSpectrum(1)
m_e = 0.51099893
ebeam = [1100.44,2142.00] 
theta = [atan2((xpos[i]**2+ypos[i]**2)**0.5,582.) for i in range(1728)]
alpha = [[(ebeam[i]-m_e)/(ebeam[i]+m_e)*cos(x)**2 for x in theta] for i in range(2)]
e_moller = [[m_e*(1+y)/(1-y) for y in x] for x in alpha]
print [len(x) for x in l]
for i,x in enumerate(l):
  per = periods_prod[i]
  gev = int(per[0]>=1362)
  print i,len(x)
  for j in x:
    name = module_names[j]
    h = f.Get(str(per[0])+'-'+str(per[-1])).Get('he_{0}_calib_{1}'.format(name,i))
    h.Draw()
    print x.index(j),per[0],'-',per[-1],name
    st = raw_input()
    _ = s.Search(h)
    if s.GetNPeaks()!=0 and st=='': 
      p = s.GetPositionX()[0]
      r = p/e_moller[gev][j]
      g = l2[i][j][1]/r**0.5
    else: 
      print 'no peak found for {0} in period {1}'.format(name,(i,per))
      p = -1.
      r = -1.
      g = 0.25
    print g,p,r
    print
    l2[i][j] = [name,g,p,r]

    
l = []  
for i,x in enumerate(module_names):
 print i,x
 h = f.Get(x).Get('he_'+x)
 if x[0]=='G': h.GetXaxis().SetRangeUser(40.,500.)
 elif x[0]=='W': h.GetXaxis().SetRangeUser(60.,2000.)
 _ = s.Search(h)
 h.Draw()
 st = raw_input() 
 if st=='': l.append(s.GetPositionX()[0])
 else:
  h.GetXaxis().SetRangeUser(0.,500.)
  _ = s.Search(h)
  h.Draw()
  st = raw_input()
  if st=='': l.append(s.GetPositionX()[0])
  else: 
   _ = s.Search(h)
   h.Draw()
   st = raw_input()
   if st=='': l.append(s.GetPositionX()[0])
   else: l.append(0.)


lbad = ['G130','G131','G132','G133','G134','G135','G783','G784','G630','G16','G107','W891','W835','W527','W528','W560','W563','W594','W597','W629','W630']


# produce calibration files for production periods
ipass = 0
for x in periods_prod:
  in1 = os.environ['WORK']+'/calib_files/prod/pass{0}/calib_{1}-{2}.txt'.format(ipass,x[0],x[-1])
  in2 = os.environ['WORK']+'/calib_results/prod/pass{0}/results_{1}-{2}.txt'.format(ipass,x[0],x[-1])
  l1 = readfile(in1,[str,float,float,float,float,float])
  l2 = readfile(in2,[str,float,float,float,float,float,float])
  l3 = []
  for i in range(1728):
    if abs(l2[i][2]-l2[i][6])/l2[i][2]<0.2: l3.append([l1[i][0],l1[i][1]/l2[i][6]**0.5]+l1[i][2:])
    elif abs(l2[i][2]-l2[i][4])/l2[i][2]<0.2: l3.append([l1[i][0],l1[i][1]/l2[i][4]**0.5]+l1[i][2:])
    elif 0.2<l2[i][2]<2.5: l3.append([l1[i][0],l1[i][1]/l2[i][2]**0.5]+l1[i][2:])
    else:
      print i,module_names[i]
      l3.append(l1[i])
  writelist(os.environ['WORK']+'/calib_files/prod/pass{0}/calib_{1}-{2}.txt'.format(ipass+1,x[0],x[-1]),l3) 


f = TFile('/work/hallb/prad/mlevilla/stability/gains_LMS.root')
lmod = {}
for x in module_names:
  gped_mean = f.Get(x).Get('gped_'+x)
  gped_sigma = f.Get(x).Get('gped_sigma_'+x)
  gled = f.Get(x).Get('gled_'+x)
  if x == 'G176': runs = [int(gled.GetX()[i]) for i in range(gled.GetN())]
  for i,y in enumerate(runs):
    if gled.GetY()[i] < gped_mean.GetY()[i]+10*gped_sigma.GetY()[i] or gled.GetY()[i]>8000: 
      if x in lmod: lmod[x].append(y)
      else: lmod[x] = [y]
  

runs = sum_dict_lists(runs_list)
ped_mean = readfile(home_folder+'/pedestal_mean.txt',[float]*len(runs),cols=range(1,len(runs)+1),header=[int]*len(runs),tr=True,out=dict)
ped_sigma = readfile(home_folder+'/pedestal_sigma.txt',[float]*len(runs),cols=range(1,len(runs)+1),header=[int]*len(runs),tr=True,out=dict)
lmod = []
for x in l:
  if x in [1136,1243,1342,1362,1483,1515,1516]: 
    lmod.append([])
    continue
  print x
  lrun = []
  f = TFile('/work/hallb/prad/mlevilla/adc_histos/adc_peak_'+str(x)+'.root')
  for i in range(1728):
    name = module_names[i]
    h = f.Get('hadc_maxe_'+name)
    if h.GetEntries()<10 or h.GetMean()<ped_mean[x][i]+10*ped_sigma[x][i]: lrun.append(i)
  _ = gROOT.GetListOfFiles().Remove(f)
  lmod.append(lrun)


f = TFile(os.environ['WORK']+'/calib_results/calib/pass11/results2.root')
g = f.pwo.Get('gsigma_pwo')
hs = [TH1F('hs_'+str(i),'',100,0.02,0.05) for i in range(6)]
for x in module_names:
  if x[0]=='G': continue
  gmod = f.Get(x).Get('gsigma_'+x)
  xmod = [gmod.GetX()[i] for i in range(gmod.GetN())]
  for i in range(6):
    if 350+100*i not in xmod: continue
    _ = hs[i].Fill(gmod.GetY()[xmod.index(350+100*i)])

sigma = []
for i in range(6):
  peak = hs[i].GetBinCenter(hs[i].GetMaximumBin())
  hs[i].Fit('gaus','r','',peak-0.002,peak+0.002)
  sigma.append(hs[i].GetFunction('gaus').GetParameter(2))

g2 = TGraphErrors(6,array('f',[g.GetX()[i] for i in range(6)]),array('f',[g.GetY()[i] for i in range(6)]),array('f',[0. for i in range(6)]), array('f',[(g.GetErrorY(i)**2+sigma[i]**2)**0.5 for i in range(6)]))


# efficiency from tree
folders = [os.environ['WORK']+'/calib_results/calib/pass12/',os.environ['WORK']+'/calib_results/fsnake/pass0/']
lnames = ['_0','_t40','_t40_xy20','_t40_xy20_s5','_t40_xy20_s4','_t40_xy20_s3','_t40_xy20_s2']
files = [[x+'stats'+y+'.txt' for y in lnames] for x in folders]
files2 = [[x+'stats_cleaned'+y+'.txt' for y in lnames] for x in folders]
nbin = 17
l = [[[readfile(y,[float]*nbin,cols=range(1+i*nbin,1+(i+1)*nbin)) for i in range(3)] for y in x] for x in files]
l2 = [[[readfile(y,[float]*nbin,cols=range(1+i*nbin,1+(i+1)*nbin)) for i in range(3)] for y in x] for x in files2]
eff = [[[[(y[0][i][j]+y[1][i][j])/(y[0][i][j]+y[1][i][j]+y[2][i][j]) if y[0][i][j]+y[1][i][j]+y[2][i][j]!=0. else 0. for j in range(nbin)] for i in range(1733)] for y in x] for x in l]
eff2 = [[[[(y[0][i][j]+y[1][i][j])/(y[0][i][j]+y[1][i][j]+y[2][i][j]) if y[0][i][j]+y[1][i][j]+y[2][i][j]!=0. else 0. for j in range(nbin)] for i in range(1733)] for y in x] for x in l2]
stat = [[[[1./(y[0][i][j]+y[1][i][j]+y[2][i][j])*((y[0][i][j]+y[1][i][j])*(1-(y[0][i][j]+y[1][i][j])/(y[0][i][j]+y[1][i][j]+y[2][i][j])))**0.5 if y[0][i][j]+y[1][i][j]+y[2][i][j]!=0. else 0. for j in range(nbin)] for i in range(1733)] for y in x] for x in l]
stat2 = [[[[1./(y[0][i][j]+y[1][i][j]+y[2][i][j])*((y[0][i][j]+y[1][i][j])*(1-(y[0][i][j]+y[1][i][j])/(y[0][i][j]+y[1][i][j]+y[2][i][j])))**0.5 if y[0][i][j]+y[1][i][j]+y[2][i][j]!=0. else 0. for j in range(nbin)] for i in range(1733)] for y in x] for x in l2]
h = [[[[TH1F('h'+str(i)+str(j)+str(k)+str(m),';eff_{mod}',2000,0.,1.) for m in range(nbin)] for k in range(5)] for j,y in enumerate(x)] for i,x in enumerate(eff)]
h2 = [[[[TH1F('h2'+str(i)+str(j)+str(k)+str(m),';eff_{mod}',2000,0.,1.) for m in range(nbin)] for k in range(5)] for j,y in enumerate(x)] for i,x in enumerate(eff)]

for i in range(len(eff)):
  for j in range(len(eff[i])):
    for k,x in enumerate(module_names):
      for m in range(nbin):
        b = [x[0]=='G',x[0]=='W',x in names_trans_lg,x in names_trans_pwo,x in names_center]
        for n in range(5):
          if b[n]:
            _ = h[i][j][n][m].Fill(eff[i][j][k][m])
            _ = h2[i][j][n][m].Fill(eff[i][j][k][m])

syst = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in eff]
syst2 = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in eff]
mean = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in eff]
mean2 = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in eff]
for i in range(len(eff)):
  for j in range(len(eff[i])):
    for k in range(5):
      for m in range(nbin):
        print i,j,k,m
        ffit = TF1('ffit','gaus',0.6,1.)
        h[i][j][k][m].GetXaxis().SetRangeUser(0.6,1.)
        peak = h[i][j][k][m].GetBinCenter(h[i][j][k][m].GetMaximumBin())
        rms = h[i][j][k][m].GetStdDev()
        intg = h[i][j][k][m].Integral()
        print i,j,k,m,intg
        if intg>150:
          ffit.SetParameter(1,peak)
          ffit.SetParLimits(1,peak-rms,peak+rms)
          ffit.SetParLimits(2,0.5*rms,rms)
          _ = h[i][j][k][m].Fit(ffit,'qr','',peak-rms,peak+rms)
          syst[i][j][k][m] = ffit.GetParameter(2)  if ffit.GetParameter(1)<1 else rms
          mean[i][j][k][m] = ffit.GetParameter(1)  if ffit.GetParameter(1)<1 else peak
        else: 
          syst[i][j][k][m] = rms
          mean[i][j][k][m] = peak
        h2[i][j][k][m].GetXaxis().SetRangeUser(0.6,1.)
        peak = h2[i][j][k][m].GetBinCenter(h2[i][j][k][m].GetMaximumBin())
        rms = h2[i][j][k][m].GetStdDev()
        intg = h2[i][j][k][m].Integral()
        print i,j,k,m,intg
        if intg>150:
          ffit.SetParameter(1,peak)
          ffit.SetParLimits(1,peak-rms,peak+rms)
          ffit.SetParLimits(2,0.5*rms,rms)
          _ = h2[i][j][k][m].Fit(ffit,'qr','',peak-rms,peak+rms)
          syst2[i][j][k][m] = ffit.GetParameter(2) if ffit.GetParameter(1)<1 else rms
          mean2[i][j][k][m] = ffit.GetParameter(1) if ffit.GetParameter(1)<1 else peak
        else: 
          syst2[i][j][k][m] = rms
          mean2[i][j][k][m] = peak
        
err = [[[[(stat[i][j][1728+k][m]**2+syst[i][j][k][m]**2)**0.5 for m in range(nbin)] for k in range(5)] for j in range(len(stat[i]))] for i in range(len(stat))]
err2 = [[[[(stat2[i][j][1728+k][m]**2+syst2[i][j][k][m]**2)**0.5 for m in range(nbin)] for k in range(5)] for j in range(len(stat2[i]))] for i in range(len(stat2))]

g = [[[tgraph([225+50*i+2*j for i in range(nbin)],y,[0. for i in range(nbin)],yerr) for y,yerr in zip(a,aerr)] for j,(a,aerr) in enumerate(zip(b,berr))] for b,berr in zip(mean,err)]
g2 = [[[tgraph([225+50*i+2*j for i in range(nbin)],y,[0. for i in range(nbin)],yerr) for y,yerr in zip(a,aerr)] for j,(a,aerr) in enumerate(zip(b,berr))] for b,berr in zip(mean2,err2)]
colors = [kBlack,kBlue,kGreen-6,kGreen+2,kMagenta,kRed-7,kRed]
i,k = 0,0
c = custom(g2[i][0][k],color=colors[0],first=True,xtitle='E_{#gamma} (MeV)',ytitle='Trigger Efficiency',grid=True)
for j in range(4,7): custom(g2[i][j][k],c,color=colors[j])   


#resolution from tree
folders = [os.environ['WORK']+'/calib_results/calib/pass12/',os.environ['WORK']+'/calib_results/fsnake/pass0/']
#lnames = ['_0','_t40','_t40_xy20','_t40_xy20_s5','_t40_xy20_s4']
lnames = ['_t40_xy20_s4']
files = [[x+'sigma_tr2'+y+'.txt' for y in lnames] for x in folders]
dfiles = [[x+'dsigma_tr2'+y+'.txt' for y in lnames] for x in folders]
nbin = 9
l = [[readfile(y,[float]*nbin,cols=range(1,1+nbin)) for y in x] for x in files]
dl = [[readfile(y,[float]*nbin,cols=range(1,1+nbin)) for y in x] for x in dfiles]
sigma = [[[[y[i][j] for j in range(nbin)] for i in range(1733)] for y in x] for x in l]
stat = [[[[y[i][j] for j in range(nbin)] for i in range(1733)] for y in x] for x in dl]
h = [[[[TH1F('h'+str(i)+str(j)+str(k)+str(m),';reso_{mod}',2000,0.02,0.3) for m in range(nbin)] for k in range(5)] for j,y in enumerate(x)] for i,x in enumerate(sigma)]

for i in range(len(sigma)):
  for j in range(len(sigma[i])):
    for k,x in enumerate(module_names):
      for m in range(nbin):
        b = [x[0]=='G',x[0]=='W',x in names_trans_lg,x in names_trans_pwo,x in names_center]
        for n in range(5):
          if b[n]:
            _ = h[i][j][n][m].Fill(sigma[i][j][k][m])

syst = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in sigma]
mean = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in sigma]
for i in range(len(sigma)):
  for j in range(len(sigma[i])):
    for k in range(5):
      for m in range(nbin):
        ffit = TF1('ffit','gaus',0.02,0.3)
        peak = h[i][j][k][m].GetBinCenter(h[i][j][k][m].GetMaximumBin())
        rms = h[i][j][k][m].GetStdDev()
        intg = h[i][j][k][m].Integral()
        print i,j,k,m,intg
        if intg>150:
          ffit.SetParameter(1,peak)
          ffit.SetParLimits(1,peak-rms,peak+rms)
          ffit.SetParLimits(2,0.5*rms,rms)
          _ = h[i][j][k][m].Fit(ffit,'qr','',peak-rms,peak+rms)
          syst[i][j][k][m] = ffit.GetParameter(2) 
          mean[i][j][k][m] = ffit.GetParameter(1) 
        else: 
          syst[i][j][k][m] = rms
          mean[i][j][k][m] = peak

        
err = [[[[(stat[i][j][1728+k][m]**2+syst[i][j][k][m]**2)**0.5 for m in range(nbin)] for k in range(5)] for j in range(len(stat[i]))] for i in range(len(stat))]
errcenter = [[[stat[i][j][1732][m] for m in range(nbin)] for j in range(len(stat[i]))] for i in range(len(stat))]

g = [[[tgraph([250+100*i for i in range(nbin)],y,[0. for i in range(nbin)],yerr) for y,yerr in zip(a,aerr)] for j,(a,aerr) in enumerate(zip(b,berr))] for b,berr in zip(mean,err)]
gcenter = [[tgraph([250+100*i for i in range(nbin)],a[4],[0. for i in range(nbin)],aerr) for (a,aerr) in zip(b,berr)] for b,berr in zip(mean,errcenter)]
colors = [kBlack,kBlue,kGreen+2,kMagenta,kRed]
i,k = 0,0
c = custom(g[i][0][k],color=colors[0],first=True,xtitle='E (MeV)',ytitle='#sigma(E)/E',grid=True)
for j in range(1,5): custom(g[0][j][0],c,color=colors[j])   

#frac{#sigma(E)}{E} = #frac{}{#sqrt{E (GeV)}}

#non linearity from tree
folders = [os.environ['WORK']+'/calib_results/calib/pass13/',os.environ['WORK']+'/calib_results/fsnake/pass1/']
lnames = ['_t40_xy20']
lnames = ['_0','_t40','_t40_xy20','_t40_xy20_s5','_t40_xy20_s4','_t40_xy20_s3','_t40_xy20_s2']
files = [[x+'mean_tr2'+y+'.txt' for y in lnames] for x in folders]
dfiles = [[x+'dmean_tr2'+y+'.txt' for y in lnames] for x in folders]
nbin = 9
l = [[readfile(y,[float]*nbin,cols=range(1,1+nbin)) for y in x] for x in files]
dl = [[readfile(y,[float]*nbin,cols=range(1,1+nbin)) for y in x] for x in dfiles]
means = [[[[y[i][j] for j in range(nbin)] for i in range(1733)] for y in x] for x in l]
stat = [[[[y[i][j] for j in range(nbin)] for i in range(1733)] for y in x] for x in dl]
h = [[[[TH1F('h'+str(i)+str(j)+str(k)+str(m),';means_{mod}',2000,0.5,1.5) for m in range(nbin)] for k in range(5)] for j,y in enumerate(x)] for i,x in enumerate(means)]

for i in range(len(means)):
  for j in range(len(means[i])):
    for k,x in enumerate(module_names):
      for m in range(nbin):
        b = [x[0]=='G',x[0]=='W',x in names_trans_lg,x in names_trans_pwo,x in names_center]
        for n in range(5):
          if b[n]:
            _ = h[i][j][n][m].Fill(means[i][j][k][m])

syst = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in means]
mean = [[[[0. for m in range(nbin)] for k in range(5)] for y in x] for x in means]
for i in range(len(means)):
  for j in range(len(means[i])):
    for k in range(5):
      for m in range(nbin):
        print i,j,k,m
        ffit = TF1('ffit','gaus',0.5,1.5)
        peak = h[i][j][k][m].GetBinCenter(h[i][j][k][m].GetMaximumBin())
        rms = h[i][j][k][m].GetStdDev()
        intg = h[i][j][k][m].Integral()
        print i,j,k,m,intg
        if intg>150:
          ffit.SetParameter(1,peak)
          ffit.SetParLimits(1,peak-rms,peak+rms)
          ffit.SetParLimits(2,0.5*rms,rms)
          _ = h[i][j][k][m].Fit(ffit,'qr','',peak-rms,peak+rms)
          syst[i][j][k][m] = ffit.GetParameter(2) 
          mean[i][j][k][m] = ffit.GetParameter(1) 
        else: 
          syst[i][j][k][m] = rms
          mean[i][j][k][m] = peak
        
err = [[[[(stat[i][j][1728+k][m]**2+syst[i][j][k][m]**2)**0.5 for m in range(nbin)] for k in range(5)] for j in range(len(stat[i]))] for i in range(len(stat))]

g = [[[tgraph([250+100*i for i in range(nbin)],y,[0. for i in range(nbin)],yerr) for y,yerr in zip(a,aerr)] for j,(a,aerr) in enumerate(zip(b,berr))] for b,berr in zip(mean,err)]
colors = [kBlack,kBlue,kGreen+2,kMagenta,kRed]
i,k = 0,0
c = custom(g[i][0][k],color=kBlack,first=True,xtitle='E (MeV)',ytitle='<E/E_{#gamma}>',grid=True)

# do graph from reso.py for special zones
folders = ['/work/hallb/prad/mlevilla/calib_results/'+x for x in ['calib/pass13/','fsnake/pass0/']]
names = ['tr2','tr2_t40','tr2_t40_xy20','tr2_t40_xy20_n3']
sig = [3,2]
region = [get_regions('edge'),get_regions('edge2'),get_regions('lg')-get_regions('edge2')-get_regions('trans_lg'),get_regions('pwo')-get_regions('trans_pwo')-get_regions('inner')]
l = [[[readfile(x+'d_'+y+'/special_sigma_'+y+'_s'+str(z)+'.txt',[float]*9,cols=range(1,10)) for z in sig] for y in names] for x in folders]
dl = [[[readfile(x+'d_'+y+'/special_dsigma_'+y+'_s'+str(z)+'.txt',[float]*9,cols=range(1,10)) for z in sig] for y in names] for x in folders]
lmod = [[[readfile(x+'d_'+y+'/sigma_'+y+'_s'+str(z)+'.txt',[float]*9,cols=range(1,10)) for z in sig] for y in names] for x in folders]

lsyst = [[[[[(sum([(l[ifo][ina][isig][ir][ibin]-lmod[ifo][ina][isig][imod][ibin])**2 for imod in range(1728) if module_names[imod] in region[ir] and lmod[ifo][ina][isig][imod][ibin]!=0])/sum([1 for imod in range(1728) if module_names[imod] in region[ir] and lmod[ifo][ina][isig][imod][ibin]!=0])/len(region[ir])+dl[ifo][ina][isig][ir][ibin]**2)**0.5 for ibin in range(9)] for ir in range(len(region))] for isig in range(len(sig))] for ina in range(len(names))] for ifo in range(len(folders))]

gstat = [[[[tgraph([250+100*i for i,z in enumerate(y) if z!=0 and i not in [0,8]],[z for i,z in enumerate(y) if z!=0 and i not in [0,8]],[0. for i,z in enumerate(y) if z!=0 and i not in [0,8]],[dz for i,(z,dz) in enumerate(zip(y,dy)) if z!=0 and i not in [0,8]]) for y,dy in zip(x,dx)] for x,dx in zip(w,dw)] for w,dw in zip(v,dv)] for v,dv in zip(l,dl)]
gsyst = [[[[tgraph([250+100*i for i,z in enumerate(y) if z!=0 and i not in [0,8]],[z for i,z in enumerate(y) if z!=0 and i not in [0,8]],[0. for i,z in enumerate(y) if z!=0 and i not in [0,8]],[dz for i,(z,dz) in enumerate(zip(y,dy)) if z!=0 and i not in [0,8]]) for y,dy in zip(x,dx)] for x,dx in zip(w,dw)] for w,dw in zip(v,dv)] for v,dv in zip(l,lsyst)]

i,j,k,m = 0,3,1,2

ffit = TF1('ffit','(([0]/(x*0.001)**0.5)**2+([1]/(x*0.001))**2+([2]/(x*0.001)**2)**2)**0.5',300,1000)
ffit2 = TF1('ffit2','(([0]/(x*0.001)**0.5)**2+([1]/(x*0.001))**2+([2]/(x*0.001)**2)**2)**0.5',700,1000)
ffit2.SetLineColor(ROOT.kBlue)
for i in range(len(folders)):
  for j in range(len(names)):
    for k in range(len(sig)):
      for m in range(len(region)):
        gsyst[i][j][k][m].Fit(ffit,'0'); gsyst[i][j][k][m].Fit(ffit2,'0r');
        if m==3: insert = [[0.15,0.15,0.55,0.30,'#frac{#sigma(E)}{E} = #frac{'+str(round(ffit.GetParameter(0),3))+'}{#sqrt{E GeV}}',ffit.GetLineColor()]]
        else:
          insert = []
          if ffit.GetParameter(0)>50*ffit.GetParameter(1):
            insert.append([0.15,0.15,0.55,0.30,'#frac{#sigma(E)}{E} = #frac{'+str(round(ffit.GetParameter(0),3))+'}{#sqrt{E GeV}}',ffit.GetLineColor()])
          else:
            insert.append([0.15,0.15,0.55,0.30,'#frac{#sigma(E)}{E} = #sqrt{#(){#frac{'+str(round(ffit.GetParameter(0),3))+'}{#sqrt{E}}}^{2}+#(){#frac{'+str(round(abs(ffit.GetParameter(1)),3))+'}{E}}^{2}}',ffit.GetLineColor()])
          if ffit2.GetParameter(2)>50*ffit2.GetParameter(1):
            insert.append([0.51,0.76,0.91,0.91,'#frac{#sigma(E)}{E} = #frac{'+str(round(ffit2.GetParameter(2),3))+'}{E^{2} GeV^{2}}',ffit2.GetLineColor()])
          else:
            insert.append([0.51,0.76,0.91,0.91,'#frac{#sigma(E)}{E} = #sqrt{#(){#frac{'+str(round(ffit2.GetParameter(1),3))+'}{E}}^{2}+#(){#frac{'+str(round(ffit2.GetParameter(2),3))+'}{E^{2}}}^{2}}',ffit2.GetLineColor()])
        c,t = custom(gsyst[i][j][k][m],color=ROOT.kBlack,first=True,xtitle='E (MeV)',ytitle='#sigma(E)/E',xrg=[310,990],yrg=[min(l[i][j][k][m])-0.1*(max(l[i][j][k][m])-min(l[i][j][k][m])),max(l[i][j][k][m])+0.3*(max(l[i][j][k][m])-min(l[i][j][k][m]))],grid=True,insert=insert)
        ffit.Draw('same')
        if m!=3: ffit2.Draw('same');
        c.Print('img/161111/reso_special_{0}_{1}_{2}_{3}.pdf'.format(['cpp','primex'][i],names[j],sig[k],['edge','edge2','lg_center','pwo_center'][m]))

# do graph from reso.py for transition regions
folders = ['/work/hallb/prad/mlevilla/calib_results/'+x for x in ['calib/pass13/','fsnake/pass0/']]
names = ['tr2','tr2_t40','tr2_t40_xy20','tr2_t40_xy20_n3']
sig = [3,2]
l = [[[[readfile(x+'d_'+y+'/trans_reso_'+y+'_s'+str(z)+'.txt',[float],cols=[1],rows=range(i*16,(i+1)*16)) for i in range(4)] for z in sig] for y in names] for x in folders]
dl = [[[[readfile(x+'d_'+y+'/trans_reso_'+y+'_s'+str(z)+'.txt',[float],cols=[2],rows=range(i*16,(i+1)*16)) for i in range(4)] for z in sig] for y in names] for x in folders]
g = [[[[tgraph([-3.75+0.5*i for i,z in enumerate(y) if z!=0],[z for z in y if z!=0],[0. for z in y if z!=0],[dz for z,dz in zip(y,dy) if z!=0]) for y,dy in zip(x,dx)] for x,dx in zip(w,dw)] for w,dw in zip(v,dv)] for v,dv in zip(l,dl)]
for i in range(len(folders)):
  for j in range(len(names)):
    for k in range(len(sig)):
      for m in range(4):
        c,t = custom(g[i][j][k][m],color=ROOT.kBlack,first=True,xtitle='Distance from the transition region (cm)',ytitle='#sigma(E)/E',xrg=[-4.5,4.5],yrg=[min(l[i][j][k][m])-0.2*(max(l[i][j][k][m])-min(l[i][j][k][m])),max(l[i][j][k][m])+0.2*(max(l[i][j][k][m])-min(l[i][j][k][m]))],grid=True)
        c.Print('img/161111/reso_trans_{0}_{1}_{2}_{3}.pdf'.format(['cpp','primex'][i],names[j],sig[k],['top','right','bottom','left'][m]))


# draw inner modules
folders = ['/work/hallb/prad/mlevilla/calib_results/'+x for x in ['calib/pass13/','fsnake/pass0/']]
names = ['tr2','tr2_t40','tr2_t40_xy20']
modules = ['W527','W528','W560','W563','W594','W597','W629','W630']
imod = [2,3,5,8,9,12,14,15]
f = [[TFile(x+'/d_'+y+'/histo_'+y+'.root') for y in names] for x in folders]
h = [[[y.Get(z).Get('hr_'+z+'_3') for z in modules] for y in x] for x in f]
c = TCanvas('c','',0,0,800,500)
c.Divide(4,4)
i = 1
for j in range(len(modules)):
  p = c.cd(imod[j])
  p,t = custom(h[i][2][j],c=p,color=ROOT.kBlue,first=True,xtitle='E_{cluster}/E_{#gamma}',ytitle='',xrg=[0.5,2.0],grid=True)

gmodules = ['W527','W563','W594']
l = [readfile(x+'d_tr2_t40_xy20/inner_sigma_tr2_t40_xy20_s2.txt',[str]+[float]*9) for x in folders]
l2 = [[l[i][[y[0] for y in l[i]].index(x)][1:] for x in gmodules] for i in range(2)]
dl = [readfile(x+'d_tr2_t40_xy20/inner_dsigma_tr2_t40_xy20_s2.txt',[str]+[float]*9) for x in folders]
dl2 = [[dl[i][[y[0] for y in dl[i]].index(x)][1:] for x in gmodules] for i in range(2)]
g = [[tgraph([250+100*i for i,z in enumerate(y) if z!=0 and i<4],[z for i,z in enumerate(y) if z!=0 and i<4],[0. for i,z in enumerate(y) if z!=0 and i<4],[dz for i,(z,dz) in enumerate(zip(y,dy)) if z!=0 and i<4]) for y,dy in zip(x,dx)] for x,dx in zip(l2,dl2)]
ffit = TF1('ffit','(([0]/(x*0.001)**0.5)**2+([1]/(x*0.001))**2+([2]/(x*0.001)**2)**2)**0.5',300,1000)
for i in range(2):
  for j in range(3):
    g[i][j].Fit(ffit)
    print i,j ffit.GetParameter(0),ffit.GetParameter(1),ffit.GetParameter(2)
    insert = [[0.15,0.15,0.55,0.30,'#frac{#sigma(E)}{E} = #frac{'+str(round(ffit.GetParameter(0),3))+'}{#sqrt{E GeV}}',ffit.GetLineColor()]]
    c,t = custom(g[i][j],color=ROOT.kBlack,first=True,xtitle='E (MeV)',ytitle='#sigma(E)/E',xrg=[310,990],yrg=[0.02,0.1],grid=True,insert=insert)
    c.Print('img/161111/reso_'+['cpp','primex'][i]+'_'+gmodules[j]+'.pdf')

# draw energy resolution vs distance from cluster center
files = ['/work/hallb/prad/mlevilla/calib_results/'+x+'/pass'+y+'/d_tr2_t40_xy20/inside_reso_tr2_t40_xy20_s2.txt' for x,y in zip(['calib','fsnake'],['13','0'])]
l = [readfile(x,[float,float],cols=[1,2]) for x in files]
g = [[tgraph([3*(2-j/2)*k-22.5*(2-j/2) for k in range(16)],[l[i][k][0] for k in range(16*j,16*(j+1))],dy=[(l[i][k][1]**2+[[6.3e-4,1.7e-4],[6.3e-4,5.7e-5]][i][j/2]**2)**0.5 for k in range(16*j,16*(j+1))]) for j in range(4)] for i in range(2)]


# correcting inner modules gains
mods = ['W526','W527','W528','W529','W560','W563','W594','W597','W628','W629','W630','W631']
f = TFile('/work/hallb/prad/mlevilla/calib_results/calib/pass14/inner.root')
h1 = [f.Get('hr1_'+mods[i]) for i in range(12)]
h2 = [f.Get('hr2_'+mods[i]) for i in range(12)]
linner = []
i = 0
ffit = TF1('ffit','[0]*exp(-[1]*x)+gaus(2)',0.5,5)
ffit.SetParLimits(2,100,300)
ffit.SetParLimits(3,0.9,1.1)
ffit.SetParLimits(4,0.02,0.1)
for i in range(12): 
  h2[i].Fit(ffit,'r')
  h2[i].Fit(ffit,'r')
  _ = raw_input()

linner = [x.GetFunction('ffit').GetParameter(3) for x in h2]
l = readfile('/work/hallb/prad/mlevilla/calib_files/calib/pass14/calib_889-979.txt',[str]+[float]*5)

for i in range(1728):
  if module_names[i] in mods: 
    l[i][1] = round(l[i][1]/linner[mods.index(module_names[i])],3)
  else: 
    l[i][1] = round(l[i][1],3)
                
writelist('/work/hallb/prad/mlevilla/calib_files/calib/pass14/calib_889-979.txt',l,ljust=9)

# changing inner module gains
mods = ['W526','W527','W528','W529','W560','W563','W594','W597','W628','W629','W630','W631']
ipass = [0,13]
method = ['fsnake','calib']
columns = [6,5]
factor = 0.25
for i,m,c in zip(ipass,method,columns):
  l = readfile('/work/hallb/prad/mlevilla/calib_files/'+m+'/pass'+str(i)+'/calib_889-979_original.txt',[str]+[float]*c)
  for x in mods:
    l[module_names.index(x)][1] = l[module_names.index(x)][1]*factor
  writelist('/work/hallb/prad/mlevilla/calib_files/'+m+'/pass'+str(i)+'/calib_889-979.txt',l,ljust=9)

# checking distributions for inner modules with tails
names = ['025','05','10','15','20']
mods = ['W526','W527','W528','W529','W560','W563','W594','W597','W628','W629','W630','W631']
f = [TFile('/work/hallb/prad/mlevilla/calib_results/fsnake/pass0/inner_'+x+'.root') for x in names]
hr0 = [[x.Get(y).Get('hr0_'+y) for y in mods] for x in f]
hr2 = [[x.Get(y).Get('hr2_'+y) for y in mods] for x in f]
mean0 = [[fit_gaussian(hr0[i][j],10.,1.,0.,2.)[0] for i in range(5)] for j in range(12)]
mean2 = [[0. for i in range(5)] for j in range(12)]
ffit = TF1('ffit','(x>[1] ? [0]*exp(-log(x-[1])^2/2/[2]^2) : 0.) + gaus(3)',0.,5.)
ffit.SetParLimits(1,-1.,2.)
ffit.SetParLimits(2,0.01,10.)
ffit.SetParLimits(4,0.2,1.5)
ffit.SetParLimits(5,0.01,0.5)

for i in range(5):
  for j in range(12):
    peak = hr2[i][j].GetMaximum()
    ffit.SetParLimits(0,peak/2.,peak)
    hr2[i][j].Fit(ffit)
    mean2[j][i] = ffit.GetParameter(1)+1
    
# x-xpos and y-ypos sigma and mean for each run
runs =  sum_dict_lists(runs_list,'snake')
ymean,xmean,ysigma,xsigma = [],[],[],[]
for x in runs:
  f = TFile('/work/hallb/prad/mlevilla/calib_results/calib/pass13/trees/tree_'+str(x)+'.root')
  f.hycal.Draw("y-ypos>>htemp(1000,-300,300)","abs(tg-tcl[0])<40 && E>300")
  ymean.append(int(raw_input()))
  ysigma.append(int(raw_input()))


# spatial resolution 'deconvolution'
f0 = TFile('/work/hallb/prad/mlevilla/photon_beam_profile.root')
t0 = f0.tree
nbin0, r0 = 500,50
h0x = [TH1F('h0x'+str(j),'',nbin0,-r0,r0) for j in range(8)]
h0y = [TH1F('h0y'+str(j),'',nbin0,-r0,r0) for j in range(8)]
for j in range(8): 
  t0.Draw('x+5817*tx/sqrt(1-tx^2-ty^2)>>hx'+str(j),'E>0.2+0.1*'+str(j)+' && E<0.3+0.1*'+str(j))
  t0.Draw('y+5817*ty/sqrt(1-tx^2-ty^2)>>hy'+str(j),'E>0.2+0.1*'+str(j)+' && E<0.3+0.1*'+str(j))

def fpy(x,par,j,axis):
  if axis=='x':
    return par[0]*sum([hx[j][int(nbin0/2./r0*(i+r0))+1]*math.exp(-0.5*((x[0]-i-par[1])/par[2])**2) for i in range(-r0,r0)])/sum([hx[j][int(nbin0/2./r0*(i+r0))+1] for i in range(-r0,r0)])  
  elif axis=='y': 
    return par[0]*sum([hy[j][int(nbin0/2./r0*(i+r0))+1]*math.exp(-0.5*((x[0]-i-par[1])/par[2])**2) for i in range(-r0,r0)])/sum([hy[j][int(nbin0/2./r0*(i+r0))+1] for i in range(-r0,r0)])

[[14.3, 5.8], [14.1, 5.7], [14.0, 5.6], [13.8, 5.6], [13.7, 5.6], [13.6, 5.5], [13.5, 5.5], [13.3, 5.5]]

fconv = TF1Convolution('exp(-0.5*(x-[0])^2/13.8^2)/(1+x^2/5.65^2)','[0]*exp(-0.5*(x-[1])^2/[2]^2)',-40,40)
fconv2 = TF1('fconv2',fconv,-40,40,3)
def fpy(x,par): return fconv2.EvalPar(x,par)*(1+erf(par[3]*(x[0]-par[1])))
fconv3 = TF1('fconv3',fpy,-40,40,4)

run = 966
f1 = TFile('/work/hallb/prad/mlevilla/calib_results/calib/pass13/trees/tree_'+str(run)+'.root')
t1 = f1.hycal
nbin1,r1 = 200,40
lcor = readfile('/work/hallb/prad/mlevilla/transporter_time_position/transporter_time_position_'+str(run)+'.txt',[long,float,float])
h1x = [TH1F('h1x'+str(j),'',nbin1,-r1,r1) for j in range(8)]
h1y = [TH1F('h1y'+str(j),'',nbin1,-r1,r1) for j in range(8)]
# lcor968 = readfile('/work/hallb/prad/mlevilla/transporter_time_position/transporter_time_position_968.txt',[long,float,float])
# lcor969 = readfile('/work/hallb/prad/mlevilla/transporter_time_position/transporter_time_position_969.txt',[long,float,float])
# lcor970 = readfile('/work/hallb/prad/mlevilla/transporter_time_position/transporter_time_position_970.txt',[long,float,float])
# h1x = [[TH1F('h1x'+str(k)+'_'+str(j),'',nbin1,-r1,r1) for j in range(8)] for k in range(20)]
for ev in progress(t1,n=t1.GetEntries(),show=3):
  if ev.n_cl!=1: continue
  if ev.n_tag!=1: continue
  if ev.trigger not in [1,2]: continue
  if ev.status[0]%2==0: continue
  if abs(ev.tg[0]-ev.tcl[0])>40: continue
  if ev.nh[0]<4: continue
  if ev.Eg[0]>=1000 or ev.Eg[0]<200: continue
  if abs(1-ev.E[0]/ev.Eg[0])>0.1: continue
  if abs(ev.x[0])>330: continue
  # if ev.ypos<-360: [xpos,ypos] = correct_w_time(lcor970,ev.time)
  # elif ev.ypos<=-338.614: [xpos,ypos] = correct_w_time(lcor969,ev.time)
  # else: [xpos,ypos] = correct_w_time(lcor968,ev.time)  
  [xpos,ypos] = correct_w_time(lcor,ev.time) 
  if abs(ypos+276.4)>0.1: continue
  _ = h1x[int((ev.Eg[0]-200)/100)].Fill(ev.x[0]-xpos)
  _ = h1y[int((ev.Eg[0]-200)/100)].Fill(ev.y[0]-ypos)
  # _ = h1x[int(-(ev.y[0]+310)/5)][int((ev.Eg[0]-200)/100)].Fill(ev.x[0]-xpos)

sigma_x, sigma_y, err_x, err_y = [],[],[],[]
for k in range(20):
  sigma_x2, err_x2, = [],[]
  for j in range(8):
    def fpyjx(x,par): return fpy(x,par,j,'x') 
    ffitx = TF1('ffitx',fpyjx,-30,30,3)
    ffitx.SetParLimits(2,0.1,20)
    ffitx.SetParameter(2,3.)
    h1x[k][j].Fit(ffitx,'L')
    sigma_x2.append(ffitx.GetParameter(2))
    err_x2.append(ffitx.GetParError(2))
  sigma_x.append(sigma_x2)
  err_x.append(err_x2)
  # def fpyjy(x,par): return fpy(x,par,j,'y') 
  # ffity = TF1('ffity',fpyjy,-30,30,3)
  # ffity.SetParLimits(2,0.1,20)
  # ffity.SetParameter(2,3.)
  # h1y[j].Fit(ffity,'L')
  # sigma_y.append(ffity.GetParameter(2))
  # err_y.append(ffity.GetParError(2))

sigma_r = [(0.5*(x**2+y**2))**0.5 for x,y in zip(sigma_x,sigma_y)]
err_r = [(0.5*(x**2+y**2))**0.5 for x,y in zip(err_x,err_y)]

gx = tgraph([0.25+0.1*j for j in range(8)],sigma_x,dy=err_x)
gy = tgraph([0.25+0.1*j for j in range(8)],sigma_y,dy=err_y)
gr = tgraph([0.25+0.1*j for j in range(8)],sigma_r,dy=err_r)

ffit = TF1('ffit','[0]+[1]/sqrt(x)+[2]/x')

gx = [tgraph([0.25+0.1*j for j in range(8)],sigma_x[k],dy=err_x[k]) for k in range(20)]
lx = []
for k in range(20):
  gx[k].Fit(ffit)
  lx.append(ffit.Eval(1.))

gtrans = tgraph([-312.5-k*5. for k in range(20)],lx)

gx.Fit(ffit)
value = gx.Eval(1.)
lx.append(value)
c,t = custom(gx,first=True,color=ROOT.kBlue,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{x} (mm)',margin=[0.14,0.12],insert=[[0.5,0.7,0.9,0.9,'#sigma_{x}(1 GeV) = '+str(round(value,2))+' mm',ROOT.kRed]])
c.Print('img/161202/sigmax_'+str(run)+'_primex.pdf')
gy.Fit(ffit)
value = gy.Eval(1.)
c,t = custom(gy,first=True,color=ROOT.kBlue,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{y} (mm)',margin=[0.14,0.12],insert=[[0.5,0.7,0.9,0.9,'#sigma_{y}(1 GeV) = '+str(round(value,2))+' mm',ROOT.kRed]])
c.Print('img/161202/sigmay_'+str(run)+'_primex.pdf')
gr.Fit(ffit)
value = gr.Eval(1.)
c,t = custom(gr,first=True,color=ROOT.kBlue,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{r} (mm)',margin=[0.14,0.12],insert=[[0.5,0.7,0.9,0.9,'#sigma_{r}(1 GeV) = '+str(round(value,2))+' mm',ROOT.kRed]])
c.Print('img/161202/sigmar_'+str(run)+'_primex.pdf')


# xpos and ypos correction
runs =  sum_dict_lists(runs_list,'snake')
for x in runs:
  f = TFile('/work/hallb/prad/mlevilla/calib_results/calib/pass13/tree_'+str(x)+'.root')
  t = f.hycal
  l = []
  for i,e in progress(enumerate(t),show=3,n=t.GetEntries()):
    if i==0 or t.xpos!=l[-1][1] or t.ypos!=l[-1][2]:
      l.append([t.time,t.xpos,t.ypos])
  print x,len(l)
  writelist('transporter_time_position_'+str(x)+'.txt',l)


# efficiency plots
methods = [['calib','13'],['fsnake','0']]
folders = ['/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/'.format(m[0],m[1]) for m in methods]
cuts = ['0','t40','t40_xy','t40_xy_s4','t40_xy_s3','t40_xy_s2','t40_xy_s3_n2','t40_xy_s3_n3']
cuts = ['t40_xy_s3']
f = [[TFile(fo+'/d_'+c+'/eff.root') for c in cuts] for fo in folders]
regions = ['lg_center','pwo_center']
regions = ['edge','lg_center','trans_pwo','pwo_center']
g = [[[y.Get(re).Get('g0_'+re) for re in regions] for y in x] for x in f]
2 #sigma_{E/E_{#gamma}}

# photon beam profile 

f0 = TFile('/work/hallb/prad/mlevilla/photon_beam_profile.root')
t0 = f0.tree
nbin0, r0 = 1000,100
h0 = TH1F('h0','',nbin0,-r0,r0)

for ev in progress(t0,n=t0.GetEntries()):
  xtemp = ev.x+7500*ev.tx/(1-ev.tx**2-ev.ty**2)**0.5
  ytemp = ev.y+7500*ev.ty/(1-ev.tx**2-ev.ty**2)**0.5
  if (xtemp**2+ytemp**2)>50**2: continue
  _ = h0.Fill(ev.x+13500*ev.tx/(1-ev.tx**2-ev.ty**2)**0.5)

def fconv(x,par): return par[0]*sum([h0[int(nbin0/2./r0*(i+r0))+1]*math.exp(-0.5*((x[0]-i-par[1])/par[2])**2) for i in range(-r0,r0)])

ffitx = TF1('ffitx',fconv,-30,30,3)
ffitx.SetParameters(1.,0.,3.)
ffitx.SetParLimits(1,-10,10)
ffitx.SetParLimits(2,1.5,15.)
h1x[3].Fit(ffitx,'RL','',-20,20)

ftest = TF1('ftest','[0]*exp(-0.5*(x-[1])^2/[2]^2)/(1+(x-[1])^2/[3]^2)+[4]',-30,30)
h0.Fit(ftest,'L')


ftest2 = TF1('ftest2','[0]*exp(-0.5*(x-[1])^2/[2]^2)*(1+erf([3]*(x-[1])))',-30,30)

fconv = TF1Convolution('exp(-0.5*x^2/31^2)/(1+x^2/13.6^2)*(1+erf(-0.133))','[0]*exp(-0.5*(x-[1])^2/[2]^2)',-10,10)
fconv2 = TF1('fconv2',fconv,-30,30,3)


# cross check physics calibration
namesw = readfile('/work/hallb/prad/mlevilla/finalconstant/period5_subperiod2/cal_period52.dat',[str],cols=[0])
alphaw = readfile('/work/hallb/prad/mlevilla/finalconstant/period5_subperiod2/cal_period52.dat',[float],cols=[3])
gainsw = readfile('/work/hallb/prad/mlevilla/finalconstant/period5_subperiod2/cal_period52.dat',[float],cols=[1])
gainsw = [gainsw[namesw.index(x)] for x in module_names_prim]
alphaw = [alphaw[namesw.index(x)] for x in module_names_prim]
ratio = readfile('/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/d_{2}/ep_{3}_{4}_{2}.txt'.format('fprod',0,1,'1238_1341','all'),[float],cols=[1])
alpha = readfile('/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/d_{2}/alpha_{3}_{4}_{2}.txt'.format('fprod',0,1,'1238_1341','all'),[float],cols=[1])

from viewer import *
valpha = hcviewer('alpha',zmin=0.2,zmax=2.)
valpha.fill([[x/y if y!=0. else 0. for x,y in zip(alpha,alphaw)]],[module_names_prim])
valpha.draw()



# normalized efficicency map
method,ipass = 'fsnake',0
folder = '/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/d_t40_xy_s3/'.format(method,ipass)
eff = readfile(folder+'eff_modules.txt',[float]*17,cols=range(1,18))
stat1 = readfile(folder+'stat_tr1.txt',[float]*17,cols=range(1,18))
stat2 = readfile(folder+'stat_tr2.txt',[float]*17,cols=range(1,18))
stat = [[x1+y1 for x1,y1 in zip(x2,y2)] for x2,y2 in zip(stat1,stat2)]
mean_eff = [sum([x[i] for i in range(5,15) if x[i]>0])/len([x[i] for i in range(5,15) if x[i]>0]) if len([x[i] for i in range(5,15) if x[i]>0])!=0 else 0. for x in eff]
sum_stat = [sum([x[i] for i in range(5,15)]) for x in stat]
sum_stat0 = prad_to_primex(sum_stat,-1)
mean_eff0 = prad_to_primex(mean_eff,0)
zones = [range(1,25)+range(31,55)+range(61,85)+range(91,115),[x for x in range(1001) if ((x not in (range(1,25)+range(31,55)+range(61,85)+range(91,115))) and (module_names_primex[x] in module_names))],range(1001,1478),range(1478,1749),range(1749,1919),range(1919,2157)]

h = [TH1F('h'+str(i),'',1000,0.9,1.1) for i in range(6)]
for i in range(2157):
  for j in range(6):
    if i in zones[j]: _ = h[j].Fill(mean_eff0[i])

zone_means = [0. for i in range(6)]
zone_sigmas = [0. for i in range(6)]
for i in range(6):
  h[i].Draw()
  s = raw_input()
  r = h[i].Fit('gaus','rS')
  zone_means[i] = r.GetParams()[1]
  zone_sigmas[i] = r.GetParams()[2]

choice_mean = 0.995
choice_sigma = 1e-3
mean_eff1 = []
for i in range(2157):
  belong = False
  for j in range(6):
    if i in zones[j]:
      belong = True
      #mean_eff1.append(mean_eff0[i]*choice/zone_means[j])
      mean_eff1.append(choice_mean+choice_sigma*(mean_eff0[i]-zone_means[j])/zone_sigmas[j])
      break
  if not belong: mean_eff1.append(0)

correct_av_lg = [10,101,102,130,131,132,133,134,135,136,159,393,394,783,784,793,794,795,813,820,832,900]
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
    if module_names_primex[j] not in module_names: continue
    l.append(j)
  if i in not_working_lg+not_working_pwo: mean_eff1[i] = 0.
  else: mean_eff1[i] = sum([mean_eff1[j] for j in l])/len(l)

g = tgraph(range(2157),mean_eff1)
g.Draw('ap')
v = hcviewer('eff',False,0.95,1.0)
v.fill([primex_to_prad(mean_eff1)])
v.draw()


#efficciency by zones
method,ipass = 'fsnake',0
folder = '/work/hallb/prad/mlevilla/calib_results/{0}/pass{1}/d_t40_xy_s3/'.format(method,ipass)
eff = readfile(folder+'eff_modules.txt',[float]*17,cols=range(1,18))
eff = prad_to_primex(eff,-1)
stat1 = readfile(folder+'stat_tr1.txt',[float]*17,cols=range(1,18))
stat2 = readfile(folder+'stat_tr2.txt',[float]*17,cols=range(1,18))
stat5 = readfile(folder+'stat_tr5.txt',[float]*17,cols=range(1,18))
stat = [[x1+y1+z1 for x1,y1,z1 in zip(x2,y2,z2)] for x2,y2,z2 in zip(stat1,stat2,stat5)]
zones = [range(1,25)+range(31,55)+range(61,85)+range(91,115),[x for x in range(1001) if ((x not in (range(1,25)+range(31,55)+range(61,85)+range(91,115))) and (module_names_primex[x] in module_names))],range(1001,1478),range(1478,1749),range(1749,1919),range(1919,2157)]
exclude = [10,101,102,130,131,132,133,134,135,136,159,393,394,783,784,793,794,795,813,820,832,900,1393,1419,1420,1470,1471,1472,1527,1528,1560,1563,1629,1630,1637,1752,1824,1835,1891]
heff = [[TH1F('eff'+str(i)+'_'+str(j),'',1000,0.91,1.01) for j in range(17)] for i in range(len(zones))]
hstat = [TH1F('stat'+str(i),'',17)]
for i in range(2157):
  if i in exclude: continue
  for j in range(len(zones)):
    if i in zones[j]:
      for k in range(17): _ = h[j][k].Fill(eff[i][k])
    break
              

# moller yield
from moller import *
from energy_loss import *
from xsect import *
from geo3D import *
run = 1342
[names,ecalib,alpha] = readfile('/work/hallb/prad/mlevilla/finalconstant/period5_subperiod5/cal_period55.dat',[str,float,float],cols=[0,2,3],tr=True)
thbin = array('d',[0.5+0.05*i for i in range(111)])
opthbin = array('d',[2+0.05*i for i in range(111)])
q2bin = array('d',[2*1.1*moller_energy(th*degrad,1.1)*(1-cos(th*degrad)) for th in thbin])
f = TFile('/work/hallb/prad/mlevilla/calib_results/fprod_gem/pass1/trees/tree_primex_gem_'+str(run)+'.root')
t = f.hycal
fout = TFile('histo_yields3_'+str(run)+'.root','recreate')
hth = [[TH1F('hth'+str(i)+'_'+str(j),';#theta (deg);',110,thbin) for i in range(3)] for j in range(3)]
hth2 = [[TH1F('hth2'+str(i)+'_'+str(j),';#theta (deg);',110,thbin) for i in range(3)] for j in range(3)]
hopth1 = [[TH1F('hopth1'+str(i)+'_'+str(j),';#theta (deg);',110,opthbin) for i in range(3)] for j in range(3)]
hopth2 = [[TH1F('hopth2'+str(i)+'_'+str(j),';#theta (deg);',110,opthbin) for i in range(3)] for j in range(3)]
hq2 = [[TH1F('hq2'+str(i)+'_'+str(j),';Q^{2} (GeV^{2});',110,q2bin) for i in range(3)] for j in range(3)]
hphidiff = [[TH1F('hphidiff_'+str(i)+'_'+str(j),';#Delta#phi (deg);',1000,-180,180) for i in range(3)] for j in range(2)]
helas = [[TH1F('helas'+str(i)+'_'+str(j),';(#sum E_{rec})/E_{beam}-1;',1000,-1,1) for i in range(3)]  for j in range(2)]
hratio = [[TH1F('hratio_'+str(i)+'_'+str(j),';E_{rec}/E_{theo};',1000,0,3) for i in range(3)] for j in range(3)]
hz = [[TH1F('hz'+str(i)+'_'+str(j),';z_{vertex} (m);',1000,-30,30) for i in range(3)] for j in range(3)]
hr = [[TH1F('hr'+str(i)+'_'+str(j),';r_{vertex} (mm);',1000,-300,300) for i in range(3)] for j in range(3)]
hd = [[TH1F('hd'+str(i)+'_'+str(j),';d_{vertex} (mm);',1000,-500,500) for i in range(3)] for j in range(3)]
        
for _ in progress(t,show=2,n=t.GetEntries(),modulo=10000):
  #preselection
  indexes = []
  for j in range(t.n_cl):
    # if t.nh[j]<2 or t.E[j]<50: continue
    if any([t.zgem[j][k]<5000 for k in range(2)]): continue
    if t.nh[j]>15: continue
    indexes.append(j)
  # loop on cluster pairs
  if len(indexes)<2: continue
  for i1 in range(len(indexes)):
    for i2 in range(i1+1,len(indexes)):
      # variables
      i = [indexes[i1],indexes[i2]]
      cgem = [[[t.xgem[j][k],t.ygem[j][k],t.zgem[j][k]] for k in range(2) if t.zgem[j][k]>5000] for j in i]
      for j in range(2):
        if len(cgem[j])==1: cgem[j]==cgem[j][0]
        else:
          cgem[j] = [proj(cgem[j][k],[0,0,1,-0.5*(cgem[j][0][2]+cgem[j][1][2])]) for k in range(2)]
          cgem[j] = [0.5*(cgem[j][0][k]+cgem[j][1][k]) for k in range(3)]
      if cgem[0]==cgem[1]: continue
      rgem = [(c[0]**2+c[1]**2)**0.5 for c in cgem]
      chycal = [[t.xhycal[j],t.yhycal[j],t.zhycal[j]] for j in i]
      rhycal = [(c[0]**2+c[1]**2)**0.5 for c in chycal]
      zhycal = ((m_e+t.Ebeam/1000.)*rhycal[0]*rhycal[1]/2./m_e)**0.5
      ids = [t.id[j] for j in i]
      mnames = ['W'+str(j-1000) if j>1000 else 'G'+str(j) for j in ids]
      idms = [names.index(na) for na in mnames]
      E = [t.E[j]/(1+alpha[k]*(t.E[j]-ecalib[k])/1000.) for j,k in zip(i,idms)]
      # quantities
      elas = sum(E)/t.Ebeam-1.
      vertex = clst_l2l(chycal[0],cgem[0],chycal[1],cgem[1])
      cvertex = mult(2,plus(vertex[0],vertex[1]))
      dvertex = dist(diff(vertex[0],vertex[1]))
      if cgem[0][2]!=cgem[1][2]: cgem1p = proj(cgem[1],[0,0,1,-cgem[0][2]])
      else: cgem1p = cgem[1]
      phidiff = angle(cgem[0],cgem1p,[0,0,cgem[0][2]])%(2*pi)-pi)/degrad
      opth1 = angle(cgem[0],cgem[1])/degrad
      opth2 = angle(cgem[0],cgem[1],cvertex)/degrad
      if opth1==None or opth2==None: print 'opth error'; continue
      for j in range(2):
        theta = angle(cgem[j])
        theta2 = angle(cgem[j],vertex[j],vertex[j][:2]+[cgem[j][2]])
        q2 = 2*t.Ebeam*E[j]*1e-6-2*m_e**2-2*((t.Ebeam*1e-3)**2-m_e**2)**0.5*((E[j]*1e-3)**2-m_e**2)**0.5*cos(theta)
        if theta==None or theta2==None: break
        exp_em_energy = 1000*moller_energy(theta,t.Ebeam/1000.)
        exp_em_energy -= energy_loss(theta,exp_em_energy)
        ratio = E[j]/exp_em_energy
        k = int(ids[j]>1000)
        #  histo filling 
        for k2 in [k,2]:
          _ = hratio[0][k2].Fill(ratio)
          _ = hphidiff[0][k2].Fill(phidiff)
          _ = helas[0][k2].Fill(elas)
          _ = hz[0][k2].Fill(vertex[j][2]/1000.)
          _ = hr[0][k2].Fill(dist(vertex[j][:2]))
          _ = hd[0][k2].Fill(dvertex)
          _ = hopth1[0][k2].Fill(opth1)
          _ = hopth2[0][k2].Fill(opth2)
          _ = hth[0][k2].Fill(theta/degrad)
          _ = hth2[0][k2].Fill(theta2/degrad)
          _ = hq2[0][k2].Fill(q2)
        if abs(phidiff)<3: 
          _ = helas[1][k].Fill(elas)
          _ = helas[1][2].Fill(elas)
        if abs(elas)<0.1: 
          _ = hphidiff[1][k].Fill(phidiff)
          _ = hphidiff[1][2].Fill(phidiff)
        if abs(elas)>0.1: break
        for k2 in [k,2]:
          _ = hratio[1][k2].Fill(ratio)
          _ = hz[1][k2].Fill(vertex[j][2]/1000.)
          _ = hr[1][k2].Fill(dist(vertex[j][:2]))
          _ = hd[1][k2].Fill(dvertex)
          _ = hopth1[1][k2].Fill(opth1)
          _ = hopth2[1][k2].Fill(opth2)
          _ = hth[1][k2].Fill(theta/degrad)
          _ = hth2[1][k2].Fill(theta2/degrad)
          _ = hq2[1][k2].Fill(q2)
        if abs(phidiff)>3: break
        for k2 in [k,2]: 
          _ = hratio[2][k2].Fill(ratio)
          _ = hz[2][k2].Fill(vertex[j][2]/1000.)
          _ = hr[2][k2].Fill(dist(vertex[j][:2]))
          _ = hd[2][k2].Fill(dvertex)
          _ = hopth1[2][k2].Fill(opth1)
          _ = hopth2[2][k2].Fill(opth2)
          _ = hth[2][k2].Fill(theta/degrad)
          _ = hth2[2][k2].Fill(theta2/degrad)
          _ = hq2[2][k2].Fill(q2)

fout.cd()
for x in [hphidiff,helas,hratio,hz,hr,hd,hth,hq2,hopth1,hopth2]:
  for y in x:
    for z in y:
      z.Write()

fout.Close()

runs = [1342,1345]
bcharge = [34789.5,3996.96]
f = [TFile('histo_yields2_'+str(run)+'.root') for run in runs]


h,n = hr,3
for i in range(2):
  for j in range(n):
    if j==0: h[j][i].Draw()
    else: h[j][i].Draw('same')
  s = raw_input()



# epics tree

l = [int(x[10:14]) if x[10]=='1' else int(x[10:13]) for x in os.listdir('/work/hallb/prad/mlevilla/stability') if 'gains_LMS_' in x and '.root' in x and 'old' not in x]
l.sort()
fout = TFile('/work/hallb/prad/mlevilla/stability/epics_tree.root','recreate')
f0 = TFile('/work/hallb/prad/mlevilla/stability/gains_LMS_1288.root')
names = [x.GetName() for x in f0.epics.GetListOfKeys()]

emptyp = get_runs_between(1059,1516,'empty')
prodp = get_runs_between(1059,1516,'prod')
emptycp = get_runs_between(1059,1516,'total_empty')
carbonp = get_runs_between(1059,1516,'carbon')

lv = [['run','F','',1],['type','I','',1],['nev','I','',1],['iev','I','nev',2500]]+[[x[1:].replace(':','_'),'F','nev',2500] for x in names]
#lv = [['run','F','',1],['iev','I','',1]]+[[x,'F','',1] for x in names]
fout.cd()
t,v = newtree('tree',lv,'tree')
for i,x in enumerate(l):
  print i,x
  f = TFile('/work/hallb/prad/mlevilla/stability/gains_LMS_'+str(x)+'.root')
  v['run'][0] = x
  if x<1059: v['type'][0]=0
  elif 1059<=x<=1345:
    if x in emptycp: v['type'][0]=1
    elif x in emptyp: v['type'][0]=2
    elif x in carbonp: v['type'][0]=3
    elif x in prodp: v['type'][0]=4
  elif 1362<=x<=1516:
    if x in emptycp: v['type'][0]=5
    elif x in emptyp: v['type'][0]=6
    elif x in carbonp: v['type'][0]=7
    elif x in prodp: v['type'][0]=8
  n = f.epics.Get(names[0]).GetN()
  v['nev'][0] = n
  xx = f.epics.Get(names[0]).GetX()
  yy = [f.epics.Get(y).GetY() for y in names]
  for i in range(n):
    v['iev'][i] = int(xx[i])
    for j,y in enumerate(names):
      v[y[1:].replace(':','_')][i] = yy[j][i]
  _ = t.Fill()
  f.Close()


fout.Close()

# zhycal and zgem from kinematics

run = 1345
f = TFile('/work/hallb/prad/mlevilla/calib_results/fprod_gem/pass1/trees/tree_primex_gem_'+str(run)+'.root')
t = f.hycal
hzhycal = [TH1F('hzhycal'+str(i),';z_{vertex} (mm);',1000,4000,6500) for i in range(3)]
hzgem1 = TH1F('hzgem1',';z_{vertex} (mm);',1000,4000,6500)
hzgem2 = TH1F('hzgem2',';z_{vertex} (mm);',1000,4000,6500)

for _ in progress(t,show=2,n=t.GetEntries(),modulo=10000):
  #preselection
  indexes = []
  for i in range(t.n_cl):
    if t.nh[i]<2 or t.nh[i]>15: continue
    if t.E[i]<50: continue
    indexes.append(i)
  for i1 in range(len(indexes)-1):
    for i2 in range(i1+1,len(indexes)):
      # variables
      i = [indexes[i1],indexes[i2]]
      elas = sum([t.E[j] for j in i])/t.Ebeam-1
      chycal = [[t.xhycal[j],t.yhycal[j],t.zhycal[j]] for j in i]
      rhycal = [(c[0]**2+c[1]**2)**0.5 for c in chycal]
      zhycal = ((m_e+t.Ebeam/1000.)*rhycal[0]*rhycal[1]/2./m_e)**0.5
      phidiff_hycal = angle(chycal[0][:2]+[0.],chycal[1][:2]+[0.])%(2*pi)-pi
      if abs(elas)<0.1 and abs(phidiff_hycal/degrad)<5: 
        if t.id[indexes[i1]]<1000 and t.id[indexes[i2]]<1000: _ = hzhycal[0].Fill(zhycal)
        else: _ = hzhycal[1].Fill(zhycal)
      cgem = [[[t.xgem[2*j+k],t.ygem[2*j+k],t.zgem[2*j+k]] for k in range(2)] for j in i]
      if cgem[0][0][2]>0 and cgem[1][0][2]>0:
        rgem1 = [(c[0][0]**2+c[0][1]**2)**0.5 for c in cgem]
        zgem1 = ((m_e+t.Ebeam/1000.)*rgem1[0]*rgem1[1]/2./m_e)**0.5
        phidiff_gem1 = angle(cgem[0][0][:2]+[0.],cgem[1][0][:2]+[0.])%(2*pi)-pi
        if abs(elas)<0.1 and abs(phidiff_gem1/degrad)<5: _ = hzgem1.Fill(zgem1)
      if cgem[0][1][2]>0 and cgem[1][1][2]>0:
        rgem2 = [(c[1][0]**2+c[1][1]**2)**0.5 for c in cgem]
        zgem2 = ((m_e+t.Ebeam/1000.)*rgem2[0]*rgem2[1]/2./m_e)**0.5
        phidiff_gem2 = angle(cgem[0][1][:2]+[0.],cgem[1][1][:2]+[0.])%(2*pi)-pi
        if abs(elas)<0.1 and abs(phidiff_gem2/degrad)<5: _ = hzgem2.Fill(zgem2)


# pressure and temperature
runs = get_runs_between(1059,1516)
lf = []
for run in runs:
  l = readfile('/work/hallb/prad/replay_EnS/EventSelect_{0}.txt'.format(str(run).zfill(6)),[str,str])
  if l==[]: l = [[]]
  if type(l[0])!=list: l = [l]
  if l!=[[]] and l[0][0][0]=='*': l[0][0] = l[0][0][1:]
  l = [[int(x[0]),int(x[1])] for x in l if x!=[]] 
  f = TFile('/work/hallb/prad/mlevilla/stability/gains_LMS_{0}.root'.format(run))
  gT = f.epics.Get('gTGT:PRad:Cell_Gas_T')
  gP = f.epics.Get('gTGT:PRad:Cell_P')
  n = gT.GetN()
  x = [gT.GetX()[i] for i in range(n)]
  x2 = []
  for i in range(n):
    good = True
    for j in range(len(l)):
      if l[j][0]<x[i]<l[j][1]: 
        good = False
        break
    if good: x2.append(x[i])
  yT = [gT.GetY()[i] for i in range(n) if x[i] in x2]
  yP = [gP.GetY()[i] for i in range(n) if x[i] in x2]
  gT2 = tgraph(x2,yT)
  gP2 = tgraph(x2,yP)
  Tmean = round(gT2.GetMean(2),2)
  Tsigma = round(gT2.GetRMS(2),2)
  Pmean = round(gP2.GetMean(2),3)
  Psigma = round(gP2.GetRMS(2),3)
  thickness = '%.3e'%(1.93e16*4.*Pmean/Tmean) if Pmean>0 else 0
  dthickness = '%.3e'%(float(thickness)*((Tsigma/Tmean)**2+(Psigma/Pmean)**2)**0.5)
  lf.append([run,Tmean,Tsigma,Pmean,Psigma,thickness,dthickness])

writelist('target_info.txt',lf,ljust=10)


# moller acceptance

from xsect import *
from moller import *
load_names()
load_positions()
hall = [TH1F('hall'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hcut = [TH1F('hcut'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hcut2 = [TH1F('hcut2'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hdcut = [TH1F('hdcut'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hdcut2 = [TH1F('hdcut2'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hacc = [TH1F('hacc'+str(i),';#theta (deg);',1000,0.,10.) for i in range(2)]
hycal_eff = readfile('map_eff_normalized.txt',[str,float],out=dict)
hycal_eff = {x:hycal_eff[y] for x,y in module_names.items()}
gem_eff0 = readfile('sector_efficiency.txt',[float,float],cols=[2,3])
gem_eff = {x:y[0] for x,y in zip(module_names.keys(),gem_eff0)}
gem_deff = {x:y[1] for x,y in zip(module_names.keys(),gem_eff0)}
gem_acc = readfile('geo_acceptance.txt',[float,float,float],cols=[0,1,4])
edge = get_regions('inner')+get_regions('outer')
n=1000000

def gem_accept(theta):
  for x in gem_acc:
    if x[0]<=theta<x[1]: return x[2]
  return 0.

for i in progress(n,precision=2):
  theta = random.uniform(0.5*degrad,8.5*degrad)
  ctheta = complementary_angle(theta,1.1)
  phi = random.uniform(0.,2*pi)
  x = tan(theta)*cos(phi)*5720.
  y = tan(theta)*sin(phi)*5720.
  xc = tan(ctheta)*cos(phi+pi)*5720.
  yc = tan(ctheta)*sin(phi+pi)*5720.
  theta/=degrad
  ctheta/=degrad
  _=hall[0].Fill(theta)
  _=hall[0].Fill(ctheta)
  _=hall[1].Fill(theta)
  m = which_module(x,y)
  mc = which_module(xc,yc)
  gacc = gem_accept(theta)
  cgacc = gem_accept(ctheta)
  # gacc = 1.
  # cgacc = 1.
  if m==mc==-1 or (m==-1 and mc in edge) or (mc==-1 and m in edge) or (m in edge and mc in edge): continue
  if mc==-1 or (mc in edge):
    _=hcut[1].Fill(theta,hycal_eff[m]*gem_eff[m]*gacc)
    if gem_deff[m]<0.1:
      _=hcut2[1].Fill(theta,hycal_eff[m]*gem_eff[m]*gacc)
      _=hdcut[1].Fill(theta,(1.5e-3)**2+gem_deff[m]**2)
  if mc!=-1 and m!=-1 and (m not in edge) and (mc not in edge):
    _=hcut[0].Fill(theta,hycal_eff[m]*gem_eff[m]*hycal_eff[mc]*gem_eff[mc]*gacc*cgacc)
    _=hcut[0].Fill(ctheta,hycal_eff[m]*gem_eff[m]*hycal_eff[mc]*gem_eff[mc]*gacc*cgacc)
    _=hcut[1].Fill(theta,hycal_eff[m]*gem_eff[m]*gacc)
    if gem_deff[m]<0.1:
      _=hcut2[1].Fill(theta,hycal_eff[m]*gem_eff[m]*gacc)
      _=hdcut[1].Fill(theta,(1.5e-3)**2+gem_deff[m]**2)
    if gem_deff[m]<0.1 and gem_deff[mc]<0.1:
      _=hcut2[0].Fill(theta,hycal_eff[m]*gem_eff[m]*hycal_eff[mc]*gem_eff[mc]*gacc*cgacc)
      _=hdcut[0].Fill(theta,(1.5e-3)**2+gem_deff[m]**2+gem_deff[mc]**2)
      _=hcut2[0].Fill(ctheta,hycal_eff[m]*gem_eff[m]*hycal_eff[mc]*gem_eff[mc]*gacc*cgacc)
      _=hdcut[0].Fill(ctheta,(1.5e-3)**2+gem_deff[m]**2+gem_deff[mc]**2)
for i in range(2):
  hdcut2[i].Divide(hdcut[i],hcut2[i])
  for j in range(1000): hdcut2[i][j+1]=hdcut2[i][j+1]**0.5
  hacc[i].Divide(hcut[i],hall[i])
  for j in range(1000):
    if hall[i][j+1]!=0:
      hacc[i].SetBinError(j+1,(hacc[i][j+1]*(1-hacc[i][j+1])/hall[i][j+1])**0.5)

l = []
for i in range(1000):
  l.append([round(hacc[0].GetBinCenter(i+1),3),round(hacc[0][i+1],4),round(hdcut[0][i+1],4),round(hacc[1][i+1],4),round(hdcut[1][i+1],4),round(hacc[0].GetBinLowEdge(i+1),3),round(hacc[0].GetBinLowEdge(i+2),3)])

writelist('acceptance_theta.txt',l,ljust=8)


# comparison cross section
from xsect import *
f = TFile(work+'/moller_xs/moller_xs_1238_1341.root')
h = f.hth_all_cross_eff
thetas = [h.GetBinCenter(i+1) for i in range(200)]
ev = [MollerEvent(theta1=x*degrad,Ei=1.1) for x in thetas]
dsigma_domega = [e.dsigma_domega()*389.379e-6 for e in ev]
g = tgraph(thetas,dsigma_domega)
ratio = [x/y if y!=0 else -1 for x,y in zip([h[i+1] for i in range(200)],dsigma_domega)]
gratio = tgraph(thetas,ratio)


# plot cross_section
from moller import *
gStyle.SetOptStat(0)
plotdir = 'img/170317/'
l = ['sigma4_phi5']
l2 = [[1238,1341]]
f = [[TFile(work+'/moller_xs/{0}/moller_xs_{1}_{2}.root'.format(x,y[0],y[1])) for y in l2] for x in l]
# acceptance
fy = f[0][0]
for x in ['1arm','2arm']:
  c = TCanvas()
  h = fy.acceptance.Get('hacc_'+x)
  (c,t) = custom(h,c=c,xrg=[0.,7.],yrg=[0.3,1.1],margin=[0.12,0.12],offset=[0.8,0.9])
  c.Print(plotdir+'acceptance_'+x+'.pdf')

# selection
fy = f[0][0]
c = TCanvas()
c.SetLogz()
for i in range(2):
  for j in range([2,4][i]):
    h = fy.selection.Get('h2cut_'+str(i+1)+'arm_'+str(j))
    (c,t) = custom(h,c=c,yrg=[0.,1500.],xrg=[0.5,8.],first=True)
    c.Update()
    raw_input()
    c.Print(plotdir+'h2cut_'+str(i+1)+'arm_'+str(j)+'_1GeV.pdf')

# cross_sections
name,abb = 'theta','th'
xrg = [0.5,5.5]
(c,t) = custom(fy.Get(name).Get('h'+abb+'_prod_1arm_all'),margin=[0.13,0.12],offset=[0.8,0.9],xrg=xrg,yrg=[0.001,10],log=1,width=2)
(c,t) = custom(fy.Get(name).Get('h'+abb+'_prod_2arm_all'),first=False,color=ROOT.kRed,c=c,width=2)
(c,t) = custom(fy.Get(name).Get('h'+abb+'_empty_1arm_all'),first=False,color=ROOT.kBlue,c=c,width=2)
(c,t) = custom(fy.Get(name).Get('h'+abb+'_empty_2arm_all'),first=False,color=ROOT.kGreen+2,c=c,width=2,legend=[0.6,0.7,0.9,0.9])
c.Print(plotdir+'yields_'+name+'_1GeV.pdf')
(c,t) = custom(fy.Get(name).Get('h'+abb+'_ratio_1arm_all'),margin=[0.13,0.12],offset=[0.8,0.9],xrg=xrg,yrg=[0.001,0.1],width=2)
(c,t) = custom(fy.Get(name).Get('h'+abb+'_ratio_2arm_all'),first=False,color=ROOT.kRed,c=c,width=2,legend=[0.6,0.7,0.9,0.9])
c.Print(plotdir+'ratiobg_'+name+'_1GeV.pdf')
(c,t) = custom(fy.Get(name).Get('g'+abb+'exp_1arm'),margin=[0.13,0.12],offset=[0.8,0.7],xrg=xrg,yrg=[0.05,5.],width=2,log=1,style=[0,0,1],title=';#theta (deg);#frac{d#sigma}{d#Omega} (b)')
(c,t) = custom(fy.Get(name).Get('g'+abb+'exp_2arm'),first=False,color=ROOT.kRed,c=c,width=2,style=[0,0,1])
(c,t) = custom(fy.Get(name).Get('g'+abb+'theo_2arm'),first=False,color=ROOT.kBlue,c=c,width=2,legend=[0.6,0.7,0.9,0.9],style=[0,0,7])
c.Print(plotdir+'cross_'+name+'_1GeV.pdf')
(c,t) = custom(fy.Get(name).Get('g'+abb+'exp_1arm_eff'),margin=[0.13,0.12],offset=[0.8,0.7],xrg=[0.5,6.],yrg=[0.05,5.],width=2,log=1,style=[0,0,1],title=';#theta (deg);#frac{d#sigma}{d#Omega} (b)')
(c,t) = custom(fy.Get(name).Get('g'+abb+'exp_2arm_eff'),first=False,color=ROOT.kRed,c=c,width=2,style=[0,0,1])
(c,t) = custom(fy.Get(name).Get('g'+abb+'theo_2arm'),first=False,color=ROOT.kBlue,c=c,width=2,legend=[0.6,0.7,0.9,0.9],style=[0,0,7])
c.Print(plotdir+'cross_eff'+name+'_1GeV.pdf')
(c,t) = custom(fy.Get(name).Get('gratio_'+abb+'_1arm'),margin=[0.13,0.12],offset=[0.8,0.8],xrg=[0.5,6.],yrg=[0.5,1.1],width=2,style=[0,0,1],title=';#theta (deg);#frac{d#sigma_{exp}}{d#sigma_{theo}}')
(c,t) = custom(fy.Get(name).Get('gratio_'+abb+'_2arm'),first=False,color=ROOT.kRed,c=c,width=2,style=[0,0,1],legend=[0.6,0.7,0.9,0.9])
c.Print(plotdir+'ratiotheo_'+name+'_1GeV.pdf')
(c,t) = custom(fy.Get(name).Get('gratio_'+abb+'_1arm_eff'),margin=[0.13,0.12],offset=[0.8,0.8],xrg=[0.5,6.],yrg=[0.5,1.1],width=2,style=[0,0,1],title=';#theta (deg);#frac{d#sigma_{exp}}{d#sigma_{theo}}')
(c,t) = custom(fy.Get(name).Get('gratio_'+abb+'_2arm_eff'),first=False,color=ROOT.kRed,c=c,width=2,style=[0,0,1],legend=[0.6,0.7,0.9,0.9])
c.Print(plotdir+'ratiotheo_eff_'+name+'_1GeV.pdf')

#pulls
for i in range(6):
  h = fy.pulls.Get('pth_2arm_'+str(i))
  (c,t)= custom(h,margin=[0.13,0.12],offset=[0.8,0.8],title=';pulls;')
  raw_input()
  c.Print(plotdir+'pulls_2arm_1GeV_'+str(i)+'.pdf')

# launch jobs
for x in [3,4,5]:
  for y in [[1238,1345],[1443,1516],[1417,1417]]:
    print x,y
    os.system('jsub_runs.py -min {0} -max {1} -c moller_xs.py -m island -o /work/hallb/prad/mlevilla/moller_xs/sigma{2}_phi5_2 -sigmacut {2}'.format(y[0],y[1],x))

for x in [3,4,5]:
  for y in [[1238,1341],[1443,1516]]:
    print x,y
    _ = os.system('moller_xs.py -m island -o /work/hallb/prad/mlevilla/moller_xs/sigma{2}_phi5 -c {0} {1}'.format(y[0],y[1],x))

# re title trees
method = 'island'
indir = work+'/tree/'+method+'/'
outdir = work+'/trees/'+method+'/'
l = os.listdir(indir)
l2 = os.listdir(outdir)
l3 = list(set(l)-set(l2))
l3.sort()
print len(l3)
for i,x in enumerate(l3):
  print i,x
  idx1 = x.index(method)
  idx2 = x.index('.root')
  run = x[idx1+len(method)+1:idx2]
  f = TFile(indir+x)
  fout = TFile(outdir+x,'recreate')
  t = f.event
  t2 = t.CloneTree()
  t2.SetTitle(method+'_'+run)
  fout.cd()
  t2.Write()
  fout.Close()
  f.Close()
  
# fit zvertex double gaussian

l = [work+'/test_zcarbon2_'+x+'_3.root' for x in ['1288','empty','1345']]
f = [TFile(x) for x in l]
h = [x.hz for x in f]
ffit = TF1('ffit','gaus+gaus(3)+[6]',-100,100)
for i in [0,3]: ffit.SetParLimits(i,0,10000)

for i in [1,4]: ffit.SetParLimits(i,-20,20)

ffit.SetParLimits(2,10,40)
ffit.SetParLimits(5,30,300)
s = []
for x in h:
  s.append(x.Fit(ffit,'slr','',-300,300))

live_charge = readfile('/u/home/pcrad/backup/beam_charge.dat',[int,float],cols=[0,3],start=1,out=dict)
hdiff = TH1F('hdiff','',400,-500,500)
hdiff.Add(h[0],h[1],1.,live_charge[1288]/sum(live_charge[i] for i in [1289,1294,1306,1312,1317,1318,1324,1326,1327,1329,1330,1339]))


# acceptance from monte carlo
fout = TFile('geo_acceptance.root','recreate')
f = [TFile(work+'/pradsim/output/'+x+'_elastic_11.root') for x in ['moller','ep']]
frec = [TFile(work+'/pradsim/output/'+x+'_elastic_11_rec.root') for x in ['moller','ep']]
T = [x.T for x in f] 
Trec = [x.T for x in frec] 
reso_gem = 0.008
reso_vd = 0.014
h = listofhisto('h',[['m','ep'],range(3),range(3)],200,0.,10.,fout=fout,title=';#theta (deg);acceptance')

for k,t in enumerate(T):
  for k2,x in enumerate(progress(t,precision=1,n=t.GetEntries())):
    _ = Trec[k].GetEntry(k2)
    n = getattr(t,'GEN.N')
    theta_gen = [atan2((getattr(t,'GEN.Out.X')[i]**2+getattr(t,'GEN.Out.Y')[i]**2)**0.5,getattr(t,'GEN.Out.Z')[i]) for i in range(n)]
    E_theo = [getattr(t,'GEN.Out.P')[i] for i in range(n)]
    # gem acceptance
    n_gem = getattr(t,'GEM.N')
    idx_gem = [i for i in range(n_gem) if getattr(t,'GEM.PTID')[i]==0]
    tid_gem = [getattr(t,'GEM.TID')[i] for i in idx_gem]
    b_gem = [(i+1 in tid_gem) for i in range(n)]
    theta_gem = [atan2((getattr(t,'GEM.In.X')[i]**2+getattr(t,'GEM.In.Y')[i]**2)**0.5,getattr(t,'GEM.In.Z')[i]+3000) for i in idx_gem]
    # vd acceptance
    n_vd = getattr(t,'VD.N')
    idx_vd = [i for i in range(n_vd) if getattr(t,'VD.PTID')[i]==0]
    tid_vd = [getattr(t,'VD.TID')[i] for i in idx_vd]
    x_vd = [getattr(t,'VD.In.X')[i] for i in idx_vd]
    y_vd = [getattr(t,'VD.In.Y')[i] for i in idx_vd]
    b_vd = [(i+1 in tid_vd) and (41.5<abs(x_vd[tid_vd.index(i+1)])<524.4 or 41.5<abs(y_vd[tid_vd.index(i+1)])<524.4) for i in range(n)]
    theta_vd = [atan2((getattr(t,'VD.In.X')[i]**2+getattr(t,'VD.In.Y')[i]**2)**0.5,getattr(t,'VD.In.Z')[i]+3000) for i in idx_vd]
    # histo filling
    if all(b_gem):
      for i in range(n):
        _ = h[k][0][2].Fill(degrees(theta_gen[i]))
    if all(b_vd):
      for i in range(n):
        _ = h[k][1][2].Fill(degrees(theta_gen[i]))
    if all(b_gem) and all(b_vd):
      for i in range(n):
        _ = h[k][2][2].Fill(degrees(theta_gen[i]))
    for i in range(n):
      if b_gem[i]: _ = h[k][0][1].Fill(degrees(theta_gen[i]))
      if b_vd[i]: _ = h[k][1][1].Fill(degrees(theta_gen[i]))
      if b_gem[i] and b_vd[i]: _ = h[k][2][1].Fill(degrees(theta_gen[i]))
    for i in range(n):
      _ = h[k][0][0].Fill(degrees(theta_gen[i]))
  for i in range(2): 
    for j in range(3): h[k][j][i+1].Divide(h[k][0][0])

l = [[round(0.5+(j-10)*0.05,2),round(0.5+(j-9)*0.05,2)]+[round(h[i+1][j+1],3) for i in range(2)] for j in range(10,160)]
writelist('geo_acceptance_new.txt',l,ljust=[5,5,7,7])

# pressure in the pipes

l = readfile('../vcg2.txt',cols=[0,1,3,5,7,9,11,13],start=1)
names = ['VCG2H00ATr','VCG2H01ATr','VCG2H02ATr','VCG2C21Tr','VCG2C21ATr','VCG2C24ATr']
ltime = []
lv = [[] for x in names]
for x in l:
  if x[2]=='<undefined>': continue
  t = float(x[0][8:])+float(x[1][:2])/24.+float(x[1][3:5])/24./60.
  if t>22.3: continue
  ltime.append(t)
  for i,y in enumerate(x[2:]):
    lv[i].append(float(y)*1e6)

g = [tgraph(ltime,x,name='g'+y,title=';time (days);'+y+' (#muTorr)') for x,y in zip(lv,names)]


# yields variation
sigma = [3,4,5,6]
dphi = [3,5,7,10]
zcut = [120,150,200,500]
runs = [1345,1498,1501]
bins = [1,11,18,25,51,71,171,201]

l = [[[[[0 for u in range(len(bins)-1)] for w in zcut] for x in dphi] for y in sigma] for z in [0,1]+runs] 

for i0,r in enumerate(runs):
  for i1,x in enumerate(sigma):
    for i2,y in enumerate(dphi):
      for i3,z in enumerate(zcut):
        f = TFile(work+'/yields/yields_{0}_{1}_{2}_{3}.root'.format(r,x,y,z))
        f2 = TFile(work+'/acceptance/acceptance_{0}_{1}_{2}.root'.format(x,y,z))
        h = f.yields.Get('hth_2arm_all')
        for i4 in range(len(bins)-1):
          l[i0+2][i1][i2][i3][i4] = int(h.Integral(bins[i4],bins[i4]-1))
        for i5 in range(2):
          h = f2.Get('h_{0}{0}_2arm_both_cut'.format(i5+1))
          for i4 in range(len(bins)-1):
            l[i5][i1][i2][i3][i4] = int(h.Integral(bins[i4],bins[i4]-1))
            

l2 = [[[[[l[i0+2][i1][i2][i3][i4]/l[[0,1,1][i0]][i1][i2][i3][i4] if l[[0,1,1][i0]][i1][i2][i3][i4]!=0 else 0 for i4 in range(len(bins)-1)] for i3 in range(4)] for i2 in range(4)] for i1 in range(4)] for i0 in range(3)]         

gsigma = [[[tgraph(sigma,[100.*(sum([l2[i0][i1][i2][i3][i4] for i4 in range(2,5)])/float(sum([l2[i0][0][i2][i3][i4] for i4 in range(2,5)]))-1) for i1 in range(4)],name='gsigma_{0}_{1}_{2}'.format(r,y,z),title='run='+str(r)+', |#Delta#phi|<'+str(y)+', |z_{vertex}|<'+str(z)+';#sigma_{elasticity};% yield') for i3,z in enumerate(zcut)] for i2,y in enumerate(dphi)] for i0,r in enumerate(runs)]
gphi = [[[tgraph(dphi,[100.*(sum([l2[i0][i1][i2][i3][i4] for i4 in range(2,5)])/float(sum([l2[i0][i1][0][i3][i4] for i4 in range(2,5)]))-1) for i2 in range(4)],name='gphi_{0}_{1}_{2}'.format(r,y,z),title='run='+str(r)+', #sigma_{elasticity}='+str(y)+', |z_{vertex}|<'+str(z)+';#Delta#phi cut;% yield') for i3,z in enumerate(zcut)] for i1,y in enumerate(sigma)] for i0,r in enumerate(runs)]
gz = [[[tgraph(zcut,[100.*(sum([l2[i0][i1][i2][i3][i4] for i4 in range(2,5)])/float(sum([l2[i0][i1][i2][0][i4] for i4 in range(2,5)]))-1) for i3 in range(4)],name='gz_{0}_{1}_{2}'.format(r,y,z),title='run='+str(r)+', #sigma_{elasticity}='+str(y)+', |#Delta#phi|<'+str(z)+';z_{vertex} cut;') for i2,z in enumerate(dphi)] for i1,y in enumerate(sigma)] for i0,r in enumerate(runs)]


# create tree with beamline pressure and time and run
l = readfile('../vcg2.txt',cols=[0,1,3,5,7,9,11,13],start=1)
l = [[int(x[0][8:]),int(x[1][:2])*60+int(x[1][3:5])]+[float(x[j])*1e6 for j in range(2,8)] for x in l if x[2]!='<undefined>']

l2 = readfile('../times_run.txt')
l2 = [[int(x[0]),int(x[1][8:]),int(x[2][:2])*60+int(x[2][3:5]),int(x[3][:2])*60+int(x[3][3:5])+1] for x in l2 if len(x)==4 and x[0][0]=='1' and x[2]!='' and x[3]!='']
names = ['VCG2H00AuTr','VCG2H01AuTr','VCG2H02AuTr','VCG2C21uTr','VCG2C21AuTr','VCG2C24AuTr']
lv = [['date','I','',1],['time','I','',1],['run','I','',1]]+[[x,'F','',1] for x in names]
t,v = newtree('pressure',lv,'pressure in beamline (uTorr)')
for x in l:
  v['date'][0] = x[0]
  v['time'][0] = x[1]
  v['run'][0] = 0
  for i in range(6):
    v[names[i]][0] = x[2+i]
  for y in l2:
    if x[0]!=y[1]: continue
    if x[1]<y[2] or x[1]>y[3]: continue
    v['run'][0] = y[0]
    break
  _ = t.Fill()


# simulation tree plots
elas,dphi,zcut = 4,10,200
f = [TFile('{0}_{1}_{2}_{3}.root'.format(x,elas,dphi,zcut)) for x in [work+'/acceptance/acceptance_elastic_11_moller',work+'/acceptance/acceptance_11_moller',work+'/yields/yields_1288']]
helas1 = [x.selection.Get('helas1_3') for x in f]
hdphi = [x.selection.Get('hdphi_3') for x in f]
hz = [x.selection.Get('hz_3') for x in f]

for x in [helas1,hdphi,hz]:
  for i,y in enumerate(x):
    y.Scale(1./y.Integral())
    y.SetLineColor([ROOT.kBlue,ROOT.kGreen+2,ROOT.kRed][i])

# acceptance plot elastic/real
elas,dphi,zcut = 4,10,200
f = [TFile('{0}_{1}_{2}_{3}.root'.format(x,elas,dphi,zcut)) for x in [work+'/acceptance/acceptance_'+y+'11_moller' for y in ['elastic_','']]]
h = [x.acceptance.Get('hacc_2arm') for x in f]


# fit of module profiles
from ROOT import TF2
ffit = TF2('ffit','[0]/2/pi*((atan((x+0.5)/[1])+atan((y+0.5)/[1])+atan((x+0.5)*(y+0.5)/[1]/sqrt([1]^2+(x+0.5)^2+(y+0.5)^2)))-(atan((x-0.5)/[1])+atan((y+0.5)/[1])+atan((x-0.5)*(y+0.5)/[1]/sqrt([1]^2+(x-0.5)^2+(y+0.5)^2)))-(atan((x+0.5)/[1])+atan((y-0.5)/[1])+atan((x+0.5)*(y-0.5)/[1]/sqrt([1]^2+(x+0.5)^2+(y-0.5)^2)))+(atan((x-0.5)/[1])+atan((y-0.5)/[1])+atan((x-0.5)*(y-0.5)/[1]/sqrt([1]^2+(x-0.5)^2+(y-0.5)^2))))',-5,5,-5,5)
ffit.SetParameters(1.,0.1)

f = TFile(work+'/module_profile_1288.root')
h = [[f.Get('hp_{0}_{1}'.format(x,y)) for x in range(6)] for y in [0,3]]
p = [[x.ProjectionX("_px",100,100) for x in y] for y in h] 

ffitx = TF1('ffitx','[0]/2/pi*((atan((x+0.5)/[1])+atan(0.5/[1])+atan((x+0.5)*0.5/[1]/sqrt([1]^2+(x+0.5)^2+0.5^2)))-(atan((x-0.5)/[1])+atan(0.5/[1])+atan((x-0.5)*0.5/[1]/sqrt([1]^2+(x-0.5)^2+0.5^2)))-(atan((x+0.5)/[1])+atan((-0.5)/[1])+atan((x+0.5)*(-0.5)/[1]/sqrt([1]^2+(x+0.5)^2+(-0.5)^2)))+(atan((x-0.5)/[1])+atan((-0.5)/[1])+atan((x-0.5)*(-0.5)/[1]/sqrt([1]^2+(x-0.5)^2+(-0.5)^2))))',-5,5)

i = 0
h[i].Fit(ffit)
p = h[i].ProjectionX("px",100,101)
p.Scale(0.5)
ffitx.SetParameters(ffit.GetParameter(0),ffit.GetParameter(1))
p.Draw()

llg = readfile('prof_lg.dat',[int,int,float,float])
lpwo = readfile('prof_pwo.dat',[int,int,float,float])
h0 = listofhisto('h0',['pwo','lg'],1000,-5.005,4.995,1000,-5.005,4.995,title=';x_{f};y_{f}')
for i,l in enumerate([lpwo,llg]):
  for x in l:
    a1,b1,a2,b2,c,d = x[0]+501, x[1]+501, -x[0]+501, -x[1]+501, x[2], x[3]
    if x[0]==x[1]==0: l2 = [(a1,b1)]
    elif x[0]==0: l2 = [(a1,b1),(a1,b2),(b1,a1),(b2,a1)]
    elif x[1]==0: l2 = [(a1,b1),(a2,b1),(b1,a1),(b1,a2)]
    elif x[0]==x[1]: l2 = [(a1,b1),(a1,b2),(a2,b1),(a2,b2)]
    else: l2 = [(a1,b1),(a1,b2),(a2,b1),(a2,b2),(b1,a1),(b2,a1),(b1,a2),(b2,a2)]
    for (u,v) in l2:
      j = h0[i].GetBin(u,v)
      _=h0[i].SetBinContent(j,c)
      _=h0[i].SetBinError(j,d)

p0 = [h0[i].ProjectionX("_px",501,501) for i in range(2)]


# fit time between tagger and hycal
def triangle(x,par):
  if x[0]-par[0]<par[3]: return -par[1]*abs(par[3])+par[2]
  elif x[0]-par[0]>par[4]: return -par[1]*abs(par[4])+par[2]
  else: return -par[1]*abs(x[0]-par[0])+par[2]

ffit = TF1('ffit',triangle,-900,-700,5)
ffit.SetParLimits(3,-40,-15)
ffit.SetParLimits(4,15,40)
ffit.SetParLimits(0,-820,-780)

runs = get_runs_between(889,979,'calib')
l = []
for x in runs:
  print x
  f = TFile('/work/hallb/prad/mlevilla/trees_calib/island/tree_island_'+str(x)+'.root')
  t = f.event
  l2 = []
  htemp = TH1F("htemp","",1000,-900,-700)
  for i in [1,2,5]:
    t.Draw("tg-tcl>>htemp","n_cl==1 && n_tag==1 && E>100 && abs(xhycal-xpos)<25 && abs(yhycal-ypos)<25 && trigger=="+str(i))
    ffit.SetParameter(0,-800)
    ffit.SetParameter(1,100)
    ffit.SetParameter(2,100)
    ffit.SetParameter(3,-25)
    ffit.SetParameter(4,25)
    htemp.Fit(ffit,'ww')
    l2.append([x,ffit.GetParameter(0),ffit.GetParameter(3),ffit.GetParameter(4)])
  l.append(l2)

# comparison between sector of LG
names = ['top','right','bottom','left']
f = [TFile(work+'/yields2/lg_'+x+'/cross_section_1431_1490.root') for x in names]
g = [x.theta.Get('gthexp_eff_ep') for x in f]
thbin = [x for x in g[0].GetX()]
meang = [sum(x.GetY()[i] for x in g)/len(g) for i in range(g[0].GetN())]
ratios = [[y/m for y,m in zip(x.GetY(),meang)] for x in g]
dratios = [[y/m for y,m in zip(x.GetEY(),meang)] for x in g]
colors = [ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kGreen+2]
gratios = [tgraph([w+0.03*i for w in thbin],x,dy=y,name='ratio_exp_eff_ep',title=';#theta (deg);(d#sigma/d#Omega)_{ep, sector}/mean',color=c,xrg=[3,6]) for x,y,c,i in zip(ratios,dratios,colors,range(4))]


# ratio of cross section on lg and pwo
thbin = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
rbin = [5636*tan(th*degrad) for th in thbin]
phibin = [0 if r<353.09 else acos(353.09/r) for r in rbin]
a = [4*(phi*r**2-353.09*r*sin(phi)) for r,phi in zip(rbin,phibin)]
a[-1] -= 4*(rbin[-1]-499.34)**2/2
rda = [(a[i+1]-a[i])/(pi*(rbin[i+1]**2-rbin[i]**2)) for i in range(len(a)-1)]
rda[-1] -= 4*(rbin[-1]-499.34)**2/2/(pi*(rbin[-1]**2-rbin[-2]**2))

names = ['crystal','lg']
f = [TFile(work+'/yields2/'+x+'/cross_section_1431_1490.root') for x in names]
g = [x.theta.Get('gthexp_eff_ep') for x in f]
theta = [0.706, 0.722, 0.738, 0.755, 0.774, 0.794, 0.815, 0.838, 0.864, 0.892, 0.922, 0.956, 0.993, 1.033, 1.079, 1.129, 1.182, 1.238, 1.301, 1.375, 1.463, 1.571, 1.707, 1.89, 2.102, 2.348, 2.637, 2.936, 3.236, 3.535, 3.834, 4.134, 4.433, 4.733, 5.036]
sigma = [[y for y in x.GetY()] for x in g]
dsigma = [[y for y in x.GetEY()] for x in g]
sigma_scale = [[sigma[i][j]/[1-rda[j],rda[j] if rda[j]!=0 else 1][i] for j in range(len(theta))] for i in range(2)]
ratio = [sigma_scale[1][j]/sigma_scale[0][j] for j in range(len(theta))]
gscale = [tgraph(theta,sigma_scale[i],dy=dsigma[i],title=';#theta (deg);d#sigma_{ep}/d#Omega') for i in range(2)]
gratio= tgraph(theta,ratio,dy=[ratio[j]*((dsigma[0][j]/sigma[0][j])**2+(dsigma[1][j]/sigma[1][j])**2)**0.5 if dsigma[1][j]!=0 else ratio[j]*dsigma[0][j]/sigma[0][j] for j in range(len(theta))])
gtheo = f[0].theta.Get('gom_2_ep_elas')


# get only interesting variables
fout = TFile('short_tree_0.root','recreate')
tout,vout = tree_init('T',[('x','f'),('y','f'),('z','f'),('E','f'),('theta','f'),('theta0','f'),('E0','f'),('Q2','f'),('Q20','f')],'',10)

f1 = TFile('/work/hallb/prad/mlevilla/ep_1GeV_15e4_theta04_75_FF5_CO1_2000_0.root')
f2 = TFile('/work/hallb/prad/mlevilla/ep_1GeV_15e4_theta04_75_FF5_CO1_2000_0_rec.root')
t1 = f1.T
t2 = f2.T

for i in progress(t1.GetEntries(),show=3):
  _=t1.GetEntry(i)
  _=t2.GetEntry(i)
  if not (getattr(t1,'GUN.N')==1 and getattr(t2,'GEM.N')==1 and getattr(t2,'HC.N')==1): continue
  if getattr(t2,'GEM.Z')[0]<0: continue
  vout['x'][0]=getattr(t2,'GEM.X')[0]
  vout['y'][0]=getattr(t2,'GEM.Y')[0]
  vout['z'][0]=getattr(t2,'GEM.Z')[0]+3000-89
  vout['E'][0]=getattr(t2,'HC.P')[0]
  vout['theta'][0]=atan2((vout['x'][0]**2+vout['y'][0]**2)**0.5,vout['z'][0])
  vout['theta0'][0]=getattr(t1,'GUN.Theta')[0]
  vout['E0'][0]=getattr(t1,'GUN.P')[0]
  vout['Q2'][0]=2*vout['E'][0]*1100.44/1e6*(1-cos(vout['theta'][0]))
  vout['Q20'][0]=2*vout['E0'][0]*1100.44/1e6*(1-cos(vout['theta0'][0]))
  _=tout.Fill()

fout.cd()
tout.Write()
fout.Close()


# fit resolution from snake scan
f = TFile(work+'/histo_calib/test2/reso_test2.root')
load_names()
lreso = []
for x in module_names.values():
  names = f.modules.Get(x).GetListOfKeys()
  if 'greso_'+x not in names: 
    lreso.append([x,0,0,0,0,0,0,0,0])
  else:
    ffit = TF1('ffit','sqrt(([0]/sqrt(x*0.001))^2+([1]/(x*0.001))^2+([2])^2)',250,1050)
    greso = f.modules.Get(x).Get('greso_'+x)
    greso.Fit(ffit,'r','',250,1050)
    [reso1,reso2,reso3] = [ffit.GetParameter(i) for i in range(3)]
    [dreso1,dreso2,dreso3] = [ffit.GetParError(i) for i in range(3)]
    reso = (reso1**2+reso2**2+reso3**2)**0.5
    dreso = 0
    if reso1!=0: reso+=(dreso1/reso1)**2
    if reso2!=0: reso+=(dreso2/reso2)**2
    if reso3!=0: reso+=(dreso3/reso3)**2
    if dreso>0: dreso = reso*dreso**0.5
    if reso<0.015 or reso>0.3: reso = 0
    lreso.append([x,reso1,reso2,reso3,dreso1,dreso2,dreso3,reso,dreso])

writelist('reso_test.txt',lreso,ro=5,ljust=8)


# geometrical acceptance of gem
from selection import *
from random import uniform
thbin = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
h = listofhisto('h',range(3),len(thbin)-1,xbin=thbin)
for i in progress(10000000,show=3):
  x = uniform(0,500)
  y = uniform(0,500)
  th = atan2((x**2+y**2)**0.5,5198)/degrad
  if x<50 and y<50: continue
  _=h[0].Fill(th)
  if gem_spacers(x,y): continue
  _=h[1].Fill(th)
    
h[2].Divide(h[1],h[0])
for i in range(len(thbin)-1):
  h[2].SetBinError(i+1,((h[1][i+1]+1)*(h[1][i+1]+2)/(h[0][i+1]+2)/(h[0][i+1]+3)-(h[1][i+1]+1)**2/(h[0][i+1]+2)**2)**0.5)


# combine 2 graphs
gStyle.SetOptFit(0)
f1 = TFile('/work/hallb/prad/mlevilla/simyields3/final_cross_section_1GeV.root')
f2 = TFile('/work/hallb/prad/mlevilla/simyields3/final_cross_section_2GeV.root')
g1 = f1.q2.Get('gq2_ge_ee3')
g2 = f2.q2.Get('gq2_ge_ee3')
x1 = [x for x in g1.GetX()]
x2 = [x for x in g2.GetX()]
y1 = [x for x in g1.GetY()]
y2 = [x for x in g2.GetY()]
ey1 = [x for x in g1.GetEY()]
ey2 = [x for x in g2.GetEY()]
x3 = x1+x2
y3 = y1+y2
ey3 = ey1+ey2
for i in range(len(x3)):
  xmin = min(x3[i:])
  imin = x3.index(xmin)
  x3[i],x3[imin] = x3[imin],x3[i]
  y3[i],y3[imin] = y3[imin],y3[i]
  ey3[i],ey3[imin] = ey3[imin],ey3[i]

g3 = tgraph(x3,y3,dy=ey3,title=';Q^{2} (GeV^{2});G_{E}')
ffit_q2 = TF1('ffit_q2','[0]-[1]^2*x/6/0.197327^2',0,0.05)
g3.Fit(ffit_q2)


# slices of E vs theta
f1 = TFile('/work/hallb/prad/mlevilla/yields5/cross_section_1443_1490.root')
h12 = f1.selection.Get('h2cut1_prod_ep_1')
h1 = [h12.ProjectionY('h1_py'+str(i),10*i,10*(i+1)) for i in range(10)]
for x in h1: 
  x.SetFillStyle(0)
  x.SetLineColor(ROOT.kRed)

f2 = TFile('/work/hallb/prad/mlevilla/simyields5/cross_section_sim_2GeV_m.root')
h22 = f2.selection.Get('h2cut1_ep_1')
h2 = [h22.ProjectionY('h2_py'+str(i),10*i,10*(i+1)) for i in range(10)]
for i,x in enumerate(h2): 
  x.SetFillStyle(0)
  x.SetLineColor(ROOT.kBlue)
  if h2[i].Integral()!=0: x.Scale(h1[i].Integral()/h2[i].Integral())


# gun theta distributions
fout = TFile('gun_theta_distribution.root','recreate')
names = [['/work/hallb/prad/mlevilla/new_evgen_files/ep/1GeV/ep_1GeV_1e6_th04_75_CO1_1000_FF5','/work/hallb/prad/mlevilla/new_evgen_files/ep/2GeV/ep_2GeV_1e6_th04_75_CO1_2000_FF5'],['/work/hallb/prad/mlevilla/new_evgen_files/moller/1GeV/moller_1GeV_1e6_th04_75_CO1_1000','/work/hallb/prad/mlevilla/new_evgen_files/moller/2GeV/moller_2GeV_1e6_th04_35_CO1_2000']]

h = listofhisto('hth',[['ep','ee'],['1GeV','2GeV'],['elas','rad']],10000,0,10)
h2 = listofhisto('hEth',[['ep','ee'],['1GeV','2GeV'],['elas','rad']],1000,0,10,1000,0,2500)

for i1 in range(2):
  for i2 in range(2):
    for i3 in range(300):
      print i1,i2,i3
      name = names[i1][i2]+'_'+str(i3)+'_e-.dat' if i1==0 else names[i1][i2]+'_'+str(i3)+'.dat'
      [l_E_e1,l_th_e1,l_E_e2,l_th_e2,l_E_g,l_th_g] = readfile(name,[float]*6,cols=[0,1,3,4,6,7],tr=True)
      for j in range(1000000):
        if l_E_e2[j]!=l_E_e2[j]: continue
        for k in range(2):
          _=h[i1][i2][k].Fill(l_th_e1[j]/degrad)
          _=h2[i1][i2][k].Fill(l_th_e1[j]/degrad,l_E_e1[j])
          if i1==1: 
            _=h[i1][i2][k].Fill(l_th_e2[j]/degrad)
            _=h[i1][i2][k].Fill(l_th_e2[j]/degrad,l_E_e2[j])
        if l_E_g[j]>0: 
          _=h[i1][i2][1].Fill(l_th_g[j]/degrad)
          _=h[i1][i2][1].Fill(l_th_g[j]/degrad,l_E_g[j])

# center
f = TFile(workf+'/yields5/yields_1163.root')
t = f.clusters
pair = 1
h = listofhisto('h',['x','y'],200,-50,50)

for i,_ in enumerate(progress(t,n=t.GetEntries())):
  if t.event!=3: continue
  if t.theta_gem[0]<0.7 or t.theta_gem[0]>3.5 or t.theta_gem[1]<0.7 or t.theta_gem[1]>3.5: continue
  if pair:
    x1,y1,x2,y2 = t.x_gem[0],t.y_gem[0],t.x_gem[1],t.y_gem[1]
    pair = 0
  else:
    x3,y3,x4,y4 = t.x_gem[0],t.y_gem[0],t.x_gem[1],t.y_gem[1]
    pair = 1
    _=h[0].Fill(((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)))
    _=h[1].Fill(((x1*y2-y1*x2)*(y3-x4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)))
    


# phi dependency
runs = get_runs_between(1443,1490,'prod')
fout = TFile('temp.root','recreate')
h2 = listofhisto('h2',range(34),16,-180,180)
for j,r in enumerate(runs):
  f = TFile(workf+'/yields4/yields_'+str(r)+'.root')
  t = f.clusters
  print r
  for _ in progress(t,n=t.GetEntries(),show=3):
    if t.event!=0: continue
    if t.theta_hycal[0]<0.8 or t.theta_hycal[0]>=4.2: continue
    phi = atan2(t.y_hycal[0],t.x_hycal[0])*180/pi
    i = int((t.theta_hycal[0]-0.8)/0.1)
    _=h2[i].Fill(phi)


for i in range(17):
  h2[i].Scale(16./h2[i].GetEntries())


# event generator yields
dirw = '/work/hallb/prad/xiongw/simulation/PRadSim/farm_script/input/'
dirm = '/work/hallb/prad/mlevilla/new_evgen_files/ep3/2GeV/'

thbin = [0.4+0.01*i for i in range(30)]+[0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]+[5.7,6.2,6.7,7.2,7.7]
hth = listofhisto('hth',[['m','w','r'],['e','g','all']],len(thbin)-1,xbin=thbin)

for i in range(1,151):
  print i
  lw = readfile(dirw+'ep_2GeV_1e6_theta_0.5_7.5_FF5_CO_2000_JL_'+str(i)+'_e-.dat',[float,float],cols=[1,7])
  lm = readfile(dirm+'ep_2GeV_1e6_th05_75_CO1_2000_FF5_'+str(i)+'_e-.dat',[float,float],cols=[1,7])
  for x in lm: 
    _ = hth[0][0].Fill(x[0]/degrad)
    _ = hth[0][1].Fill(x[1]/degrad)
    _ = hth[0][2].Fill(x[0]/degrad)
    _ = hth[0][2].Fill(x[1]/degrad)
  for x in lw: 
    _ = hth[1][0].Fill(x[0]/degrad)
    _ = hth[1][1].Fill(x[1]/degrad)
    _ = hth[1][2].Fill(x[0]/degrad)
    _ = hth[1][2].Fill(x[1]/degrad)

for i in range(3):
  divide_binwidth(hth[0][i])
  divide_binwidth(hth[1][i])
  hth[2][i].Divide(hth[0][i],hth[1][i])

dirw = '/work/hallb/prad/xiongw/simulation/PRadSim/output/'
dirm = '/work/hallb/prad/mlevilla/new_simulation_files/ep2/2GeV/'

hth2 = listofhisto('hth2',['m','w','r'],len(thbin)-1,xbin=thbin)
cm = get_chain('T',dirm+'ep_2GeV_1e5_th05_75_CO1_2000_FF5',900)
cw = get_chain('T',dirw+'ep_2GeV_1e6_theta_0.5_7.5_FF5_CO_2000_JL',90,start=1,ext='_e-.root')

cm.Draw('GUN.Theta*180/pi>>hth2_m')
cw.Draw('GUN.Theta*180/pi>>hth2_w')

from esepp import *
theo0 = [[born((thbin[i]+(thbin[i+1]-thbin[i])*0.01*j)*degrad,2.142,flag=5) for j in range(101)] for i in range(len(thbin)-1)] 
theo1 = [sum(theo0[i][j]*sin((thbin[i]+(thbin[i+1]-thbin[i])*0.01*j)*degrad)*(thbin[i+1]-thbin[i])*0.01 for j in range(101))/(cos(thbin[i]*degrad)-cos(thbin[i+1]*degrad)) for i in range(len(thbin)-1)]

f = TFile(dataf+'/dat_differences.root')
hsim = f.hth_m_e
multiply_binwidth(hsim)
thbin = get_binning(hsim)
for i in range(len(thbin)-1):
  err = hsim.GetBinError(i+1)
  hsim[i+1] = hsim[i+1]/(cos(thbin[i]*degrad)-cos(thbin[i+1]*degrad))
  hsim.SetBinError(i+1,err/(cos(thbin[i]*degrad)-cos(thbin[i+1]*degrad)))
  
l = root_to_list(hsim)
ratio = [l[1][i]/theo1[i] for i in range(70)]
dratio = [l[3][i]/theo1[i] for i in range(70)]
g3 = tgraph(l[0],ratio,[],dratio)


# looking in at spread of uncertainty in esepp
dirm = '/work/hallb/prad/mlevilla/new_evgen_files/ep3/2GeV/'

thbin = [0.4+0.01*i for i in range(30)]+[0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]+[5.7,6.2,6.7,7.2,7.7]
nbin = len(thbin)-1
nf = 300
hth = listofhisto('hth',[['e','g'],range(nf)],len(thbin)-1,xbin=thbin)
for i in range(nf):
  print i
  lm = readfile(dirm+'ep_2GeV_1e6_th05_75_CO1_2000_FF5_'+str(i)+'_e-.dat',[float,float],cols=[1,7])
  for x in lm:
    _ = hth[0][i].Fill(x[0]/degrad)
    _ = hth[1][i].Fill(x[1]/degrad)

l0 = [[hth[0][j][k+1] for k in range(nbin)] for j in range(nf)]
lsum = [sum(l0[j][k] for j in range(nf)) for k in range(nbin)]
lmean = [sum(l0[j][k] for j in range(nf))/nf for k in range(nbin)]
l1 = [[l0[j][k]/lmean[k]-1. if lmean[k]!=0 else 0. for k in range(nbin)] for j in range(nf)] 
lp = [[(l0[j][k]-lmean[k])/lmean[k]**0.5 if lmean[k]>0 else 0. for k in range(nbin)] for j in range(nf)] 
hs = listofhisto('hs',range(nbin),100,-10,10)
for j in range(nf):
  for k in range(nbin):
    _=hs[k].Fill(lp[j][k])

lsigma = [0 for i in range(nbin)]
lerrsigma = [0 for i in range(nbin)]
for k in range(10,nbin):
  hs[k].Fit('gaus')
  g = hs[k].GetFunction('gaus')
  lsigma[k] = g.GetParameter(2)
  lerrsigma[k] = g.GetParError(2)

rsyst = [((max(1,lsigma[k])+lerrsigma[k])**2-1)**0.5 for k in range(nbin)] 
error_tot = [(1+rsyst[k]**2)**0.5 for k in range(nbin)]


# difference reso between snake and physics
lw = readfile(workf+'/resolution_curve_eSigma4_ex1.txt',[float,float,float],cols=range(1,4))
lm = readfile(workf+'/histo_calib/test3/reso_2gaus_5_3.txt',[float]*8,cols=range(1,9))
n = 210

lgm = [[(lm[i][0]**2/E+lm[i][1]**2/E**2+lm[i][2]**2)**0.5 for E in [0.1+0.01*j for j in range(n)]] for i in range(1728)]
lgw = [[(lw[i][0]**2/E+lw[i][1]**2/E**2+lw[i][2]**2)**0.5 for E in [0.1+0.01*j for j in range(n)]] for i in range(1728)]
lgm_max = [[((lm[i][0]+lm[i][3])**2/E+(lm[i][1]+lm[i][4])**2/E**2+(lm[i][2]+lm[i][5])**2)**0.5 for E in [0.1+0.01*j for j in range(n)]] for i in range(1728)]
lgm_min = [[((lm[i][0]-lm[i][3])**2/E+(lm[i][1]-lm[i][4])**2/E**2+(lm[i][2]-lm[i][5])**2)**0.5 for E in [0.1+0.01*j for j in range(n)]] for i in range(1728)]
lgm_err = [[0.5*(lgm_max[i][j]-lgm_min[i][j]) for j in range(n)] for i in range(1728)]
lchi2 = [(lgm[i][j]-lgw[i][j])**2/lgm_err[i][j]**2 if lgm_err[i][j]!=0 for j in range(n)] for i in range(1728)]
lgdiff = [[lgm[i][j]-lgw[i][j] for j in range(n)] for i in range(1728)]

lchi2 = [[(lgm[i][j]-lgw[i][j])**2/lgm[i][j]**2 if lgm[i][j]!=0 else 0 for j in range(220)] for i in range(1728)]
lschi2 = [sum(lchi2[i]) for i in range(1728)]
lschi2_ndf = [lschi2[i]/n for i in range(1728)]
lratio = [lschi2_ndf[i]**0.5 for i in range(1728)]

lfm = [TF1('fm'+str(i),'sqrt(('+str(lm[i][0])+'/sqrt(x*0.001))^2 + ('+str(lm[i][1])+'/(x*0.001))^2 + '+str(lm[i][2])+'^2)',0,2200) for i in range(1728)]
lfw = [TF1('fw'+str(i),'sqrt(('+str(lw[i][0])+'/sqrt(x*0.001))^2 + ('+str(lw[i][1])+'/(x*0.001))^2 + '+str(lw[i][2])+'^2)',0,2200) for i in range(1728)]
for i in range(1728): lfw[i].SetLineColor(ROOT.kBlue)

load_names()
lnames = module_names.keys()
g0 = tgraph(lnames,lratio,title=';#module;#sqrt{#frac{1}{ndf}#sum (snake-physics)^{2}/snake^{2}}')

syst_lg = [sum([ldiff[i][j] for i in range(576) if (lm[i]!=[0 for k in range(8)] and lw[i]!=[0 for k in range(3)] and ldiff[i][j]==ldiff[i][j])])/len([ldiff[i][j] for i in range(576) if (lm[i]!=[0 for k in range(8)] and lw[i]!=[0 for k in range(3)] and ldiff[i][j]==ldiff[i][j])]) for j in range(n)]
syst_pwo = [sum([ldiff[i][j] for i in range(576,1728) if (lm[i]!=[0 for k in range(8)] and lw[i]!=[0 for k in range(3)] and ldiff[i][j]==ldiff[i][j])])/len([ldiff[i][j] for i in range(576) if (lm[i]!=[0 for k in range(8)] and lw[i]!=[0 for k in range(3)] and ldiff[i][j]==ldiff[i][j])]) for j in range(n)]

fpwo = TF1("fpwo",'sqrt((0.2368/sqrt(x*0.001))^2)',0,2300)
flg = TF1("flg",'sqrt((0.0391/sqrt(x*0.001))^2 + 0.03264^2)',0,2300)

ratio_syst_lg = [syst_lg[j]/flg.Eval(100+10*j) for j in range(n)]
ratio_syst_pwo = [syst_pwo[j]/fpwo.Eval(100+10*j) for j in range(n)]
gratio_syst_lg = tgraph([100+10*i for i in range(210)],ratio_syst_lg,title=';E (MeV);#sigma_{syst}/reso')
gratio_syst_pwo = tgraph([100+10*i for i in range(210)],ratio_syst_pwo,title=';E (MeV);#sigma_{syst}/reso')


# merge simyields
folder = workf+'/simyields/'
name = 'simyields_2GeV_bla_'
for i in range(10):
  x = os.system('hadd -f '+folder+name+'m_'+str(i)+'.root '+' '.join([folder+name+str(j)+'.root' for j in range(30*i,30*(i+1))]))
  print x 

for i in range(300): os.remove(folder+name+str(i)+'.root')

# stability through ratio
# 2GeV
rprod2 = get_runs_between(1362,1516,'prod')
fprod2 = [TFile(workf+'/yields/yields_'+str(i)+'.root') for i in rprod2]
hprod2 = [[x.yields.Get('hth_'+y) for x in fprod2] for y in ['ep','ee1','ee2','ee3']]
iprod2 = [[sum(x[i] for i in range(7,25)) for x in y] for y in hprod2]
raprod2 = [[iprod2[0][j]/iprod2[i+1][j] for j in range(len(iprod2[0]))] for i in range(3)]
raprod2_err = [[raprod2[i][j]*(1/iprod2[0][j]+1/iprod2[i+1][j])**0.5 for j in range(len(iprod2[0]))] for i in range(3)]
gprod2 = [tgraph(rprod2,raprod2[i],[],raprod2_err[i],name='gprod_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
mraprod2 = [sum(x)/len(x) for x in raprod2]
nraprod2 = [[y1/x2 for y1 in x1] for x1,x2 in zip(raprod2,mraprod2)]
nraprod2_err = [[nraprod2[i][j]*(1/iprod2[0][j]+1/iprod2[i+1][j])**0.5 for j in range(len(iprod2[0]))] for i in range(3)]
gnprod2 = [tgraph(rprod2,nraprod2[i],[],nraprod2_err[i],name='gprod_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]

rempty2 = get_runs_between(1362,1516,'empty')
fempty2 = [TFile(workf+'/yields/yields_'+str(i)+'.root') for i in rempty2]
hempty2 = [[x.yields.Get('hth_'+y) for x in fempty2] for y in ['ep','ee1','ee2','ee3']]
iempty2 = [[sum(x[i] for i in range(7,25)) for x in y] for y in hempty2]
raempty2 = [[iempty2[0][j]/iempty2[i+1][j] for j in range(len(iempty2[0]))] for i in range(3)]
raempty2_err = [[raempty2[i][j]*(1/iempty2[0][j]+1/iempty2[i+1][j])**0.5 for j in range(len(iempty2[0]))] for i in range(3)]
gempty2 = [tgraph(rempty2,raempty2[i],[],raempty2_err[i],name='gempty_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee',color=[ROOT.kBlack,ROOT.kBlack,ROOT.kBlue]) for i in range(3)]
mraempty2 = [sum(x)/len(x) for x in raempty2]
nraempty2 = [[y1/x2 for y1 in x1] for x1,x2 in zip(raempty2,mraempty2)]
nraempty2_err = [[nraempty2[i][j]*(1/iempty2[0][j]+1/iempty2[i+1][j])**0.5 for j in range(len(iempty2[0]))] for i in range(3)]
gnempty2 = [tgraph(rempty2,nraempty2[i],[],nraempty2_err[i],name='gempty_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee',color=[ROOT.kBlack,ROOT.kBlack,ROOT.kBlue]) for i in range(3)]

# 1GeV
rprod1 = get_runs_between(1059,1345,'prod')
fprod1 = [TFile(workf+'/yields5/yields_'+str(i)+'_hfid.root') for i in rprod1]
hprod1 = [[x.yields.Get('hth_'+y) for x in fprod1] for y in ['ep','ee1','ee2','ee3']]
iprod1 = [[sum(x[i] for i in range(17,28)) for x in y] for y in hprod1]
raprod1 = [[iprod1[0][j]/iprod1[i+1][j] for j in range(len(iprod1[0]))] for i in range(3)]
raprod1_err = [[raprod1[i][j]*(1/iprod1[0][j]+1/iprod1[i+1][j])**0.5 for j in range(len(iprod1[0]))] for i in range(3)]
gprod1 = [tgraph(rprod1,raprod1[i],[],raprod1_err[i],name='gprod_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
mraprod1 = [sum(x)/len(x) for x in raprod1]
nraprod1 = [[y1/x2 for y1 in x1] for x1,x2 in zip(raprod1,mraprod1)]
nraprod1_err = [[nraprod1[i][j]*(1/iprod1[0][j]+1/iprod1[i+1][j])**0.5 for j in range(len(iprod1[0]))] for i in range(3)]
gnprod1 = [tgraph(rprod1,nraprod1[i],[],nraprod1_err[i],name='gprod_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]

rempty1 = get_runs_between(1059,1345,'empty')
fempty1 = [TFile(workf+'/yields5/yields_'+str(i)+'_hfid.root') for i in rempty1]
hempty1 = [[x.yields.Get('hth_'+y) for x in fempty1] for y in ['ep','ee1','ee2','ee3']]
iempty1 = [[sum(x[i] for i in range(17,28)) for x in y] for y in hempty1]
raempty1 = [[iempty1[0][j]/iempty1[i+1][j] for j in range(len(iempty1[0]))] for i in range(3)]
raempty1_err = [[raempty1[i][j]*(1/iempty1[0][j]+1/iempty1[i+1][j])**0.5 for j in range(len(iempty1[0]))] for i in range(3)]
gempty1 = [tgraph(rempty1,raempty1[i],[],raempty1_err[i],name='gempty_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee',color=[ROOT.kBlack,ROOT.kBlack,ROOT.kBlue]) for i in range(3)]
mraempty1 = [sum(x)/len(x) for x in raempty1]
nraempty1 = [[y1/x2 for y1 in x1] for x1,x2 in zip(raempty1,mraempty1)]
nraempty1_err = [[nraempty1[i][j]*(1/iempty1[0][j]+1/iempty1[i+1][j])**0.5 for j in range(len(iempty1[0]))] for i in range(3)]
gnempty1 = [tgraph(rempty1,nraempty1[i],[],nraempty1_err[i],name='gempty_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee',color=[ROOT.kBlack,ROOT.kBlack,ROOT.kBlue]) for i in range(3)]


# alignment with moller
f = TFile(workf+'/shtree_1301.root')
t = f.ee2
even = True
h = listofhisto('h',['x','y'],1000,-50,50)

for i,_ in progress(enumerate(t),n=t.GetEntries(),show=3):
  if any((th<0.9 or th>2.5) for th in t.theta_hycal): continue
  if even:
    x1,x2,y1,y2 = t.x_hycal[0],t.x_hycal[1],t.y_hycal[0],t.y_hycal[1]
    even = False
  else:
    x3,x4,y3,y4 = t.x_hycal[0],t.x_hycal[1],t.y_hycal[0],t.y_hycal[1]
    _=h[0].Fill(((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)))
    _=h[1].Fill(((x1*y2-y1*x2)*(y3-x4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)))


# comparing energy spectrum
runs = [1105,1136,1155]
fout = TFile('temp.root','recreate')
h = listofhisto('hE',runs,325,0.,1300,title=';E (MeV);#events')
for i,x in enumerate(runs):
  f = TFile(workf+'/yields5/yields_'+str(x)+'_hfid.root','theta_gem>1.0 && theta_gem<3.5')
  fout.cd()
  f.clusters.Draw('E*1000>>hE_'+str(x))

for i,x in enumerate(h):
  x.Scale(1./x.Integral(262,290))
  x.SetLineColor([ROOT.kGreen+2,ROOT.kRed,ROOT.kBlue][i])


# stability through ratio with bg subtraction
load_live_charge()

r2_all = get_runs_between(1362,1514,'all')
r2 = get_runs_between(1362,1514,'prod')
f2 = {i:TFile(workf+'/yields/yields_'+str(i)+'.root') for i in r2_all}
h2 = {i:[f2[i].yields.Get('hth_'+x) for x in ['ep','ee1','ee2','ee3']] for i in r2_all}
irange = [12,25]
i2 = {i:[sum(x[j] for j in range(*irange)) for x in h2[i]] for i in r2_all}
i22 = [[i2[i][j]-sum(i2[k][j] for k in get_bg_runs(i))*live_charge[i]/sum(live_charge[k] for k in get_bg_runs(i)) for j in range(4)] for i in r2]
prod2 = [[i2[i][j] for j in range(4)] for i in r2]
bg2 = [[sum(i2[k][j] for k in get_bg_runs(i))*live_charge[i]/sum(live_charge[k] for k in get_bg_runs(i) ) for j in range(4)] for i in r2]
ratiobg2 = [[bg2[i][j]/prod2[i][j] for j in range(len(prod2[i]))] for i in range(len(prod2))]
ra2 = [[i22[j][0]/i22[j][i+1] for j in range(len(i22))] for i in range(3)]
ra2_err = [[ra2[i][j]*(1/i22[j][0]+1/i22[j][i+1])**0.5 for j in range(len(i22))] for i in range(3)]
g2 = [tgraph(r2,ra2[i],[],ra2_err[i],name='gprod_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
mra2 = [sum(x)/len(x) for x in ra2]
nra2 = [[y1/x2 for y1 in x1] for x1,x2 in zip(ra2,mra2)]
nra2_err = [[nra2[i][j]*(1/i22[j][0]+1/i22[j][i+1])**0.5 for j in range(len(i22))] for i in range(3)]
gn2 = [tgraph(r2,nra2[i],[],nra2_err[i],name='gprod_2GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
p2 = [[0,10],[10,67],[67,74]]
p2i = [1390,1500]
gp2 = [tgraph([p2i[i],p2i[i]],[0,10]) for i in range(len(p2i))]
av20 = sum(ra2[2][i]/ra2_err[2][i]**2 for i in range(len(r2)))/sum(1/ra2_err[2][i]**2 for i in range(len(r2)))
av2 = [sum(nra2[2][i]/nra2_err[2][i]**2 for i in range(*p2[j]))/sum(1/nra2_err[2][i]**2 for i in range(*p2[j])) for j in range(len(p2))]
av2_err = [1/sum(1/nra2_err[2][i]**2 for i in range(*p2[j]))**0.5 for j in range(len(p2))]
av2x =  [sum(r2[i]/nra2_err[2][i]**2 for i in range(*p2[j]))/sum(1/nra2_err[2][i]**2 for i in range(*p2[j])) for j in range(len(p2))]
gav2 = tgraph(av2x,av2,dy=av2_err,name='gav_2GeV_ep_ee3',title=';#run;ep/ee')
pulls2 = TH1F('pulls2','',80,-10,10)
for i in range(len(r2)): _ = pulls2.Fill((ra2[2][i]-av20)/ra2_err[2][i])

r1_all = get_runs_between(1059,1341,'all')
r1 = get_runs_between(1059,1341,'prod')
r1.remove(1243)
r1_all.remove(1243)
f1 = {i:TFile(workf+'/yields/yields_'+str(i)+'.root') for i in r1_all}
h1 = {i:[f1[i].yields.Get('hth_'+x) for x in ['ep','ee1','ee2','ee3']] for i in r1_all}
i1 = {i:[sum(x[j] for j in range(20,28)) for x in h1[i]] for i in r1_all}
i12 = [[i1[i][j]-sum(i1[k][j] for k in get_bg_runs(i))*live_charge[i]/sum(live_charge[k] for k in get_bg_runs(i)) for j in range(4)] for i in r1]
prod1 = [[i1[i][j] for j in range(4)] for i in r1]
bg1 = [[sum(i1[k][j] for k in get_bg_runs(i))*live_charge[i]/sum(live_charge[k] for k in get_bg_runs(i) ) for j in range(4)] for i in r1]
ratiobg1 = [[bg1[i][j]/prod1[i][j] for j in range(len(prod1[i]))] for i in range(len(prod1))]
ra1 = [[i12[j][0]/i12[j][i+1] for j in range(len(i12))] for i in range(3)]
ra1_err = [[ra1[i][j]*(1/i12[j][0]+1/i12[j][i+1])**0.5 for j in range(len(i12))] for i in range(3)]
g1 = [tgraph(r1,ra1[i],[],ra1_err[i],name='gprod_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
mra1 = [sum(x)/len(x) for x in ra1]
nra1 = [[y1/x2 for y1 in x1] for x1,x2 in zip(ra1,mra1)]
nra1_err = [[nra1[i][j]*(1/i12[j][0]+1/i12[j][i+1])**0.5 for j in range(len(i12))] for i in range(3)]
gn1 = [tgraph(r1,nra1[i],[],nra1_err[i],name='gprod_1GeV_ep_ee'+str(i+1),title=';#run;ep/ee') for i in range(3)]
#p1 = [[0,2],[2,41],[41,75]]
p1 = [[0,2],[2,40],[40,74]]
p1i = [1065,1265]
gp1 = [tgraph([p1i[i],p1i[i]],[0,10]) for i in range(len(p1i))]
av10 = sum(ra1[2][i]/ra1_err[2][i]**2 for i in range(len(r1[2:])))/sum(1/ra1_err[2][i]**2 for i in range(len(r1[2:])))
av1 = [sum(nra1[2][i]/nra1_err[2][i]**2 for i in range(*p1[j]))/sum(1/nra1_err[2][i]**2 for i in range(*p1[j])) for j in range(len(p1))]
av1_err = [1/sum(1/nra1_err[2][i]**2 for i in range(*p1[j]))**0.5 for j in range(len(p1))]
av1x =  [sum(r1[i]/nra1_err[2][i]**2 for i in range(*p1[j]))/sum(1/nra1_err[2][i]**2 for i in range(*p1[j])) for j in range(len(p1))]
gav1 = tgraph(av1x,av1,dy=av1_err,name='gav_1GeV_ep_ee3',title=';#run;ep/ee')
pulls1 = TH1F('pulls1','',50,-10,10)
for i in range(len(r1[2:])): _ = pulls1.Fill((ra1[2][i]-av10)/ra1_err[2][i])


# submit extract_yields
suffix = ''
os.system('exec/extract_yields')


# energy spectrum  from trees
theta1 = [[[0.7,0.8],[2.0,2.2],[3.4,3.6],[4.5,5.2]],[[0.7,0.8],[2.0,2.2],[3.4,3.6],[4.2,4.5]]]
theta2 = [[[0.7,0.8],[2.0,2.2],[3.4,3.6],[4.5,5.2]],[[0.7,0.8],[1.0,1.1],[1.6,1.7],[2.4,2.5]]]
energy2 = [[[1.8,2.4],[1.8,2.4],[1.8,2.4],[1.8,2.4]],[[1.,2.],[0.8,1.8],[0.4,1.1],[0.2,0.7]]]
energy1 = [[[0.8,1.3],[0.8,1.3],[0.8,1.3],[0.8,1.3]],[[0.6,1.2],[0.2,0.7],[0.05,0.4],[0.,0.4]]]

fout = TFile('temp.root','recreate')
texp = get_chain('clusters',workf+'/yields5/yields',get_runs_between(1443,1453,'prod'))
texp_empty = get_chain('clusters',workf+'/yields5/yields',get_runs_between(1431,1455,'empty'))
tsim = get_chain('clusters',workf+'/simyields5/simyields_2GeV_m',range(2))

hexp = listofhisto('hexp',[['ep','ee1'],range(len(theta2[0]))],400,[[x[0] for x in y] for y in energy2],[[x[1] for x in y] for y in energy2],title=';E (GeV);')
hexp_empty = listofhisto('hexp_empty',[['ep','ee1'],range(len(theta2[0]))],400,[[x[0] for x in y] for y in energy2],[[x[1] for x in y] for y in energy2],title=';E (GeV);')
hexp_diff = listofhisto('hexp_diff',[['ep','ee1'],range(len(theta2[0]))],400,[[x[0] for x in y] for y in energy2],[[x[1] for x in y] for y in energy2],title=';E (GeV);')
hsim = listofhisto('hsim',[['ep','ee1'],range(len(theta2[0]))],400,[[x[0] for x in y] for y in energy2],[[x[1] for x in y] for y in energy2],title=';E (GeV);')

load_live_charge()
ratio_live_charge = sum(live_charge[run] for run in get_runs_between(1443,1453,'prod'))/sum(live_charge[run] for run in get_runs_between(1431,1455,'empty'))

for i in range(len(theta2[0])):
  texp.Draw('E>>hexp_ep_'+str(i),'event==0 && theta_gem>'+str(theta2[0][i][0])+' && theta_gem<'+str(theta2[0][i][1]))
  texp_empty.Draw('E>>hexp_empty_ep_'+str(i),'event==0 && theta_gem>'+str(theta2[0][i][0])+' && theta_gem<'+str(theta2[0][i][1]))
  tsim.Draw('E>>hsim_ep_'+str(i),'event==0 && theta_gem>'+str(theta2[0][i][0])+' && theta_gem<'+str(theta2[0][i][1]))
  texp.Draw('E>>hexp_ee1_'+str(i),'event==1 && theta_gem>'+str(theta2[1][i][0])+' && theta_gem<'+str(theta2[1][i][1]))
  texp_empty.Draw('E>>hexp_empty_ee1_'+str(i),'event==1 && theta_gem>'+str(theta2[1][i][0])+' && theta_gem<'+str(theta2[1][i][1]))
  tsim.Draw('E>>hsim_ee1_'+str(i),'event==1 && theta_gem>'+str(theta2[1][i][0])+' && theta_gem<'+str(theta2[1][i][1]))
  for j in range(2):
    hexp_diff[j][i].Add(hexp[j][i],hexp_empty[j][i],-ratio_live_charge)
    hexp_diff[j][i].SetLineColor(ROOT.kRed)
    hsim[j][i].SetLineColor(ROOT.kBlue)

# fit of trigger efficiency

f = TFile('/work/hallb/prad/mlevilla/histo_calib/test2_3/eff_3.root')
gpwo = f.regions.Get('gcenter_pwo')
glg = f.regions.Get('gcenter_lg')

ffit = TF1('ffit','[0]*(1-exp(-[1]*x+[2]))',200,1200)
ffit.SetParameters(0.995,0.1,0.1)
ffit.SetParLimits(0,0.994,1.0)
ffit.SetParLimits(1,0.0001,10)
ffit.SetParLimits(2,-10,0)

# map of trigger efficiency
f = TFile('/work/hallb/prad/mlevilla/histo_calib/test2_3/eff_3.root')
g = f.gmean_eff1
from viewer import *
v = hcviewer('eff',False,0.95,1.001)
v.fill([[x for x in g.GetY()]],[[module_names[int(x)] for x in g.GetX()]])
v.draw()

#[9, 259, 2, 182, 281, 4, 1, 0, 288, 276, 291, 295, 298, 3, 293, 280, 289, 275, 5, 282, 286, 120]

# relaunch fail 
l = [x for x in os.listdir('/u/home/mlevilla/.farm_out') if '.err' in x]
lep1,lep2,lmoller1,lmoller2 = [],[],[],[]

for x in l:
 x = x.replace('err','out')
 f = open('/u/home/mlevilla/.farm_out/'+x)
 l2 = f.readlines()
 i = l2[0].index('.mac')
 j = l2[0][:i].rindex('_')
 number = int(l2[0][j+1:i])
 if 'ep' in l2[0] and '1GeV' in l2[0]: lep1.append(number)
 elif 'ep' in l2[0] and '2GeV' in l2[0]: lep2.append(number)
 elif 'moller' in l2[0] and '1GeV' in l2[0]: lmoller1.append(number)
 elif 'moller' in l2[0] and '2GeV' in l2[0]: lmoller2.append(number)


for a in ['ep','moller']:
  for b in ['1','2']:
    folder = workf+'/new_simulation_files/'+a+'4/'+b+'GeV/'
    l = os.listdir(folder)
    l2 = []
    for x in l:
      j = x.index('.root')
      i = x.rindex('_')
      l2.append(int(x[i+1:j]))
 
    l3 = list(set(range(300))-set(l2))
    l3.sort()
    print a,b, ' '.join([str(x) for x in l3])

#foreach x ()
#pradsim.py -evtype ep -file $WORK/new_evgen_files/ep3/2GeV/ep_2GeV_1e6_th05_75_CO1_2000_FF5_${x} -name ep_2GeV_1e6_th05_75_CO1_2000_FF5_${x} -n 1000000 -nf 1 -o $WORK/new_simulation_files/ep4/2GeV -b -wt 1000

# launching pradsim 
#pradsim.py -evtype ep -file $WORK/new_evgen_files/ep3/2GeV/ep_2GeV_1e6_th05_75_CO1_2000_FF5 -name ep_2GeV_1e6_th05_75_CO1_2000_FF5 -n 1000000 -nf 300 -o $WORK/new_simulation_files/ep4/2GeV -b -wt 1000

# theta resolution from simulation
f = TFile(workf+'/ep_2GeV_1e6_th05_75_CO1_2000_FF5_0_test_rec.root')
t = f.T
h = []
for i in range(10):
  t.Draw('(atan2(sqrt(GEM.X[0]^2+GEM.Y[0]^2),GEM.Z[0]-89)-GUN.Theta[0])*180/pi>>h'+str(i)+'(1000,-0.1,0.1)','GUN.Theta[0]*180/pi>0.5*'+str(i+1)+' && GUN.Theta[0]*180/pi<0.5*'+str(i+2))
  h.append(f.Get('h'+str(i)))

t.Draw('(atan2(sqrt(GEM.X[0]^2+GEM.Y[0]^2),GEM.Z[0]-89)-GUN.Theta[0])*180/pi:GUN.Theta[0]*180/pi>>h2d(70,0.5,7.,1000,-0.1,0.1)','','colz')

f = TFile(workf+'/new_simulation_files/ep4/2GeV/ep_2GeV_1e6_th05_75_CO1_2000_FF5_0.root')
t = f.T
for i in range(10):
  t.Draw('(atan2(sqrt(GEM.X[0]^2+GEM.Y[0]^2),GEM.Z[0]-(-3000+88.9))-GUN.Theta[0])*180/pi>>h'+str(i)+'(1000,-0.1,0.1)','GUN.Theta[0]*180/pi>0.5*'+str(i+1)+' && GUN.Theta[0]*180/pi<0.5*'+str(i+2))


# looking at profiles
runs = [1443,1444,1445,1446]
f = [TFile(workf+'/trees_sel6/tree_island_'+str(x)+'.root') for x in runs]
ts = [x.event for x in f]
load_positions()
imod = [[1555,1589],[1556,1590],[1557,1591],[1558,1592],[1559,1593],[1564,1598],[1565,1599],[1566,1600],[1567,1601],[1568,1602]]
imod = [1558,1565]
hx = listofhisto('hx',[],200,-0.5,0.5,title=';x_{hycal};')
hxc = listofhisto('hxc',[],200,-0.5,0.5,title=';x_{hycal};')
d = cell_size[0][0]
for t in ts:
  for x in progress(t,show=3,n=t.GetEntries(),precision=1):
    for k in range(x.n_cl):
      if x.nh[k]<2: continue
      if x.E[k]<1000: continue
      for i in range(2):
        if x.id[k]!=imod[i]: continue
        t = (x.xhycal[k]-0.605-module_pos[x.id[k]][0])/cell_size[0][0]
        _ = hx.Fill(t)
        _ = hxc.Fill(correct(t))
        

        
  
for i in range(5):
  _ = hx[10+i].Add(hx[i])
  _ = hx[10+i].Add(hx[9-i])

h2x = listofhisto('h2x',range(5),100,0,0.5,title=';|x_{hycal}-x_{center}|;')
for i in range(5):
  for j in range(100):
    h2x[i][j+1] = hx[i+10][j+101]+hx[i+10][101-j]

h3x = listofhisto('h3x',range(5),100,0,0.5,title=';|x_{hycal}-x_{center}|/cell_size;')
for i in range(5):
  h3x[i].Add(h2x[i],100./sum(h2x[i][j+1] for j in range(100)))


# fit of module density distribution
f = TFile(workf+'/density_all2.root')
ffit = TF1('ffit','2*[0]*x^2*(-[3] + x^2)*([1]*x^2 + [2] + x^4) + 2*[0]*x^2*(x^2 - 0.25)*([1]*x^2 + [2] + x^4) + [0]*x*(-[3] + x^2)*(x^2 - 0.25)*(2*[1]*x + 4*x^3) + [0]*(-[3] + x^2)*(x^2 - 0.25)*([1]*x^2 + [2] + x^4)+[4]',0.,0.5)
h = f.final.Get('hx2_101_0')
h.Fit('ffit','lr','',0.015,0.5)
params = [ffit.GetParameter(i) for i in range(4)]




# looking at inelastic distributions
lpi0 = readfile('/work/hallb/prad/mlevilla/new_evgen_files/inelastic_nick/2GeV/inel_2GeV_1e4_pi0_0.dat',[float]*16)
lpip = readfile('/work/hallb/prad/mlevilla/new_evgen_files/inelastic_nick/2GeV/inel_2GeV_1e4_pip_0.dat',[float]*16)
lmoller = readfile('/work/hallb/prad/mlevilla/new_evgen_files/moller/2GeV/moller_2GeV_1e6_th04_35_CO1_2000_0.dat',[float]*9)
lep = readfile('/work/hallb/prad/mlevilla/new_evgen_files/ep_FF3/2GeV/ep_2GeV_4e5_th05_75_CO1_2000_FF3_0_e-.dat',[float]*9)

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta, ROOT.kBlack]
hth = listofhisto('hth',['ep','moller','pi0','pip','sum'],150,0.,7.5,title=';#theta (deg);')
he = listofhisto('he',['ep','moller','pi0','pip','sum'],250,0.,2500.,title=';E (MeV);')
he2 = listofhisto('he2',['ep','moller','pi0','pip','sum'],250,0.,2500.,title=';E (MeV);')
heth = listofhisto('heth',['ep','moller','pi0','pip','sum'],150,0.,7.5,250,0.,2500.,title=';E (MeV);')

for i,l in enumerate([lep,lmoller]):
  for j in range(len(l)):
    for k in [[0,2],range(3)][i]:
      if l[j][3*k]!=0.: 
        th = l[j][3*k+1]/degrad
        _ = he[i].Fill(l[j][3*k])
        _ = hth[i].Fill(l[j][3*k+1]/degrad)
        _ = heth[i].Fill(l[j][3*k+1]/degrad,l[j][3*k])
        if th>3:  _ = he2[i].Fill(l[j][3*k])

masses = [[0.51099893,134.9,939.56563,0.],[0.51099893,139.57,939,0]]
for i,l in enumerate([lpi0,lpip]):
  for j in range(len(l)):
    for k in range(4):
      if l[j][4*k+3]>0. and l[j][4*k]<=22:
        E = (sum(l[j][4*k+k2]**2 for k2 in range(1,4)))**0.5
        th = atan2((l[j][4*k+1]**2+l[j][4*k+2]**2)**0.5,l[j][4*k+3])/degrad
        _ = he[i+2].Fill(E)
        _ = hth[i+2].Fill(th)
        _ = heth[i+2].Fill(th,E)
        if th>3:  _ = he2[i+2].Fill(th)

luminosity = [449.,2*108.69,1855.,4200.]
for i in range(4):
  he[i].Scale(1./luminosity[i])
  he2[i].Scale(1./luminosity[i])
  hth[i].Scale(1./luminosity[i])
  heth[i].Scale(1./luminosity[i])
  he[4].Add(he[i])
  he2[4].Add(he2[i])
  hth[4].Add(hth[i])
  heth[4].Add(heth[i])

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kBlack]
for i in range(5):
  he[i].SetLineColor(colors[i])
  he2[i].SetLineColor(colors[i])
  hth[i].SetLineColor(colors[i])
  heth[i].SetLineColor(colors[i])


# looking at simyield/born
from esepp import *
thbin = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]

h = listofhisto('hth',['yields','ratio'],len(thbin)-1,xbin=thbin,title=';#theta (deg);')

for i in progress(300,show=3,precision=1):
  l = readfile(workf+'/new_evgen_files/ep_FF5/2GeV/ep_2GeV_1e6_th05_75_CO1_2000_FF5_'+str(i)+'_e-.dat',[float]*9)
  for x in l:
    E = x[0]
    th = x[1]
    Etheo = ep_energy_el(th,2.142)
    if abs(E-Etheo*1000)>140: continue
    _ = h[0].Fill(th/degrad)

      
thbar = [0.707, 0.723, 0.739, 0.757, 0.775, 0.795, 0.817, 0.84, 0.865, 0.893, 0.924, 0.957, 0.994, 1.035, 1.081, 1.13, 1.183, 1.24, 1.303, 1.376, 1.464, 1.572, 1.707, 1.888, 2.101, 2.344, 2.633, 2.933, 3.236, 3.536, 3.836, 4.137, 4.437, 4.737, 5.041]
thbar2 = [0.706, 0.722, 0.738, 0.755, 0.774, 0.794, 0.815, 0.838, 0.864, 0.892, 0.922, 0.956, 0.993, 1.033, 1.079, 1.129, 1.182, 1.238, 1.301, 1.375, 1.463, 1.571, 1.707, 1.89, 2.102, 2.348, 2.637, 2.936, 3.236, 3.535, 3.834, 4.134, 4.433, 4.733, 5.036]
sigma = [born(x*degrad,2.142,flag=3) for x in thbar2]
luminosity = 1123.2e6*300
yields_theo = [sigma[i]*2*pi*(cos(thbin[i]*degrad)-cos(thbin[i+1]*degrad))*luminosity for i in range(len(thbin)-1)]
yields_theo = [sigma[i]*2*pi*sin(thbar2[i]*degrad)*(thbin[i+1]-thbin[i])*degrad*luminosity*200 for i in range(len(thbin)-1)]

for i in range(35):
  h[1][i+1] = h[0][i+1]/yields_theo[i]
  h[1].SetBinError(i+1,h[0][i+1]**0.5/yields_theo[i])


# changing multiple names
indir = workf+'/simulation_files/inelastic/1GeV/'
name = 'pb_inel_gen_trgeff_'
new_name = 'inel_1GeV_1e4_th03_75_trgeff_'
ext = '_rec.root'
nf = 300
for i in range(nf): os.rename(indir+name+str(i)+ext,indir+new_name+str(i)+ext)

# comparing peter and nick's inelastic generator

c1,c2,c3,c4 = TChain('T'),TChain('T'),TChain('T'),TChain('T')
for i in range(100):
  _ = c1.Add(workf+'/simulation_files/ep_FF3/2GeV/ep_2GeV_4e5_th05_75_CO1_2000_FF3_trgeff_'+str(i)+'_rec.root')
  _ = c3.Add(workf+'/simulation_files/inelastic_nick/2GeV/inel_2GeV_1e4_pi0_trgeff_'+str(i)+'_rec.root')
  _ = c4.Add(workf+'/simulation_files/inelastic_nick/2GeV/inel_2GeV_1e4_pip_trgeff_'+str(i)+'_rec.root')
  _ = c2.Add(workf+'/simulation_files/inelastic/2GeV/inel_2GeV_1e4_th03_75_trgeff_'+str(i)+'_rec.root')

fout = TFile('temp.root','recreate')
c1.Draw('HC.P>>hep(500,0,2500)','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>1 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<2')
c2.Draw('HC.P>>hinel(500,0,2500)','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>1 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<2')         
c3.Draw('HC.P>>hpi0(500,0,2500)','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>1 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<2')
c4.Draw('HC.P>>hpip(500,0,2500)','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>1 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<2')

names = ['ep','inel','pi0','pip']
h = [fout.Get('h'+x) for x in names]
scales = [cs_dict['ep_2GeV']/4e7,cs_dict['inelastic_2GeV']/1e6,cs_dict['inelastic_2GeV_pi0']/1e6,cs_dict['inelastic_2GeV_pip']/1e6]
for i in range(4):
  _ =h[i].Scale(scales[i])

for i in range(86,181):
  h[2][i] = h[2][85]+(h[2][181]-h[2][85])*float(i-81)/float(181-85)

h.extend([TH1F('hinel2','',500,0,2500),TH1F('hsum1','',500,0,2500),TH1F('hsum2','',500,0,2500)])
for i,j in [(4,[2,3]),(5,[0,1]),(6,[0,4])]:
  for k in j: _ = h[i].Add(h[k])

colors = [ROOT.kBlue,ROOT.kGreen+2,ROOT.kCyan+1,ROOT.kMagenta+2,ROOT.kYellow-3,ROOT.kRed,ROOT.kOrange-3]
for i in range(7): 
  h[i].GetXaxis().SetTitle('E (MeV)')
  h[i].GetYaxis().SetTitle('A.U.')
  h[i].SetTitle('')
  h[i].SetLineColor(colors[i])
  h[i].SetOption('hist')

# comparing peter and nick's inelastic generator and data
runs_prod = get_runs_between(1288,1341,'prod')
runs_bg = get_runs_between(1288,1341,'empty')

h = listofhisto('hE',['prod','empty','ep','moller','boste','pi0','pip','nick','sum_boste','sum_nick'],2500,0,2500,title=';E (MeV);count/MeV')

cprod,cbg,cep,cmoller,cboste,cpi0,cpip = TChain('event'),TChain('event'),TChain('T'),TChain('T'),TChain('T'),TChain('T'),TChain('T')

for i in runs_prod: _=cprod.Add(volf+'/trees/trees_sel6/tree_island_'+str(i)+'.root')

for i in runs_bg: _=cbg.Add(volf+'/trees/trees_sel6/tree_island_'+str(i)+'.root')

cprod.Draw('E>>hE_prod','atan2(sqrt(xgem^2+ygem^2),zgem)*180/pi>3.3 && atan2(sqrt(xgem^2+ygem^2),zgem)*180/pi<3.4')
cbg.Draw('E>>hE_empty','atan2(sqrt(xgem^2+ygem^2),zgem)*180/pi>3.3 && atan2(sqrt(xgem^2+ygem^2),zgem)*180/pi<3.4')

for i in range(300):
  _=cep.Add(volf+'/simulation_files_rec/ep_FF3/1GeV/ep_1GeV_4e5_th05_75_CO1_1000_FF3_'+str(i)+'_rec.root')
  _=cmoller.Add(volf+'/simulation_files_rec/moller/1GeV/moller_1GeV_1e6_th04_75_CO1_1000_'+str(i)+'_rec.root')
  _=cboste.Add(volf+'/simulation_files_rec/inelastic_bosted/1GeV/inel_1GeV_1e4_th03_75_'+str(i)+'_rec.root')
  _=cpi0.Add(volf+'/simulation_files_rec/inelastic_markov/1GeV/inel_1GeV_1e4_pi0_'+str(i)+'_rec.root')
  _=cpip.Add(volf+'/simulation_files_rec/inelastic_markov/1GeV/inel_1GeV_1e4_pip_'+str(i)+'_rec.root')

cep.Draw('HC.P>>hE_ep','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>3.3 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<3.4')
cmoller.Draw('HC.P>>hE_moller','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>3.3 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<3.4')
cboste.Draw('HC.P>>hE_boste','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>3.3 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<3.4')
cpi0.Draw('HC.P>>hE_pi0','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>3.3 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<3.4')
cpip.Draw('HC.P>>hE_pip','GEM.Z>0 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi>3.3 && atan2(sqrt(GEM.X^2+GEM.Y^2),GEM.Z-89)*180/pi<3.4')

nevents = [300*x for x in [4e5,1e6,1e4,1e4,1e4]]
cs = [cs_dict['ep_1GeV'],cs_dict['moller_1GeV'],cs_dict['inelastic_1GeV'],cs_dict['inelastic_1GeV_pi0'],cs_dict['inelastic_1GeV_pip']]
lum_simu = [x/y for x,y in zip(nevents,cs)]
load_live_charge()
load_thickness() 
lcharge_prod = sum(live_charge[run]*1e4 for run in runs_prod) 
lcharge_empty = sum(live_charge[run]*1e4 for run in runs_bg) 
average_thick = sum(live_charge[run]*1e4*thickness[run] for run in runs_prod)/lcharge_prod
lum_prod = lcharge_prod/1e28*average_thick*6.242e9
gem_efficiency = 0.914
h[0].Add(h[1],-lcharge_empty/lcharge_prod)
h[0].Scale(1/gem_efficiency)
for i in range(5):
  h[i+2].Scale(lum_prod/lum_simu[i])

for i,j in [(7,[5,6]),(8,[2,3,4]),(9,[2,3,7])]:
  for k in j:
    _=h[i].Add(h[k])

colors = [ROOT.kBlue,ROOT.kRed,ROOT.kGreen+2,ROOT.kCyan,ROOT.kMagenta,ROOT.kYellow+2,ROOT.kRed-6,ROOT.kOrange,ROOT.kBlack]
for i in range(8): 
  h[i+2].SetLineColor(colors[i])
  h[i+2].SetMarkerColor(colors[i])


# binning extraction
runs1,runs2 = get_runs_between(1070,1341,'prod'),get_runs_between(1362,1516,'prod')
t1,t2 = TChain('clusters'),TChain('clusters')
ftemp = TFile('temp.root','recreate')
for r in runs1: _=t1.Add(workf+'/yields/yields_'+str(r)+'.root')

for r in runs2: _=t2.Add(workf+'/yields/yields_'+str(r)+'.root')

t1.Draw('theta_gem>>h1gev(10000,0,10)','event==0')
t2.Draw('theta_gem>>h2gev(10000,0,10)','event==0')

h1,h2 = ftemp.Get('h1gev'),ftemp.Get('h2gev')
xbin1,xbin2 = [0.85],[0.7]
aux1,aux2 = 0.,0.
test = []
i = 700 
while i<10000:
  while aux2==0 or aux2<1000000:
    aux2+=h2[i+1]
    i+=1
  i+=1
  xbin2.append(h2.GetBinLowEdge(i))
  test.append(aux2)
  aux2 = 0
  
# checking energy resolution fluctuation
load_live_charge()

runs_prod = get_runs_between(1443,1490,'prod')
runs_empty = get_runs_between(1443,1490,'empty')

lcharge_prod = sum(live_charge[r]*1e4 for r in runs_prod)
lcharge_empty = sum(live_charge[r]*1e4 for r in runs_empty)

hth_prod = listofhisto('hth_prod',range(150),1000,-1,1,title=';E\'/E_{th}-1;')
htx_prod = listofhisto('htx_prod',[],1000,-1,1,title=';;')
htxe_prod = listofhisto('htxe_prod',[],1000,-1,1,title=';;E\'/E_{th}-1')
hth_empty = listofhisto('hth_empty',range(150),1000,-1,1,title=';E\'/E_{th}-1;')
htx_empty = listofhisto('htx_empty',[],1000,-1,1,title=';;')
htxe_empty = listofhisto('htxe_empty',[],1000,-1,1,title=';;E\'/E_{th}-1')
hth_sub = listofhisto('hth_sub',range(150),1000,-1,1,title=';E\'/E_{th}-1;')
htx_sub = listofhisto('htx_sub',[],1000,-1,1,title=';;')
htxe_sub = listofhisto('htxe_sub',[],1000,-1,1,title=';;E\'/E_{th}-1')

for r in runs_prod:
  f = TFile(workf+'/eresolution_'+str(r)+'_5.root')
  _ = htx_prod.Add(f.Get('htx'))
  _ = htxe_prod.Add(f.Get('htxe'))
  for i in range(150):
    _ = hth_prod[i].Add(f.Get('hth_'+str(i)))
  
for r in runs_empty:
  f = TFile(workf+'/eresolution_'+str(r)+'_5.root')
  _ = htx_empty.Add(f.Get('htx'))
  _ = htxe_empty.Add(f.Get('htxe'))
  for i in range(150):
    _ = hth_empty[i].Add(f.Get('hth_'+str(i)))

htx_sub.Add(htx_prod,htx_empty,1,-lcharge_prod/lcharge_empty)
htxe_sub.Add(htxe_prod,htxe_empty,1,-lcharge_prod/lcharge_empty)
htxe_sub.Divide(htx_sub)
htxe_sub.Scale(2.142**0.5)
htxe_prod.Divide(htx_prod)
htxe_prod.Scale(2.142**0.5)

for i in range(150):
  _ = hth_sub[i].Add(hth_prod[i],hth_empty[i],1,-lcharge_prod/lcharge_empty)

ffit = TF1('2gaus','[0]*exp(-(x-[1])^2/[2]^2/2)+ (x<[1] ? [3]*1/(1+[4]*(x-[1])^2) : [5]*exp(-(x-[1])^2/[6]^2/2))',-0.3,0.3)
ffit.SetParLimits(0,1e0,1e7)
ffit.SetParLimits(1,-0.01,0.01)
ffit.SetParLimits(2,0.01,0.1)
ffit.SetParLimits(3,1e-1,1e7)
ffit.SetParLimits(5,1e-1,1e7)
ffit.SetParLimits(6,0.02,1)

ffit = TF1('cball','[0]* ((x-[1])/[2]>-[3] ? exp(-(x-[1])^2/[2]^2/2) : ([4]/abs([3]))^[4] * exp(-[3]^2/2) * (([4]/abs([3])-abs([3]))-(x-[1])/[2])^[4])',-0.3,0.3)

for i in range(8,150):
  if i<20: hth_sub[i].GetXaxis().SetRangeUser(-0.2,0.2)
  elif i<40: hth_sub[i].GetXaxis().SetRangeUser(-0.3,0.3)
  else: hth_sub[i].GetXaxis().SetRangeUser(-0.5,0.3) 
  hth_sub[i].Fit(ffit,'l')
  hth_sub[i].Fit(ffit,'l')

for i in range(4,60):
  hth_sub[i].Fit('gaus','r','',-0.025,0.025)

for i in range(30):
  hth_sub[i+60].Fit('gaus','r','',-0.05,0.05)

sigma,err = [],[]
for i in range(90):
  if hth_sub[i].GetFunction('gaus') == None: 
    sigma.append(0)
    err.append(0)
  else: 
    sigma.append(hth_sub[i].GetFunction('gaus').GetParameter(2)*2.142**0.5)
    err.append(hth_sub[i].GetFunction('gaus').GetParError(2)*2.142**0.5)

thbin = [0.525+0.05*k for k in range(60)]+[3.55+0.1*k for k in range(30)]
g = tgraph(thbin,sigma,dy=err,title=';#theta;#sigma_{E}/#sqrt{E}')

params = [[0 for i in range(7)] for j in range(8)]
err = [[0 for i in range(7)] for j in range(8)]
for i in range(8,150):
  params.append([hth_sub[i].GetFunction('2gaus').GetParameter(j) for j in range(7)])
  err.append([hth_sub[i].GetFunction('2gaus').GetParError(j) for j in range(7)])

thbin = [0.5125+0.025*k for k in range(120)]+[3.55+0.1*k for k in range(30)]
g = [tgraph(thbin,[params[i][j] for i in range(150)],[err[i][j] for i in range(150)], title = ';#theta (deg);') for j in range(7)]


# energy correction vs coordinate
f = TFile('/work/hallb/prad/mlevilla/density/density_ene_42_2GeV_correct.root')
ffit = TF1('ffit','[0]*(1+[1]*exp(-[2]/x^2))',-0.5,0.5)
ffit.SetParameters(2100,0.04,0.3)
ffit.SetParLimits(2,0.1,0.5)

lparam = [[] for i in range(112)]
for i in range(112):
    h = f.raw_correct.Get('hxc_'+str(i)+'_3')
    hE = f.raw_correct.Get('hxcE_'+str(i)+'_3')
    hE.Divide(h)
    hE.Fit(ffit,'','',-0.5,0.5)
    lparam[i] = [ffit.GetParameter(k) for k in range(3)]

# looking at chi2 before abd after energy_correction
f = TFile('ecorrect2D_fit_test_correct_2GeV.root')
region_names = ['pwo','lg_top','lg_right','lg_bottom','lg_left','lg']
h = [f.regions.Get('hE_{0}_ep'.format(x)) for x in region_names]
hc = [f.regions_correct.Get('hEc_{0}_ep'.format(x)) for x in region_names]
chi = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in h]
chic = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in hc]

ffit = h.GetFunction('ffit')
chi2 = (sum(sum((h.GetBinContent(6+i,6+j)-ffit.Eval(-0.45+0.1*i,-0.45+0.1*j))**2 for i in range(10)) for j in range(10))/100)**0.5

ffit2 = TF2('ffit2','[0]*(1+[1]*x^2+[2]*y^2+[3]*x^2*y^2+[4]*x^4+[5]*y^4+[6]*x^2*y^4+[7]*x^4*y^2+[8]*x^4*y^4)',-0.5,0.5,-0.5,0.5)
ffit2 = TF2('ffit2','[0]*(1+[1]*x^2+[2]*y^2+[3]*x^2*y^2+[4]*x^4+[5]*y^4)',-0.5,0.5,-0.5,0.5)
chi3 = (sum(sum((h.GetBinContent(6+i,6+j)-ffit2.Eval(-0.45+0.1*i,-0.45+0.1*j))**2 for i in range(10)) for j in range(10))/100)**0.5

h_lg = [f.Edist.Get('hxyE_'+str(i)+'_0') for i in range(64)]
h_pwo_ep = [f.Edist.Get('hxyE_'+str(i)+'_0') for i in range(64,353)]
h_pwo_ee = [f.Edist.Get('hxyE_'+str(i)+'_1') for i in range(64,353)]

chi2_lg = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in h_lg]
chi2_pwo_ep = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in h_pwo_ep]
chi2_pwo_ee = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in h_pwo_ee]

hc_lg = [f.Edist_correct.Get('hxyEc_'+str(i)+'_0') for i in range(64)]
hc_pwo_ep = [f.Edist_correct.Get('hxyEc_'+str(i)+'_0') for i in range(64,353)]
hc_pwo_ee = [f.Edist_correct.Get('hxyEc_'+str(i)+'_1') for i in range(64,353)]

chi2c_lg = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in hc_lg]
chi2c_pwo_ep = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in hc_pwo_ep]
chi2c_pwo_ee = [(sum(sum((x.GetBinContent(6+i,6+j)-1)**2 for i in range(10)) for j in range(10))/100)**0.5 for x in hc_pwo_ee]

# looking at chi2 before abd after energy_correction for groups
f = TFile(workf+'/energy_correction/ecorrect2D_fit_36_group_correct_2GeV.root')
region_names = ['pwo','lg_top','lg_right','lg_bottom','lg_left','lg']
h = [f.regions.Get('hE_{0}_ep'.format(x)) for x in region_names]
hc = [f.regions_correct.Get('hEc_{0}_ep'.format(x)) for x in region_names]
chi = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in h]
chic = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in hc]

h_lg = [f.Edist.Get('hxyE_ep_'+str(i)) for i in range(42)]
h_pwo_ep = [f.Edist.Get('hxyE_ep_'+str(i)) for i in range(42,392)]
h_pwo_ee = [f.Edist.Get('hxyE_ee_'+str(i)) for i in range(447)]

chi2_lg = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in h_lg]
chi2_pwo_ep = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in h_pwo_ep]
chi2_pwo_ee = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in h_pwo_ee]

hchi2_lg = TH1F('hchi2_lg',';#delta;',100,0.,0.1)
hchi2_pwo_ep = TH1F('hchi2_pwo_ep',';#delta;',100,0.,0.03)
hchi2_pwo_ee = TH1F('hchi2_pwo_ee',';#delta;',100,0.,0.03)

hc_lg = [f.Edist_correct.Get('hxyEc_ep_'+str(i)) for i in range(42)]
hc_pwo_ep = [f.Edist_correct.Get('hxyEc_ep_'+str(i)) for i in range(42,392)]
hc_pwo_ee = [f.Edist_correct.Get('hxyEc_ee_'+str(i)) for i in range(447)]

chi2c_lg = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in hc_lg]
chi2c_pwo_ep = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in hc_pwo_ep]
chi2c_pwo_ee = [(sum(sum((x.GetBinContent(1+i,1+j)-1)**2 for i in range(5)) for j in range(5))/25)**0.5 for x in hc_pwo_ee]

hchi2c_lg = TH1F('hchi2c_lg',';#delta;',100,0.,0.1)
hchi2c_pwo_ep = TH1F('hchi2c_pwo_ep',';#delta;',100,0.,0.05)
hchi2c_pwo_ee = TH1F('hchi2c_pwo_ee',';#delta;',100,0.,0.05)

for x,y in zip([chi2_lg,chi2_pwo_ep,chi2_pwo_ee,chi2c_lg,chi2c_pwo_ep,chi2c_pwo_ee],[hchi2_lg,hchi2_pwo_ep,hchi2_pwo_ee,hchi2c_lg,hchi2c_pwo_ep,hchi2c_pwo_ee]):
  for z in x: _=y.Fill(z)

l = readfile(workf+'/energy_correction/ecorrect2D_36_group_2GeV.txt',[float]*5,cols=range(1,6))
offset_lg = [x[0] for x in l[:42]]
offset_pwo_ep = [x[0] for x in l[42:392]]
offset_pwo_ee = [x[0] for x in l[392:]]
amplitude_lg = [0.25*(x[1]+x[2])+0.0625*(x[3]+x[4]+x[5]) for x in l[:42]]
amplitude_pwo_ep = [0.25*(x[1]+x[2])+0.0625*(x[3]+x[4]+x[5]) for x in l[42:392]]
amplitude_pwo_ee = [0.25*(x[1]+x[2])+0.0625*(x[3]+x[4]+x[5]) for x in l[392:]]
hoffset_lg = TH1F('hoffset_lg',';p_{0};',100,0.9,1)
hoffset_pwo_ep = TH1F('hoffset_pwo_ep',';p_{0};',100,0.95,1.05)
hoffset_pwo_ee = TH1F('hoffset_pwo_ee',';p_{0};',100,0.95,1.05)
hamplitude_lg = TH1F('hamplitude_lg',';1/4*(p_{1}+p_{2})+1/16*(p_{3}+p_{4}_p_{5});',100,-0.1,0.3)
hamplitude_pwo_ep = TH1F('hamplitude_pwo_ep',';1/4*(p_{1}+p_{2})+1/16*(p_{3}+p_{4}_p_{5});',100,-0.01,0.1)
hamplitude_pwo_ee = TH1F('hamplitude_pwo_ee',';1/4*(p_{1}+p_{2})+1/16*(p_{3}+p_{4}_p_{5});',100,-0.01,0.1)

for x,y in zip([offset_lg,offset_pwo_ep,offset_pwo_ee,amplitude_lg,amplitude_pwo_ep,amplitude_pwo_ee],[hoffset_lg,hoffset_pwo_ep,hoffset_pwo_ee,hamplitude_lg,hamplitude_pwo_ep,hamplitude_pwo_ee]):
  for z in x:
    _=y.Fill(z)

# look at elastic peak shift before/after s-shape correction
from viewer import *
lep = readfile(conff+'/groupindex_1GeV_ep.txt',[int]*2,out=dict)
lee = readfile(conff+'/groupindex_1GeV_ee.txt',[int]*2,out=dict)
f = TFile(workf+'/energy_correction/sshape_36_10_1GeV.root')
h = getlisto(f.elasticity)
hep = [x for x in h if 'ep' in x.GetName()]
hee = [x for x in h if 'ee' in x.GetName()]
shift_ep = [fit_gaus(x,sigma=1)[0] for x in hep]
shift_ee = [fit_gaus(x,sigma=1)[0] for x in hee]
shift_ep2 = [shift_ep[x] if x!=-1 else -1 for x in lep.values()] 
shift_ee2 = [shift_ee[x] if x!=-1 else -1 for x in lee.values()] 

v = hcviewer('shift',zmin=-0.03,zmax=0.03)
v.fill(shift_ee2,reset=True)
v.fill([-1,-1,-1,-1],[1561,1562,1595,1596])
v.draw()
