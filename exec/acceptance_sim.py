#!/usr/bin/env python
from params import *

largs = [('n','','','suffix'),
         ('b',False,'','batch_mode'),
         ('m','island','','clustering_method'),
         ('s',0,'','show_bar'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('phi',[10.,5.],'','deltaphi_cut'),
         ('z',[500.,150.],'','zvertex_cut'),
         ('sim','','','simulation_file_name'),
         ('ebeam',2,'','beam_energy_for_simulation'),
         ('nf',1,'','number_of_simulation_files'),
         ('f',False,'','final merging')]
print_help(largs)
[suffix,batch,method,show,outdir,limit,phicut,zcut,sim,ebeam,nf,final] = catch_args(largs)

if batch:
  jsub_fast(largs,'acceptance',['b','nf','o','s'],[['n','sim'],['_'+str(x) for x in range(nf)]],outdir)
  sys.exit()

if final:
  from misc_root import *
  f = [TFile(outdir+'/acceptance{0}_{1}.root'.format(suffix,i)) for i in range(nf)]
  lhth0 = [getlisto(x.hth0) for x in f]
  lhthall = [getlisto(x.hthall) for x in f]
  lhth1 = [getlisto(x.hth1) for x in f]
  lhth2 = [getlisto(x.hth2) for x in f]
  laccall = [getlisto(x.accall) for x in f]
  lacc1 = [getlisto(x.acc1) for x in f]
  lacc2 = [getlisto(x.acc2) for x in f]
  for l in [lhth0,lhthall,lhth1,lhth2]:
    for x in l[1:]:
      for i in range(len(l[0])):
        l[0][i].Add(x[i])
  hth0, hthall, hth, haccall, hacc =  lhth0[0], lhthall[0], [lhth1[0],lhth2[0]], lhaccall[0], [lhacc1[0],lhacc2[0]]
  for n in range(2):
    for k in range(len(cnames)):
      haccall[n][k].Divide(hthall[n][k],hth0[n])
      for m in range([nbin1,nbin2][n]): 
        temp = max(0,2*hth0[n][m+1]-hthall[n][k][m+1]) if hth0[n][m+1]<hthall[n][k][m+1] else hthall[n][k][m+1]
        haccall[n][k].SetBinError(m+1,((temp+1.)*(temp+2.)/(hth0[n][m+1]+2.)/(hth0[n][m+1]+3.)-(temp+1.)**2/(hth0[n][m+1]+2.)**2)**0.5)
      for i in range(len(enames)):
        for j in range(len(elasnames)):
          hacc[n][i][j][k].Divide(hth[n][i][j][k],hth0[n])
          for m in range([nbin1,nbin2][n]): 
            temp = max(0,2*hth0[n][m+1]-hth[n][i][j][k][m+1]) if hth0[n][m+1]<hth[n][i][j][k][m+1] else hth[n][i][j][k][m+1]
            hacc[n][i][j][k].SetBinError(m+1,((temp+1.)*(temp+2.)/(hth0[n][m+1]+2.)/(hth0[n][m+1]+3.)-(temp+1.)**2/(hth0[n][m+1]+2.)**2)**0.5)

  # writing
  fout = TFile(outdir+'/acceptance{0}_all.root'.format(suffix),'recreate')
  writelisto(hth0,fout,'hth0')
  writelisto(hthall,fout,'hthall')
  writelisto(hth[0],fout,'hth1')
  writelisto(hth[1],fout,'hth2')
  writelisto(haccall,fout,'accall')
  writelisto(hacc[0],fout,'acc1')
  writelisto(hacc[1],fout,'acc2')
  sys.exit()

print ' '.join(sys.argv)
date = subprocess.check_output(['date'])
print 'start: ',date[:-1]

from misc_root import *
from selection import *

f = TFile(sim+'_rec.root')
t = f.T
t.AddFriend('T',sim+'.root')
t.SetBranchStatus("*",0)
for x in ['GEM.X','GEM.Y','GEM.Z','HC.N','HC.X','HC.Y','HC.Z','HC.P','HC.CID','GUN.Theta']: t.SetBranchStatus(x,1)
fout = TFile(outdir+'/acceptance{0}.root'.format(suffix),'recreate')

run = 1443 if ebeam==2 else 1288
load_run_params(run)
Ebeam = beam_energy[0]/1000.
max_energy = 1500. if Ebeam<2. else 2500.

thbin1 = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
thbin2 = [0.5+0.02*i for i in range(100)]+[2.5+0.1*i for i in range(21)]
nbin1 = len(thbin1)-1
nbin2 = len(thbin2)-1

enames = ['ep','ee1','ee2','ee3']
elasnames = range(5,0,-1)
cnames = ['elas','wfid','wdead','spacer']

hth0_1 = listofhisto('hth0_1',[],nbin1,xbin=thbin1,title=';#theta (deg);')
hth0_2 = listofhisto('hth0_2',[],nbin2,xbin=thbin2,title=';#theta (deg);')
hth0 = [hth0_1,hth0_2]

hthall = listofhisto('hthall',[range(1,3),cnames],[nbin1,nbin2],xbin=[thbin1,thbin2],title=';#theta (deg);')
haccall = listofhisto('haccall',[range(1,3),cnames],[nbin1,nbin2],xbin=[thbin1,thbin2],title=';#theta (deg);')
hth = listofhisto('hth',[range(1,3),enames,elasnames,cnames],[nbin1,nbin2],xbin=[thbin1,thbin2],title=';#theta (deg);')
hacc = listofhisto('hacc',[range(1,3),enames,elasnames,cnames],[nbin1,nbin2],xbin=[thbin1,thbin2],title=';#theta (deg);')

for _ in progress(t,show=show,n=t.GetEntries(),precision=2,limit=limit):

  if getattr(t,'HC.N')!=1 or getattr(t,'GEM.N')!=1: continue
  for h in hth0: h.Fill(getattr(t,'GUN.Theta')[0]/degrad)

  [E,c,idx,lg,cid,_] = get_variables_sim(t,match=1,do_matching=1)
  theta = [[ftheta(y) if y!=[] else 0 for y in x] for x in c]
  thetadeg = [[abs(degrees(y)) for y in x] for x in theta]
  
  # all
  for j in range(len(idx)):
    th = thetadeg[j][1]
    for h in hthall: h[0].Fill(th)
    if fiducial_cut(cid[j],c[j][1],1,0,0.): continue
    for h in hthall: h[1].Fill(th)
    if fiducial_cut(cid[j],c[j][1],1,1,1.): continue
    for h in hthall: h[2].Fill(th)
    if gem_spacers(c[j][1][0],c[j][1][1]) or gem_dead(c[j][1][0],c[j][1][1]): continue
    for h in hthall: h[3].Fill(th)

  # ep 
  for j in range(len(idx)):
    th = thetadeg[j][1]
    E_ep = ep_energy_el(th*degrad,Ebeam)
    elas_ep = E[j]/E_ep/1000-1
    elasc = [0.024,0.062][lg[j]]/E_ep**0.5

    if abs(elas_ep)>6*elasc: continue

    for i in range(5):
      if abs(elas_ep)>(5-i)*elasc: continue
      for h in hth: h[0][i][0].Fill(th)
      if fiducial_cut(cid[j],c[j][1],1,0,0.): continue
      for h in hth: h[0][i][1].Fill(th)
      if fiducial_cut(cid[j],c[j][1],1,1,1.): continue
      for h in hth: h[0][i][2].Fill(th)
      if gem_spacers(c[j][1][0],c[j][1][1]) or gem_dead(c[j][1][0],c[j][1][1]): continue
      for h in hth: h[0][i][3].Fill(th)

 
  # ee1 
  for j in range(len(idx)):
    th = thetadeg[j][1]
    E_m = moller_energy(th*degrad,Ebeam)
    elas_m = E[j]/E_m/1000-1
    elasc = [0.024,0.062][lg[j]]/E_m**0.5

    for i in range(5):
      if abs(elas_m)>(5-i)*elasc: continue
      for h in hth: h[1][i][0].Fill(th)
      if fiducial_cut(cid[j],c[j][1],1,0,0.): continue
      for h in hth: h[1][i][1].Fill(th)
      if fiducial_cut(cid[j],c[j][1],0,1,1.): continue
      for h in hth: h[1][i][2].Fill(th)
      if gem_spacers(c[j][1][0],c[j][1][1]) or gem_dead(c[j][1][0],c[j][1][1]): continue
      for h in hth: h[1][i][3].Fill(th)

  # ee2/ee3
  for k in range(2):
    dEmin, jmin =10000, [-1,-1]
    for j1 in range(len(idx)-1):
      for j2 in range(j1+1,len(idx)):
        dphi = degrees((fphi(c[j1][1-k])-fphi(c[j2][1-k]))%(2*pi)-pi)
        cproj = [[c[j][1-k][m]*hycal_center[2]/c[j][1-k][2] for m in range(2)] for j in [j1,j2]]
        rj2 = [(x[0]**2+x[1]**2)**0.5 for x in cproj]
        zvertex = ((m_e+Ebeam)*rj2[0]*rj2[1]/2./m_e)**0.5-hycal_center[2]
        dE = E[j1]+E[j2]-Ebeam*1e3
        dEcut = (([0.024,0.062][lg[j1]]*E[j1]**0.5)**2+([0.024,0.062][lg[j2]]*E[j2]**0.5)**2)**0.5*1000**0.5
        if abs(dphi)>phicut[1-k] or abs(zvertex)>zcut[1-k] or abs(dE)>6*dEcut: continue

        if dE<dEmin:
          dEmin = dE
          jmin = [j1,j2]
    
    if dEmin!=10000:
      for j in jmin:
        th = thetadeg[j][1]

        E_m = moller_energy(th*degrad,Ebeam)
        elas_m = E[j]/E_m/1000-1
        elasc = [0.024,0.062][lg[j]]/E_m**0.5
        for i in range(5):
          if abs(elas_m)>(5-i)*elasc: continue
          for h in hth: h[2+k][i][0].Fill(th)
          if fiducial_cut(cid[j],c[j][1],1,0,0.): continue
          for h in hth: h[2+k][i][1].Fill(th)
          if fiducial_cut(cid[j],c[j][1],0,1,1.): continue
          for h in hth: h[2+k][i][2].Fill(th)
          if gem_spacers(c[j][1][0],c[j][1][1]) or gem_dead(c[j][1][0],c[j][1][1]): continue
          for h in hth: h[2+k][i][3].Fill(th)

# acceptance calculation
for n in range(2):
  for k in range(len(cnames)):
    haccall[n][k].Divide(hthall[n][k],hth0[n])
    for m in range([nbin1,nbin2][n]): 
      temp = max(0,2*hth0[n][m+1]-hthall[n][k][m+1]) if hth0[n][m+1]<hthall[n][k][m+1] else hthall[n][k][m+1]
      haccall[n][k].SetBinError(m+1,((temp+1.)*(temp+2.)/(hth0[n][m+1]+2.)/(hth0[n][m+1]+3.)-(temp+1.)**2/(hth0[n][m+1]+2.)**2)**0.5)
    for i in range(len(enames)):
      for j in range(len(elasnames)):
        hacc[n][i][j][k].Divide(hth[n][i][j][k],hth0[n])
        for m in range([nbin1,nbin2][n]): 
          temp = max(0,2*hth0[n][m+1]-hth[n][i][j][k][m+1]) if hth0[n][m+1]<hth[n][i][j][k][m+1] else hth[n][i][j][k][m+1]
          hacc[n][i][j][k].SetBinError(m+1,((temp+1.)*(temp+2.)/(hth0[n][m+1]+2.)/(hth0[n][m+1]+3.)-(temp+1.)**2/(hth0[n][m+1]+2.)**2)**0.5)
      
for h in unnest([hth0,hthall,hth,hacc]): h.SetOption('E')

# writing
writelisto(hth0,fout,'hth0')
writelisto(hthall,fout,'hthall')
writelisto(hth[0],fout,'hth1')
writelisto(hth[1],fout,'hth2')
writelisto(haccall,fout,'accall')
writelisto(hacc[0],fout,'acc1')
writelisto(hacc[1],fout,'acc2')

fout.Close()
sys.exit()
