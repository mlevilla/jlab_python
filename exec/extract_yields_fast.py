#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','several','suffix'),
         ('b',False,'','batch_mode'),
         ('s',0,'','show_bar'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('phi',[10.,5.],'','deltaphi_cut'),
         ('e',[6.,4.],'','energy_cut'),
         ('z',[500.,150.],'','zvertex_cut'),
         ('theta',[0.6,0.7],'','theta_cut'),
         ('gem',1,'','use_gem_coordinate'),
         ('i',workf,'','input_dir'),
         ('t',True,'','write_tree'),
         ('sim',False,'','sim_file'),
         ('rdead',1.,'','radius_for_dead_module'),
         ('fiducial',1,'','fiducial_cut'),
         ('sector',[0,1,2,3,4],'','sectors included'),
         ('nf',1,'','number_of_simulation_files'),
         ('spacer',False,'','gem_spacers')]
print_help(largs)
[lrun,suffix,batch,show,outdir,limit,phicut,elascut,zcut,thetacut,gem,indir,tree,sim,rdead,fiducial,sector,nf,spacer] = catch_args(largs)

######################
## batch processing ##
######################
if batch:
  if sim=='': jsub_fast(largs,'yields',['b','r','o','s'],[['r'],get_runs_between(lrun[0],lrun[-1])],outdir)
  else: jsub_fast(largs,'simyields',['b','nf','o','s'],[['n','sim'],['_'+str(x) for x in range(nf)]],outdir)
  sys.exit()

print ' '.join(sys.argv)
date = subprocess.check_output(['date'])
print 'start: ',date[:-1]

from misc_root import *
from selection import *

if sim: f = TFile('{0}/simyields{1}.root'.format(indir,suffix[-1]))
else: f = TFile('{0}/yields_{1}{2}.root'.format(indir,lrun[0],suffix[-1]))

t = f.clusters
tree_title = t.GetTitle()
idx = tree_title.index('_')
run = int(tree_title[idx+1:])
method = tree_title[:idx]

fout = TFile(outdir+'/yields_{0}{1}.root'.format(run,suffix[0]),'recreate') if not sim else TFile(outdir+'/simyields{0}.root'.format(suffix[0]),'recreate')

load_run_params(run)
Ebeam = beam_energy[0]/1000.
alpham = (Ebeam-m_e)/(Ebeam+m_e)
max_energy = 1500. if Ebeam<2. else 2500.
nbin=35
#thbin = [0.05*i*degrad for i in range(nbin+1)]
thbin = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
q2bin = [ep_q2(th*degrad,Ebeam) if th!=0 else 0. for th in thbin]
nbin = len(thbin)-1
sector = set(sector)

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

# yield histograms
hth = listofhisto('hth',['ep','ee1','ee2','ee3'],nbin,xbin=thbin,title=';#theta (deg);')
hq2 = listofhisto('hq2_ep',[],nbin,xbin=q2bin,title=';Q^{2} (GeV^{2});')

if tree: 
  tout,vout = tree_init('clusters',[('event','i'),('n_cl','i'),('dphi','f'),('zv','f'),('E','f[n_cl]'),('E_theo','f[n_cl]'),('theta_gem','f[n_cl]'),('theta_hycal','f[n_cl]'),('x_gem','f[n_cl]'),('y_gem','f[n_cl]'),('x_hycal','f[n_cl]'),('y_hycal','f[n_cl]'),('cid','i[n_cl]')],method+'_'+str(run),2)

for _ in progress(t,show=show,n=t.GetEntries(),precision=1,limit=limit):
  # variables
  n = t.n_cl
  dphi = t.dphi
  zv = t.zv
  E = [x for x in t.E]
  E_theo = [x for x in t.E_theo]
  theta_gem = [x for x in t.theta_gem]
  theta_hycal = [x for x in t.theta_hycal]
  x_gem = [x for x in t.x_gem]
  y_gem = [x for x in t.y_gem]
  x_hycal = [x for x in t.x_hycal]
  y_hycal = [x for x in t.y_hycal]
  cid = [x for x in t.cid]
  lg = [int(x<1000) for x in t.cid]

  # ep
  if t.event==0:
    th = theta_gem[0]
    if fiducial and fiducial_cut(cid[0],[x_gem[0],y_gem[0],5198],1,1,rdead): continue
    if spacer and (gem_spacers(x_gem[0],y_gem[0]) or gem_dead(x_gem[0],y_gem[0])): continue
    if th<thetacut[1]: continue

    # temp for lg sector
    # if not is_in_sector(cid[0],sector): continue 

    elas_ep = E[0]/E_theo[0]-1
    elasc = elascut[1]*[0.024,0.062][lg[0]]/E_theo[0]**0.5
    helas1[0][0].Fill(elas_ep)
    h2cut1[0][0].Fill(th,E[0]*1000.)
    # elasticity cuts
    if abs(elas_ep)>elasc: continue
    helas1[0][1].Fill(elas_ep)
    h2cut1[0][1].Fill(th,E[0]*1000)

    q2_ep = ep_q2(th*degrad,Ebeam)
    dq2domega_ep = ep_dq2_domega(th*degrad,Ebeam)
    # barycenter
    hbarth[0].Fill(th,th)
    hbarq2.Fill(q2_ep,q2_ep)
    # yields
    hth[0].Fill(th)
    hq2.Fill(q2_ep)
    # occupancy
    for m in range(2): hxy[0][m].Fill([x_hycal,x_gem][m][0],[y_hycal,y_gem][m][0])
    # tree
    if tree:
      vout['event'][0] = 0
      vout['n_cl'][0] = 1
      vout['cid'][0] = cid[0]
      vout['E'][0] = E[0]
      vout['dphi'][0] = 0
      vout['zv'][0] = 0
      vout['x_gem'][0] = x_gem[0]
      vout['y_gem'][0] = y_gem[0]
      vout['x_hycal'][0] = x_hycal[0]
      vout['y_hycal'][0] = y_hycal[0]
      vout['E_theo'][0] = E_theo[0]
      vout['theta_hycal'][0] = theta_hycal[0]
      vout['theta_gem'][0] = theta_gem[0]
      tout.Fill()

  # moller

  # ee1
  if t.event==1:
    th = theta_gem[0]
    if fiducial and fiducial_cut(cid[0],[x_gem[0],y_gem[0],5198],1,1,rdead): continue
    if spacer and (gem_spacers(x_gem[0],y_gem[0]) or gem_dead(x_gem[0],y_gem[0])): continue
    if th<thetacut[1]: continue

    # temp for lg sector
    # if not is_in_sector(cid[0],sector): continue 

    elas_m = E[0]/E_theo[0]-1
    elasc = elascut[1]*[0.024,0.062][lg[0]]/E_theo[0]**0.5
    helas1[1][0].Fill(elas_m)
    h2cut1[1][0].Fill(th,E[0]*1000)
    # cuts
    if abs(elas_m)>elasc: continue
    helas1[1][1].Fill(elas_m)
    h2cut1[1][1].Fill(th,E[0]*1000)
    # yields
    hth[1].Fill(th)
    # barycenter
    hbarth[1].Fill(th,th)
    # occupancy
    for m in range(2): hxy[0][m].Fill([x_hycal,x_gem][m][0],[y_hycal,y_gem][m][0])
    # tree
    if tree:
      vout['event'][0] = 0
      vout['n_cl'][0] = 1
      vout['cid'][0] = cid[0]
      vout['E'][0] = E[0]
      vout['dphi'][0] = 0
      vout['zv'][0] = 0
      vout['x_gem'][0] = x_gem[0]
      vout['y_gem'][0] = y_gem[0]
      vout['x_hycal'][0] = x_hycal[0]
      vout['y_hycal'][0] = y_hycal[0]
      vout['E_theo'][0] = E_theo[0]
      vout['theta_hycal'][0] = theta_hycal[0]
      vout['theta_gem'][0] = theta_gem[0]
      tout.Fill()
      
  # ee2, ee3
  for k in range(2):
    if t.event!=k+2: continue
    if fiducial and any(fiducial_cut(cid[j],[x_gem[j],y_gem[j],5198],1,1,rdead) for j in range(2)): continue
    if spacer and any(gem_spacers(x_gem[j],y_gem[j]) or gem_dead(x_gem[j],y_gem[j]) for j in range(2)): continue
    th = [[theta_hycal[j],theta_gem[j]] for j in range(2)]
    if any(th[j][1-k]<thetacut[1-k] for j in range(n)): continue
    dE = E[0]+E[1]-Ebeam
    elas_sum = dE/Ebeam
    dEcut = elascut[1-k]*(([0.024,0.062][lg[0]]*E[0]**0.5)**2+([0.024,0.062][lg[1]]*E[1]**0.5)**2)**0.5
    # cuts
    hdphi[k][0].Fill(dphi)
    helas2[k][0].Fill(elas_sum)
    hzvertex[k][0].Fill(zv)
    for j in range(2): h2cut2[k][0].Fill(th[j][1],E[j]*1000)
    b1 = abs(dphi)<phicut[1-k]
    b2 =  abs(zv)<zcut[1-k]
    b3 = abs(dE)<dEcut
    if b1: 
      for j in range(2): h2cut2[k][1].Fill(th[j][1],E[j]*1000)
      hzvertex[k][1].Fill(zv)
      helas2[k][1].Fill(elas_sum)
      if b2: 
        hdphi[k][1].Fill(dphi)
        helas2[k][2].Fill(elas_sum)
      if b3:
        hzvertex[k][2].Fill(zv)
        hdphi[k][2].Fill(dphi)
      if b1 and b2:
        for j in range(2): h2cut2[k][2].Fill(th[j][1],E[j]*1000)
        helas2[k][3].Fill(elas_sum)
      if b1 and b3: hzvertex[k][3].Fill(zv)
      if b2 and b3: hdphi[k][3].Fill(dphi)
      if not (b1 and b2 and b3): continue
      helas2[k][4].Fill(elas_sum)
      hdphi[k][4].Fill(dphi)
      hzvertex[k][4].Fill(zv)
      for j in range(2):
        if th[j][1]<thetacut[1]: continue
        h2cut2[k][3].Fill(th[j][1],E[j]*1000)
        elas_m = E[j]/E_theo[j]-1
        elasc = elascut[1]*[0.024,0.062][lg[j]]/(E_theo[j])**0.5
        helas1[2+k][0].Fill(elas_m)
        # cuts
        if abs(elas_m)>elasc: continue
        helas1[2+k][1].Fill(elas_m)
        h2cut2[k][4].Fill(th[j][1],E[j]*1000)
        # yields
        hth[2+k].Fill(th[j][1])
        # barycenter
        hbarth[2+k].Fill(th[j][1],th[j][1])
        # occupancy
        for m in range(2): hxy[0][m].Fill([x_hycal,x_gem][m][j],[y_hycal,y_gem][m][j])
      # tree
      if tree:
        vout['event'][0] = k+2
        vout['n_cl'][0] = 2
        vout['dphi'][0] = dphi
        vout['zv'][0] = zv
        for j in range(2):
          vout['E'][j] = E[j]
          vout['E_theo'][j] = E_theo[j]
          vout['theta_hycal'][j] = theta_hycal[j]
          vout['theta_gem'][j] = theta_gem[j]
          vout['x_hycal'][j] = x_hycal[j]
          vout['y_hycal'][j] = y_hycal[j]
          vout['x_gem'][j] = x_gem[j]
          vout['y_gem'][j] = y_gem[j]
          vout['cid'][j] = cid[j]
        tout.Fill()
        
hbarq2.Divide(hq2)
for i in range(4): hbarth[i].Divide(hth[i])
    
writelisto([helas1,helas2,hdphi,hzvertex,h2cut1,h2cut2],fout,['selection'])
writelisto(hth,fout,['yields'])
writelisto([hbarth,hbarq2],fout,['barycenter'])
writelisto(hxy,fout,['occupancy'])

if tree:
  fout.cd()
  tout.Write()
fout.Close()

date = subprocess.check_output(['date'])
print 'end: ',date[:-1]
