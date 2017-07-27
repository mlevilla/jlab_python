#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','','suffix'),
         ('b',False,'','batch_mode'),
         ('m','island','','clustering_method'),
         ('s',0,'','show_bar'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('i',workf+'/trees_sel/island','','input_dir'),
         ('t',True,'','write_tree'),
         ('count',False,'','counter'),
         ('gem',1,'','use_gem_coordinate'),
         ('phi',[10.,5.],'','deltaphi_cut'), # cut
         ('e',[6.,4.],'','energy_cut'), # cut
         ('z',[500.,150.],'','zvertex_cut'), # cut
         ('theta',[0.6,0.7],'','theta_cut'), # cut
         ('rdead',1.,'','radius_for_dead_module'), # cut
         ('fiducial',True,'','fiducial_cut'), # cut
         ('spacer',False,'','gem_spacer_and_dead_area'), # cut
         ('sim','','several','simulation_file_name'), # simulation
         ('ebeam',2,'','beam_energy_for_simulation'), # simulation
         ('nf',1,'','number_of_simulation_files'), # simulation
         ('lum',[1.,1.],'','luminosity (moller,ep)') # simulation
       ]
print_help(largs)
[lrun,suffix,batch,method,show,outdir,limit,indir,tree,cntv,gem,phicut,elascut,zcut,thetacut,rdead,fiducial,spacer,sim,ebeam,nf,luminosity] = catch_args(largs)

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
from counter import *

if sim!=['']:
  t = TChain('T')
  for x in sim: t.AddFile(x+'_rec.root')
  run = 1288 if ebeam==1 else 1443
  t.SetBranchStatus("*",0)
  for x in ['GEM.X','GEM.Y','GEM.Z','HC.N','HC.X','HC.Y','HC.Z','HC.P','HC.CID']: t.SetBranchStatus(x,1)
  fout = TFile(outdir+'/simyields{0}.root'.format(suffix),'recreate')
else:
  run = lrun[0]
  f = TFile(indir+'/tree_{0}_{1}.root'.format(method,run))
  t = f.event
  t.SetBranchStatus("*",0)
  for x in ['Ebeam','n_cl','id','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
  fout = TFile(outdir+'/yields_{0}{1}.root'.format(run,suffix),'recreate')

load_run_params(run)
Ebeam = beam_energy[0]/1000.
max_energy = 1500. if Ebeam<2. else 2500.
nbin = 35
#thbin = [0.05*i*degrad for i in range(nbin+1)]
thbin = [0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
#thbin = [0.5,0.55,0.6,0.65,0.7,0.715,0.731,0.748,0.766,0.785,0.806,0.828,0.852,0.879,0.908,0.94,0.975,1.014,1.057,1.105,1.157,1.211,1.27,1.338,1.417,1.514,1.634,1.787,2.0,2.213,2.492,2.792,3.092,3.392,3.692,3.992,4.292,4.592,4.892,5.2]
q2bin = [ep_q2(th*degrad,Ebeam) if th!=0 else 0. for th in thbin]
nbin = len(thbin)-1

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

# counter
counter = Counter([0.7+0.5*i for i in range(15)],['ep_raw','ep_match','ep_theta','ep_elas','ee1_raw','ee1_match','ee1_theta','ee1_elas','ee2_raw','ee2_match','ee2_theta','ee2_phi','ee2_zvertex','ee2_elas2','ee2_raw1','ee2_match1','ee2_theta1','ee2_elas1','ee3_raw','ee3_match','ee3_theta','ee3_phi','ee3_zvertex','ee3_elas2','ee3_raw1','ee3_match1','ee3_theta1','ee3_elas1'],cntv)

if tree: 
  tout,vout = tree_init('clusters',[('event','i'),('n_cl','i'),('dphi','f'),('zv','f'),('E','f[n_cl]'),('E_theo','f[n_cl]'),('theta_gem','f[n_cl]'),('theta_hycal','f[n_cl]'),('x_gem','f[n_cl]'),('y_gem','f[n_cl]'),('x_hycal','f[n_cl]'),('y_hycal','f[n_cl]'),('cid','i[n_cl]')],method+'_'+str(run),2)

ratio_lum = luminosity[0]/luminosity[1]+0.5
print ratio_lum

for _ in progress(t,show=show,n=t.GetEntries(),precision=2,limit=ratio_lum):

  # variables
  if sim=='': 
    [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=0,rdead=rdead,spacer=spacer)
  else: 
    [E,c,idx,lg,cid,_] = get_variables_sim(t,lincorr=0,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=0,rdead=rdead,spacer=spacer)

  theta = [[ftheta(y) if y!=[] else 0 for y in x] for x in c]
  thetadeg = [[abs(degrees(y)) for y in x] for x in theta]

  # ep
  for j in range(len(idx)):
    th = thetadeg[j][gem]
    counter.count('ep_raw',th)
    if counter.cut('ep_match',c[j][gem]!=[],th): continue
    if counter.cut('ep_theta',th>thetacut[1],th): continue
    E_ep = [ep_energy_el(theta[j][i],Ebeam) for i in range(2)]
    elas_ep = E[j]/E_ep[gem]/1000-1
    elasc = elascut[1]*[0.024,0.062][lg[j]]/E_ep[gem]**0.5
    helas1[0][0].Fill(elas_ep)
    h2cut1[0][0].Fill(th,E[j])
    # elasticity cut
    if counter.cut('ep_elas',abs(elas_ep)<elasc,th): continue
    helas1[0][1].Fill(elas_ep)
    h2cut1[0][1].Fill(th,E[j])

    q2_ep = ep_q2(theta[j][gem],Ebeam)
    dq2domega_ep = ep_dq2_domega(theta[j][gem],Ebeam)
    # barycenter
    hbarth[0].Fill(th,th)
    hbarq2.Fill(q2_ep,q2_ep)
    # yields
    hth[0].Fill(th)
    hq2.Fill(q2_ep)
    # occupancy
    for m in range(2): hxy[0][m].Fill(c[j][m][0],c[j][m][1])
    # tree
    if tree:
      vout['event'][0] = 0
      vout['n_cl'][0] = 1
      vout['cid'][0] = cid[j]
      vout['E'][0] = E[j]/1000.
      vout['dphi'][0] = 0
      vout['zv'][0] = 0
      vout['x_gem'][0] = c[j][1][0]
      vout['y_gem'][0] = c[j][1][1]
      vout['x_hycal'][0] = c[j][0][0]
      vout['y_hycal'][0] = c[j][0][1]
      vout['E_theo'][0] = E_ep[gem]
      vout['theta_hycal'][0] = thetadeg[j][0]
      vout['theta_gem'][0] = thetadeg[j][1]
      tout.Fill()
      
  # moller
  Em = [[moller_energy(y,Ebeam) if y!=0 else 0 for y in x] for x in theta]

  # ee1
  for j in range(len(idx)):
    th = thetadeg[j][gem]
    counter.count('ee1_raw',th)
    if counter.cut('ee1_match',c[j][gem]!=[],th): continue
    if counter.cut('ee1_theta',th>thetacut[1],th): continue
    elas_m = E[j]/Em[j][gem]/1000.-1
    elasc = elascut[1]*[0.024,0.062][lg[j]]/(Em[j][gem])**0.5
    helas1[1][0].Fill(elas_m)
    h2cut1[1][0].Fill(th,E[j])
    # cuts
    if counter.cut('ee1_elas',abs(elas_m)<elasc,th): continue
    helas1[1][1].Fill(elas_m)
    h2cut1[1][1].Fill(th,E[j])
    # yields
    hth[1].Fill(th)
    # barycenter
    hbarth[1].Fill(th,th)
    # occupancy
    for m in range(2): hxy[1][m].Fill(c[j][m][0],c[j][m][1])
    # tree
    if tree:
      vout['event'][0] = 1
      vout['n_cl'][0] = 1
      vout['cid'][0] = cid[j]
      vout['E'][0] = E[j]/1000.
      vout['dphi'][0] = 0
      vout['zv'][0] = 0
      vout['x_gem'][0] = c[j][1][0]
      vout['y_gem'][0] = c[j][1][1]
      vout['x_hycal'][0] = c[j][0][0]
      vout['y_hycal'][0] = c[j][0][1]
      vout['E_theo'][0] = Em[j][gem]
      vout['theta_hycal'][0] = thetadeg[j][0]
      vout['theta_gem'][0] = thetadeg[j][1]
      tout.Fill()
                   
  # ee2, ee3
  for k in range(2):
    dEmin, jmin =10000, [-1,-1]
    for j1 in range(len(idx)-1):
      for j2 in range(j1+1,len(idx)):
        ths = [thetadeg[j1][1-k],thetadeg[j2][1-k]]
        counter.count('ee'+str(k+2)+'_raw',ths)
        if counter.cut('ee'+str(k+2)+'_match',all(c[j][1-k]!=[] for j in [j1,j2]),ths): continue
        if counter.cut('ee'+str(k+2)+'_theta',all(th>thetacut[1-k] for th in ths),ths): continue
        dphi = degrees((fphi(c[j1][1-k])-fphi(c[j2][1-k]))%(2*pi)-pi)
        dE = E[j1]+E[j2]-Ebeam*1e3
        elas_sum = dE/Ebeam/1e3
        dEcut = elascut[1]*(([0.024,0.062][lg[j1]]*E[j1]**0.5)**2+([0.024,0.062][lg[j2]]*E[j2]**0.5)**2)**0.5*1000**0.5
        # rj = [(c[j][1-k][0]**2+c[j][1-k][1]**2)**0.5 for j in [j1,j2]]
        # dz = c[j1][1-k][2] - c[j2][1-k][2]
        # zvertex = ((m_e+Ebeam)*rj[0]*rj[1]/2./m_e+dz**2/4)**0.5 - c[0][1-k][2] + dz/2.
        cproj = [[c[j][1-k][m]*hycal_center[2]/c[j][1-k][2] for m in range(2)] for j in [j1,j2]]
        rj2 = [(x[0]**2+x[1]**2)**0.5 for x in cproj]
        zvertex = ((m_e+Ebeam)*rj2[0]*rj2[1]/2./m_e)**0.5-hycal_center[2]
        # cuts
        hdphi[k][0].Fill(dphi)
        helas2[k][0].Fill(elas_sum)
        hzvertex[k][0].Fill(zvertex)
        for j in [j1,j2]: h2cut2[k][0].Fill(thetadeg[j][1-k],E[j])
        b1 = abs(dphi)<phicut[1-k]
        b2 =  abs(zvertex)<zcut[1-k]
        b3 = abs(dE)<dEcut
        if b1: 
          for j in [j1,j2]: h2cut2[k][1].Fill(thetadeg[j][1-k],E[j])
          hzvertex[k][1].Fill(zvertex)
          helas2[k][1].Fill(elas_sum)
          counter.count('ee'+str(k+2)+'_phi',ths)
        if b2: 
          hdphi[k][1].Fill(dphi)
          helas2[k][2].Fill(elas_sum)
        if b3:
          hzvertex[k][2].Fill(zvertex)
          hdphi[k][2].Fill(dphi)
        if b1 and b2:
          counter.count('ee'+str(k+2)+'_zvertex',ths)
          for j in [j1,j2]: h2cut2[k][2].Fill(thetadeg[j][1-k],E[j])
          helas2[k][3].Fill(elas_sum)
        if b1 and b3: hzvertex[k][3].Fill(zvertex)
        if b2 and b3: hdphi[k][3].Fill(dphi)
        if counter.cut('ee'+str(k+2)+'_elas2',b1 and b2 and b3,ths): continue
        helas2[k][4].Fill(elas_sum)
        hdphi[k][4].Fill(dphi)
        hzvertex[k][4].Fill(zvertex)
        for j in [j1,j2]: h2cut2[k][3].Fill(thetadeg[j][1-k],E[j])
        if dE<dEmin:
          dEmin = dE
          jmin = [j1,j2]
          vmin = [elas_sum,dphi,zvertex]
    if dEmin!=10000:
      for j in jmin:
        th = thetadeg[j][gem]
        counter.count('ee'+str(k+2)+'_raw1',th)
        if counter.cut('ee'+str(k+2)+'_match1',c[j][gem]!=[],th): continue
        if counter.cut('ee'+str(k+2)+'_theta1',th>thetacut[1],th): continue
        elas_m = E[j]/Em[j][gem]/1000.-1
        elasc = elascut[1]*[0.024,0.062][lg[j]]/(Em[j][gem])**0.5
        helas1[2+k][0].Fill(elas_m)
        # cuts
        if counter.cut('ee'+str(k+2)+'_elas1',abs(elas_m)<elasc,th): continue
        helas1[2+k][1].Fill(elas_m)
        h2cut2[k][4].Fill(th,E[j])
        # yields
        hth[2+k].Fill(th)
        # barycenter
        hbarth[2+k].Fill(th,th)
        # occupancy
        for m in range(2): hxy[2+k][m].Fill(c[j][m][0],c[j][m][1])
      # tree
      if tree:
        vout['event'][0] = k+2
        vout['n_cl'][0] = 2
        vout['dphi'][0] = vmin[1]
        vout['zv'][0] = vmin[2]
        for m,j in enumerate(jmin):
          vout['cid'][m] = cid[j]
          vout['E'][m] = E[j]/1000.
          if c[j][1]:
            vout['x_gem'][m] = c[j][1][0]
            vout['y_gem'][m] = c[j][1][1]
          else:
            vout['x_gem'][m] = -10000
            vout['y_gem'][m] = -10000
          vout['x_hycal'][m] = c[j][0][0]
          vout['y_hycal'][m] = c[j][0][1]
          vout['E_theo'][m] = Em[j][gem]
          vout['theta_hycal'][m] = thetadeg[j][0]
          vout['theta_gem'][m] = thetadeg[j][1]
        tout.Fill()
     
hbarq2.Divide(hq2)
for i in range(4): hbarth[i].Divide(hth[i])
   
if sim=='': counter.write(outdir+'/yields_{0}{1}_counter.txt'.format(run,suffix))
else: counter.write(outdir+'/simyields{0}_counter.txt'.format(suffix))
 
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
