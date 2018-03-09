#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','','suffix'),
         ('m','island','','clustering_method'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('phi',[10.,5.],'','deltaphi_cut'),
         ('e',[6.,4.],'','energy_cut'),
         ('theta',[0.6,0.7],'','theta_cut'),
         ('i',workf+'/trees_sel5','','input_dir'),
         ('sim','','','simulation_file_name'),
         ('ebeam',2,'','beam_energy_for_simulation'),
         ('rdead',1.,'','radius_for_dead_module'),
         ('fiducial',True,'','fiducial_cut'),
         ('spacer',False,'','gem_spacer_and_dead_area'),
         ('nf',1,'','number_of_simulation_files'),
         ('gem',[0,1],'','use_gem_coordinate')]
aparser = ArgParser(largs,sub=1)
[lrun,suffix,method,outdir,limit,phicut,elascut,thetacut,indir,sim,ebeam,rdead,fiducial,spacer,nf,gem] = aparser.argl

######################
## batch processing ##
######################
if batch[0]:
  if sim=='': jsub_fast(aparser,'shtree'+suffix,['b','r','o','s'],[['r'],get_runs_between(lrun[0],lrun[-1])],outdir,time=wtime[0])
  else: jsub_fast(aparser,'shtree'+suffix,['b','nf','o','s'],[['n','sim'],['_'+str(x) for x in range(nf)]],outdir,time=wtime[0])
  sys.exit()

####################
## initialization ##
####################
print ' '.join(sys.argv)
date = subprocess.check_output(['date'])
print 'start: ',date[:-1]

from misc_root import *
from selection import *

if sim!='':
  f = TFile(sim+'_rec.root')
  t = f.T
  run = 1288 if ebeam==1 else 1443
  t.SetBranchStatus("*",0)
  for x in ['GEM.X','GEM.Y','GEM.Z','HC.N','HC.X','HC.Y','HC.Z','HC.P','HC.CID']: t.SetBranchStatus(x,1)
  fout = TFile(outdir+'/zvertex_{0}{1}.root'.format(suffix),'recreate')
else:
  run = lrun[0]
  f = TFile(indir+'/tree_{0}_{1}.root'.format(method,run))
  t = f.event
  t.SetBranchStatus("*",0)
  for x in ['iev','Ebeam','n_cl','id','nh','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
  fout = TFile(outdir+'/zvertex_{0}{1}.root'.format(run,suffix),'recreate')

load_run_params(run)
load_positions(run=run)
load_dead_pos(run=run)
if sim=='': load_run_calib(run)
Ebeam = beam_energy[0]/1000.

tree_ee2, v_ee2 = tree_init('ee2',[('iev','i'),('dphi','f'),('zv1','f'),('zv2','f'),('dE','f'),('n_cl','i'),('icl','i[n_cl]'),('E','f[n_cl]'),('x_hycal','f[n_cl]'),('y_hycal','f[n_cl]'),('z_hycal','f[n_cl]'),('id_hycal','i[n_cl]'),('theta_hycal','f[n_cl]'),('phi_hycal','f[n_cl]'),('x_gem','f[n_cl]'),('y_gem','f[n_cl]'),('z_gem','f[n_cl]'),('theta_gem','f[n_cl]'),('phi_gem','f[n_cl]'),('E_theo','f[n_cl]')],'ee2',2)


thetacut = [thetacut[g] for g in gem]
phicut = [phicut[g] for g in gem]

################
## event loop ##
################
for i,_ in enumerate(progress(t,show=show_progress[0],n=t.GetEntries(),precision=2,limit=limit)):
  if sim=='': 
    [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=fiducial,exclude_dead=0,match=0,rdead=rdead,spacer=spacer)
  else: 
    [E,c,idx,lg,cid,_] = get_variables_sim(t,lincorr=0,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=any(gem),rdead=rdead,spacer=spacer)

  for j1 in range(len(idx)-1):
    for j2 in range(j1+1,len(idx)):
      
      flag = is_ee2(E=[E[j1],E[j2]],Ebeam=Ebeam,c=[c[j1],c[j2]],lg=[lg[j1],lg[j2]],sigma=[elascut[1],elascut[gem[1]]],thetacut=[thetacut[gem[0]],thetacut[gem[1]]],gem=gem,phicut=phicut,zcut=[10000,10000])
      if not all(flag): continue
      v_ee2['iev'][0] = t.iev if sim=='' else i
      v_ee2['dphi'][0] = ((fphi(c[j1][gem[0]])-fphi(c[j2][gem[0]]))%(2*pi)-pi)/degrad
      cproj = [[x[gem[0]][m]*zhycal/x[gem[0]][2] for m in range(2)] for x in [c[j1],c[j2]]]
      rj1 = [(x[0]**2+x[1]**2)**0.5 for x in cproj]
      v_ee2['zv1'][0] = ((m_e+Ebeam)*rj1[0]*rj1[1]/2./m_e)**0.5-zhycal
      
      rj2 = [(c[j][gem[0]][0]**2+c[j][gem[0]][1]**2)**0.5 for j in [j1,j2]]
      dz = c[j1][gem[0]][2] - c[j2][gem[0]][2]
      v_ee2['zv2'][0] = ((m_e+Ebeam)*rj2[0]*rj2[1]/2./m_e+dz**2/4)**0.5 - c[j1][gem[0]][2] + dz/2.
      
      v_ee2['dE'][0] = E[j1]+E[j2]-Ebeam*1000.
      n_cl = 0
      for k,j in enumerate([j1,j2]):
        if not flag[k+1]: continue
        v_ee2['icl'][k] = idx[j]
        v_ee2['E'][k] = E[j]
        v_ee2['x_hycal'][k] = c[j][0][0]
        v_ee2['y_hycal'][k] = c[j][0][1]
        v_ee2['z_hycal'][k] = c[j][0][2]
        v_ee2['id_hycal'][k] = cid[j]
        v_ee2['theta_hycal'][k] = ftheta(c[j][0])/degrad
        v_ee2['phi_hycal'][k] = fphi(c[j][0])/degrad
        v_ee2['x_gem'][k] = c[j][gem[1]][0]
        v_ee2['y_gem'][k] = c[j][gem[1]][1]
        v_ee2['z_gem'][k] = c[j][gem[1]][2]
        v_ee2['theta_gem'][k] = ftheta(c[j][gem[1]])/degrad
        v_ee2['phi_gem'][k] = fphi(c[j][gem[1]])/degrad
        v_ee2['E_theo'][k] = moller_energy(v_ee2['theta_gem'][k]*degrad,Ebeam)
        n_cl+=1
      v_ee2['n_cl'][0] = n_cl
      tree_ee2.Fill()

fout.cd()
tree_ee2.Write()
fout.Close()
