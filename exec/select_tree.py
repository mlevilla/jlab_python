#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','','suffix'),
         ('m','island','','clustering_method'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('phi',[10.,5.],'','deltaphi_cut'),
         ('e',[6.,4.],'','energy_cut'),
         ('z',[500.,150.],'','zvertex_cut'),
         ('theta',[0.6,0.7],'','theta_cut'),
         ('i',workf+'/trees_sel/island','','input_dir'),
         ('sim','','','simulation_file_name'),
         ('ebeam',2,'','beam_energy_for_simulation'),
         ('rdead',1.,'','radius_for_dead_module'),
         ('fiducial',True,'','fiducial_cut'),
         ('spacer',False,'','gem_spacer_and_dead_area'),
         ('nf',1,'','number_of_simulation_files'),
         ('gem',[0,1],'','use_gem_coordinate')]
parser = ArgParser(largs)
[lrun,suffix,method,outdir,limit,phicut,elascut,zcut,thetacut,indir,sim,ebeam,rdead,fiducial,spacer,nf,gem] = parser.argl

######################
## batch processing ##
######################
if batch[0]:
  if sim=='': jsub_fast(largs,'shtree'+suffix,['b','r','o','s'],[['r'],get_runs_between(lrun[0],lrun[-1])],outdir)
  else: jsub_fast(largs,'shtree'+suffix,['b','nf','o','s'],[['n','sim'],['_'+str(x) for x in range(nf)]],outdir)
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
else:
  run = lrun[0]
  f = TFile(indir+'/tree_{0}_{1}.root'.format(method,run))
  t = f.event
  t.SetBranchStatus("*",0)
  for x in ['iev','Ebeam','n_cl','id','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
  
load_run_params(run)
Ebeam = beam_energy[0]/1000.

fout = TFile(outdir+'/shtree_{0}{1}.root'.format(run,suffix),'recreate')
tree_ep, v_ep = tree_init('ep',[('iev','i'),('icl','i'),('E','f'),('x_hycal','f'),('y_hycal','f'),('z_hycal','f'),('id_hycal','i'),('theta_hycal','f'),('phi_hycal','f'),('x_gem','f'),('y_gem','f'),('z_gem','f'),('theta_gem','f'),('phi_gem','f'),('E_theo','f'),('Q2','f')],'ep')
tree_ee1, v_ee1 = tree_init('ee1',[('iev','i'),('icl','i'),('E','f'),('x_hycal','f'),('y_hycal','f'),('z_hycal','f'),('id_hycal','i'),('theta_hycal','f'),('phi_hycal','f'),('x_gem','f'),('y_gem','f'),('z_gem','f'),('theta_gem','f'),('phi_gem','f'),('E_theo','f')],'ee1')
tree_ee2, v_ee2 = tree_init('ee2',[('iev','i'),('dphi','f'),('zv','f'),('dE','f'),('n_cl','i'),('icl','i[n_cl]'),('E','f[n_cl]'),('x_hycal','f[n_cl]'),('y_hycal','f[n_cl]'),('z_hycal','f[n_cl]'),('id_hycal','i[n_cl]'),('theta_hycal','f[n_cl]'),('phi_hycal','f[n_cl]'),('x_gem','f[n_cl]'),('y_gem','f[n_cl]'),('z_gem','f[n_cl]'),('theta_gem','f[n_cl]'),('phi_gem','f[n_cl]'),('E_theo','f[n_cl]')],'ee2',2)

thetacut = [thetacut[g] for g in gem]
phicut = [phicut[g] for g in gem]
zcut = [zcut[g] for g in gem]

################
## event loop ##
################
for i,_ in enumerate(progress(t,show=show_progress[0],n=t.GetEntries(),precision=2,limit=limit)):
  if sim=='': 
    [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=any(gem),rdead=rdead,spacer=spacer)
  else: 
    [E,c,idx,lg,cid,_] = get_variables_sim(t,lincorr=0,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=any(gem),rdead=rdead,spacer=spacer)

  for j in range(len(idx)):
    if is_ep(E=E[j],Ebeam=Ebeam,c=c[j][gem[1]],lg=lg[j],sigma=elascut[1],thetacut=thetacut[1]): 
      v_ep['iev'][0] = t.iev if sim=='' else i
      v_ep['icl'][0] = idx[j]
      v_ep['E'][0] = E[j]
      v_ep['x_hycal'][0] = c[j][0][0]
      v_ep['y_hycal'][0] = c[j][0][1]
      v_ep['z_hycal'][0] = c[j][0][2]
      v_ep['id_hycal'][0] = cid[j]
      v_ep['theta_hycal'][0] = ftheta(c[j][0])/degrad
      v_ep['phi_hycal'][0] = fphi(c[j][0])/degrad
      if c[j][1]:
        v_ep['x_gem'][0] = c[j][1][0]
        v_ep['y_gem'][0] = c[j][1][1]
        v_ep['z_gem'][0] = c[j][1][2]
        v_ep['theta_gem'][0] = ftheta(c[j][1])/degrad
        v_ep['phi_gem'][0] = fphi(c[j][1])/degrad
      else:
        v_ep['x_gem'][0] = -10000
        v_ep['y_gem'][0] = -10000
        v_ep['z_gem'][0] = -10000
        v_ep['theta_gem'][0] = -10000
        v_ep['phi_gem'][0] = -10000
      v_ep['E_theo'][0] = ep_energy_el(v_ep['theta_gem'][0]*degrad,Ebeam)
      v_ep['Q2'][0] = ep_q2(v_ep['theta_gem'][0]*degrad,Ebeam)
      tree_ep.Fill()

    if is_ee1(E=E[j],Ebeam=Ebeam,c=c[j][gem[1]],lg=lg[j],sigma=elascut[gem[1]],thetacut=thetacut[1]): 
      v_ee1['iev'][0] = t.iev if sim=='' else i
      v_ee1['icl'][0] = idx[j]
      v_ee1['E'][0] = E[j]
      v_ee1['x_hycal'][0] = c[j][0][0]
      v_ee1['y_hycal'][0] = c[j][0][1]
      v_ee1['z_hycal'][0] = c[j][0][2]
      v_ee1['id_hycal'][0] = cid[j]
      v_ee1['theta_hycal'][0] = ftheta(c[j][0])/degrad
      v_ee1['phi_hycal'][0] = fphi(c[j][0])/degrad
      v_ee1['x_gem'][0] = c[j][gem[1]][0]
      v_ee1['y_gem'][0] = c[j][gem[1]][1]
      v_ee1['z_gem'][0] = c[j][gem[1]][2]
      v_ee1['theta_gem'][0] = ftheta(c[j][gem[1]])/degrad
      v_ee1['phi_gem'][0] = fphi(c[j][gem[1]])/degrad
      v_ee1['E_theo'][0] = moller_energy(v_ee1['theta_gem'][0]*degrad,Ebeam)
      tree_ee1.Fill()

  for j1 in range(len(idx)-1):
    for j2 in range(j1+1,len(idx)):
      
      flag = is_ee2(E=[E[j1],E[j2]],Ebeam=Ebeam,c=[c[j1],c[j2]],lg=[lg[j1],lg[j2]],sigma=[elascut[1],elascut[gem[1]]],thetacut=[thetacut[gem[0]],thetacut[gem[1]]],gem=gem,phicut=phicut,zcut=zcut)
      if not flag[0]: continue
      v_ee2['iev'][0] = t.iev if sim=='' else i
      v_ee2['dphi'][0] = ((fphi(c[j1][gem[0]])-fphi(c[j2][gem[0]]))%(2*pi)-pi)/degrad
      cproj = [[x[gem[0]][m]*zhycal/x[gem[0]][2] for m in range(2)] for x in [c[j1],c[j2]]]
      rj = [(x[0]**2+x[1]**2)**0.5 for x in cproj]
      v_ee2['zv'][0] = ((m_e+Ebeam)*rj[0]*rj[1]/2./m_e)**0.5-zhycal
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
tree_ep.Write()
tree_ee1.Write()
tree_ee2.Write()
fout.Close()
