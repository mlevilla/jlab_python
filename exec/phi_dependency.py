#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','','suffix'),
         ('b',False,'','batch_mode'),
         ('m','island','','clustering_method'),
         ('s',0,'','show_bar'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('phi',[10.,5.],'','deltaphi_cut'),
         ('e',4.,'','energy_cut'),
         ('z',[500.,150.],'','zvertex_cut'),
         ('theta',[0.6,0.7],'','theta_cut'),
         ('i',workf+'/trees_sel/island','','input_dir'),
         ('sim','','','simulation_file_name'),
         ('ebeam',2,'','beam_energy_for_simulation'),
         ('rdead',1.,'','radius_for_dead_module'),
         ('fiducial',True,'','fiducial_cut'),
         ('spacer',False,'','gem_spacer_and_dead_area'),
         ('nf',1,'','number_of_simulation_files')]
print_help(largs)
[lrun,suffix,batch,method,show,outdir,limit,phicut,elascut,zcut,thetacut,indir,sim,ebeam,rdead,fiducial,spacer,nf] = catch_args(largs)

######################
## batch processing ##
######################
if batch:
  if sim=='': jsub_fast(largs,'phidep',['b','r','o','s'],[['r'],get_runs_between(lrun[0],lrun[-1])],outdir)
  else: jsub_fast(largs,'phidep',['b','nf','o','s'],[['n','sim'],['_'+str(x) for x in range(nf)]],outdir)
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
  for x in ['Ebeam','n_cl','id','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
  
load_run_params(run)
Ebeam = beam_energy[0]/1000.

fout = TFile(outdir+'/phidep_{0}{1}.root'.format(run,suffix),'recreate')
hphi = listofhisto('hphi',range(35),16,-180,180)

################
## event loop ##
################
for _ in progress(t,show=show,n=t.GetEntries(),precision=2,limit=limit):
  
  if sim=='': 
    [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=0,rdead=rdead,spacer=spacer)
  else: 
    [E,c,idx,lg,cid,_] = get_variables_sim(t,lincorr=0,exclude_edge=fiducial,exclude_dead=(rdead!=0),match=0,rdead=rdead,spacer=spacer)

  for j in range(len(idx)):
    if  not is_ep(E=E[j],Ebeam=Ebeam,c=c[j][0],lg=lg[j],sigma=elascut,thetacut=0.8): continue
    phi = fphi(c[j][0])/degrad
    theta = ftheta(c[j][0])/degrad
    itheta = int((theta-0.8)/0.1)
    if itheta>=35: continue
    hphi[itheta].Fill(phi)

#############
## writing ##
#############
for h in hphi: 
  if h.GetEntries()!=0: h.Scale(16./h.GetEntries())

writelisto(hphi,fout)
fout.Close()
