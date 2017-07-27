#!/usr/bin/env python
from params import *

largs = [('r',[-1,'-1'],'','run_number'),
         ('i',work+'/trees_calib/island','','input_dir'),
         ('o',work,'','output_folder'),
         ('m','island','','clustering_method'),
         ('n','','','suffix'),
         ('s',0,'','show_bar'),
         ('l',1.0,'','limit'),
         ('b',False,'','batch_mode'),
         ('a',False,'','merge runs'),
         ('d','gem','','coordinate_reconstruction_detector')]
print_help(largs)
[lrun,indir,outdir,method,suffix,show,limit,batch,hadd,det] = [catch_arg(x,y,z) for x,y,z,c in largs]

lrun[1] = lrun[0] if lrun[1][0]=='-' else int(lrun[1])

#######################################################################
################ batch jobs ###########################################
#######################################################################
if (not hadd) and batch: 
  input_files = ' '.join([indir+'/tree_'+method+'_'+str(run)+'.root' for run in get_runs_between(lrun[0],lrun[1],'prod')])
  jsub(project='prad',track='analysis',jobname='histo{0}_{1}_{2}'.format(suffix,lrun[0],lrun[1]),command=os.path.abspath('histos_from_tree_prod.py'),options='-o . -n '+suffix+' -d '+det,memory='2 GB',input_files=input_files,input_data='tree.root',output_data='histo'+suffix+'*.root',output_template=outdir+'/@OUTPUT_DATA@')
  sys.exit()
  

#######################################################################
################ merging histograms ###################################
#######################################################################
if hadd:
  if batch:
    jsub(project='prad',track='analysis',jobname='hadd_{0}{1}_{2}_{3}'.format(suffix,lrun[0],lrun[1]),command='hadd -f {0}/histo{1}_{2}_{3}_{4}.root '.format(outdir,suffix,method,run[0],run[1])+' '.join(['{0}/histo{1}_{2}_{3}.root'.format(outdir,suffix,method,run) for x in get_runs_between(lrun[0],lrun[1],'prod')]),memory='4 GB')
  else: 
    os.system('hadd -f {0}/histo{1}_{2}_{3}_{4}.root '.format(outdir,suffix,method,run[0],run[1])+' '.join(['{0}/histo{1}_{2}_{3}.root'.format(outdir,suffix,method,run) for x in get_runs_between(lrun[0],lrun[1],'prod')]))
  sys.exit()
              

#######################################################################
################ loop on the tree #####################################
#######################################################################
from misc_root import *
from selection import *

print ' '.join(sys.argv)
load_names()

if lrun[0]==-1:
  f = TFile('tree.root')
  t = f.event
  tree_title = t.GetTitle()
  idx = tree_title.index('_')
  method = tree_title[:idx]
  run = int(tree_title[idx+1:])
else:  
  run = lrun[0]
  f = TFile(indir+'/tree_'+method+'_'+str(run)+'.root')
  t = f.event

load_run_params(run)
max_energy = 1500 if run<1362 else 2500
gem = 1 if det=='gem' else 0

fout = TFile('{0}/histo{1}_{2}.root'.format(outdir,suffix,run),'recreate')
hr = {x:listofhisto('hr_'+y,['ep','ee1','ee2'],1000,0.,5.,title=';E_{cluster}/E_{#gamma};') for x,y in module_names.items()}
he = {x:listofhisto('he_'+y,['ep','ee1','ee2'],1000,0.,max_energy,title=';E_{cluster} (MeV);') for x,y in module_names.items()}

t.SetBranchStatus("*",0)
for x in ['Ebeam','n_cl','id','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
n = t.GetEntries()

for (i,x) in enumerate(progress(t,n=n,show=show,precision=2)):
  
  # preselection
  [E,c,idx,Ebeam] = get_variables(t,lincorr=lincorr,mgem=1,match=gem)
  if gem and has_br(c): continue
  if any([t.id[j] in [1561,1562,1595,1596] for j in idx]): continue

  # 2 arm moller events
  if len(idx)==2:
    elas = elasticity(c,E,Ebeam,gem)
    dphi = degrees(deltaphi(c,gem))
    helas_moller2[0].Fill(elas)
    hphidiff[0].Fill(dphi)
    if abs(elas)<0.1 and abs(dphi)<5:
      for j in range(2):
        hr_moller2[t.id[idx[j]]].Fill(elasticity_moller([c[j]],[E[j]],Ebeam,gem)+1)
        he_moller2[t.id[idx[j]]].Fill(E[j])
      helas_moller2[1].Fill(elas)
      hphidiff[1].Fill(dphi)

  # 1 arm moller events
  for j in range(len(idx)):
    elas = elasticity_moller([c[j]],[E[j]],Ebeam,gem)
    helas_moller1[0].Fill(elas)
    if abs(elas)>0.1: continue
    hr_moller1[t.id[idx[j]]].Fill(elas+1)
    he_moller1[t.id[idx[j]]].Fill(E[j])
    helas_moller1[1].Fill(elas)

  # ep events
  for j in range(len(idx)):
    elas = elasticity([c[j]],[E[j]],Ebeam,gem)
    helas_ep[0].Fill(elas)
    if abs(elas)>0.1: continue
    hr_ep[t.id[idx[j]]].Fill(elasticity_ep_el([c[j]],[E[j]],Ebeam,gem)+1)
    he_ep[t.id[idx[j]]].Fill(E[j])
    helas_ep[1].Fill(elas)
 
# writing histograms
fout.cd()
for x in [helas_moller1,helas_moller2,helas_ep,hphidiff]:
  for i in range(2): 
    x[i].Write()

for x,y in module_names.items(): 
  u = fout.mkdir(y)
  u.cd()
  hr_ep[x].Write()
  hr_moller1[x].Write()
  hr_moller2[x].Write()
  he_ep[x].Write()
  he_moller1[x].Write()
  he_moller2[x].Write()
fout.Close()

