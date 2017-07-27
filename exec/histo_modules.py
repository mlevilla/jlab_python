#!/usr/bin/env python

from params import *

largs = [('r',-1,'several','run_number'),
         ('i',workf+'/trees_calib/island','','input_dir'),
         ('o',workf,'','output_folder'),
         ('m','island','','clustering_method'),
         ('n','','','suffix'),
         ('s',0,'','show_bar'),
         ('l',1.0,'','limit'),
         ('b',False,'','batch_mode'),
         ('a',False,'','merge runs'),
         ('nbin',18,'','number_of_energy_bins'),
         ('sigma',-1,'','cut on sigma of elasticity'),
         ('c',False,'','correct mean elasticity')]
print_help(largs)
[lrun,indir,outdir,method,suffix,show,limit,batch,hadd,nbinE,sigma,correct] = catch_args(largs)

######################################################################################
################ batch jobs ##########################################################
######################################################################################
if batch: 
  input_files = ' '.join([indir+'/tree_'+method+'_'+str(run)+'.root' for run in get_runs_between(lrun[0],lrun[-1],'calib')])
  args = del_args(largs,['b','r'])
  jsub(project='prad',track='analysis',jobname='histo{0}_{1}_{2}'.format(suffix,lrun[0],lrun[-1]),command=os.path.abspath('histo_modules.py'),options=args,memory='1024 MB',input_files=input_files,input_data='tree.root',output_data='histo_*'+suffix+'*.root',output_template=outdir+'/@OUTPUT_DATA@',os='centos7',other_files='{0}/mean{1}_2.txt {0}/sigma{1}_2.txt'.format(outdir,suffix))
  sys.exit()

######################################################################################
################ merging histograms ##################################################
######################################################################################
if hadd:
  if batch:
    jsub(project='prad',track='analysis',jobname='hadd'+suffix,command='hadd -f {0}/histo{1}.root {0}/histo_8*{1}.root {0}/histo_9*{1}.root'.format(outdir,suffix),memory='4 GB',os='centos7')
  else:
    os.system('hadd -f {0}/histo{1}.root {0}/histo_8*{1}.root {0}/histo_9*{1}.root'.format(outdir,suffix))
  sys.exit()

######################################################################################
################ loop on the tree ####################################################
######################################################################################
from misc_root import *

print ' '.join(sys.argv)

load_names()
load_positions()

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

load_tagger_timings(run)

lmodule_run = readfile(conff+'/module_run_calib.txt',raw=True)
lmodule_run = [[int(x) for x in y[1:]] for y in lmodule_run if int(y[0])==run][0]
print lmodule_run

fout = TFile('{0}/histo_{1}{2}{3}.root'.format(outdir,run,suffix,('_'+str(int(sigma)))*int(sigma!=-1)),'recreate')

hr = {x:listofhisto('hr_'+module_names[x],[[1,2,5],range(nbinE)],1000,0.,5.,title=';E_{cluster}/E_{#gamma};') for x in lmodule_run}
he = {x:listofhisto('he_'+module_names[x],[[1,2,5]],1200,0.,1200.,title=';E_{cluster} (MeV);') for x in lmodule_run}
heg = {x:listofhisto('heg_'+module_names[x],[[1,2,5]],1200,0.,1200.,title=';E_{#gamma} (MeV);') for x in lmodule_run}

ht = listofhisto('ht',[1,2,5],1000,-100.,100.,title=';t_{hycal}-t_{tagger} (ns);')
hnh = listofhisto('hnh',[1,2,5],20,0,20,title=';n_{hits};')
hid = listofhisto('hid',[1,2,5],2156,0,2156,title=';id_{module};')
hx = listofhisto('hx',[1,2,5],1000,-100,100,title=';x_{cluster}-x_{trans};')
hy = listofhisto('hy',[1,2,5],1000,-100,100,title=';y_{cluster}-y_{trans};')
hr_all = listofhisto('hr',[1,2,5],1000,0.,5.,title=';E_{cluster}/E_{#gamma};')
he_all = listofhisto('he',[1,2,5],1200,0.,1200.,title=';E_{cluster} (MeV);')
heg_all = listofhisto('heg',[1,2,5],1200,0.,1200.,title=';E_{#gamma} (MeV);')

hrx = listofhisto('hrx',[range(8),range(20),range(nbinE)],1000,0.,5.,title=';E_{cluster}/E_{#gamma};')
hegx = listofhisto('hegx',[range(8),range(20),range(nbinE)],1000,0.,5.,title=';E_{#gamma} (MeV);')

tr = {1:0,2:1,5:2}
trans = [353.09,352.75]

if sigma!=-1 or correct:
  mean_peak = readfile(outdir+'/mean{0}_2.txt'.format(suffix),[str]+[float]*nbinE,out=dict)
  sigma_peak = readfile(outdir+'/sigma{0}_2.txt'.format(suffix),[str]+[float]*nbinE,out=dict)
  mean_peak = {x:mean_peak[y] for x,y in module_names.items()}
  sigma_peak = {x:sigma_peak[y] for x,y in module_names.items()}

for _ in progress(t,n=t.GetEntries(),show=show,precision=1,limit=limit):

  if t.n_cl!=1 or t.n_tag!=1: continue
  if t.trigger not in [1,2,5]: continue
  if t.id[0] not in lmodule_run: continue
  j = int((t.E[0]-200.)/(900./nbinE))
  if j>nbinE-1: continue
  if t.Eg[0]==0: continue
  
  if t.E[0]<100: continue
  if t.nh[0]<2: continue
  if abs(t.xhycal[0]-t.xpos)>15: continue
  if abs(t.yhycal[0]-t.ypos)>15: continue
  if abs(t.tg[0]-t.tcl[0]-tagger_timings[t.trigger][0])>min(abs(tagger_timings[t.trigger][i]) for i in range(1,3)): continue

  if sigma!=-1 and abs(t.E[0]/t.Eg[0]-mean_peak[t.id[0]][j])>sigma*sigma_peak[t.id[0]][j]: continue
  
  E = t.E[0] if (not correct or mean_peak[t.id[0]][j]==0) else t.E[0]/mean_peak[t.id[0]][j]
  itr = tr[t.trigger]
  hr[t.id[0]][itr][j].Fill(E/t.Eg[0])
  he[t.id[0]][itr].Fill(E)
  heg[t.id[0]][itr].Fill(t.Eg[0])
  
  ht[itr].Fill(t.tg[0]-t.tcl[0]-tagger_timings[t.trigger][0])
  hnh[itr].Fill(t.nh[0])
  hx[itr].Fill(t.xhycal[0]-t.xpos)
  hy[itr].Fill(t.yhycal[0]-t.ypos)
  hid[itr].Fill(t.id[0])
  hr_all[itr].Fill(E/t.Eg[0])
  he_all[itr].Fill(E)
  heg_all[itr].Fill(t.Eg[0])
                         
  k1,k2 = [],[]
  if abs(t.xhycal[0])<trans[0] and abs(t.yhycal[0]-trans[1])<50.: 
    k1.append(0)
    k2.append(int((t.yhycal[0]-trans[1]+50.)/5.))
  if abs(t.yhycal[0])<trans[1] and abs(t.xhycal[0]+trans[0])<50.:
    k1.append(1)
    k2.append(int(-(t.xhycal[0]+trans[0]-50.)/5.))
  if abs(t.xhycal[0])<trans[0] and abs(t.yhycal[0]+trans[1])<50.:
    k1.append(2)
    k2.append(int(-(t.yhycal[0]+trans[1]-50.)/5.))
  if abs(t.yhycal[0])<trans[1] and abs(t.xhycal[0]-trans[0])<50.: 
    k1.append(3)
    k2.append(int((t.xhycal[0]-trans[0]+50.)/5.))
                         
  if t.id[0]<1000: 
    k1.append(4)
    k2.append(int((t.xhycal[0]-module_pos[t.id[0]][0]*10+48.)/6.))
  if t.id[0]<1000: 
    k1.append(5)
    k2.append(int((t.yhycal[0]-module_pos[t.id[0]][1]*10+48.)/6.))
  if t.id[0]>1000: 
    k1.append(6)
    k2.append(int((t.xhycal[0]-module_pos[t.id[0]][0]*10+24.)/3.))
  if t.id[0]>1000: 
    k1.append(7)
    k2.append(int((t.yhycal[0]-module_pos[t.id[0]][1]*10+24.)/3.))

  for j1,j2 in zip(k1,k2):
    if 0<=j2<20: 
      hrx[j1][j2][j].Fill(E/t.Eg[0])
      hegx[j1][j2][j].Fill(t.Eg[0])
                         

# writing histograms
writelisto([ht,hx,hy,hnh,hid,hr_all,he_all,heg_all],fout,'selection')
writelisto([hrx,hegx],fout,'transition')
for x in hr.keys(): 
  if sum(he[x][i].GetEntries() for i in range(3))>100: 
    writelisto([hr[x],he[x],heg[x]],fout,module_names[x])

fout.Close()
f.Close()

sys.exit()
