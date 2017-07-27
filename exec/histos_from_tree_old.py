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
         ('a',False,'','merge runs')]
print_help(largs)
[run,indir,outdir,method,suffix,show,limit,batch,hadd] = [catch_arg(x,y,z) for x,y,z,c in largs]

if ipass=='' and indir=='': sys.exit('not enough arguments: '+sys.argv)
if ipass=='' and '/pass' in indir: 
  ind = indir.index('/pass')
  ipass = indir[ind+5:]
if indir=='': indir = os.environ['WORK']+'/calib_results/'+method+'/pass'+ipass

######################################################################################
################ batch jobs ##########################################################
######################################################################################
if batch: 
  jsub(project='prad',track='analysis',jobname='histo'+name+'_'+run,command=os.path.abspath('histos_from_tree.py'),options=' '.join(sys.argv[1:]).replace('-b',''),memory='2 GB')
  sys.exit()

######################################################################################
################ merging histograms ##################################################
######################################################################################
if hadd:
  jsub(project='prad',track='analysis',jobname='hadd'+name,command='hadd -f '+indir+'/d'+name+'/histo'+name+'.root '+indir+'/d'+name+'/histo'+name+'_8*.root '+indir+'/d'+name+'/histo'+name+'_9*.root',memory='4 GB')
  sys.exit()

######################################################################################
################ loop on the tree ####################################################
######################################################################################
from ROOT import TFile, TH1F

print ' '.join(sys.argv)

fout = TFile(indir+'/d'+name+'/'+'histo'+name+'_'+run+'.root','recreate')
hr = [[TH1F('hr_'+x+'_'+str(j),';E_{cluster}/E_{#gamma};',1000,0.,5.) for j in range(9)] for x in module_names]
he = [[TH1F('he_'+x+'_'+str(j),';E_{cluster} (MeV);',1200,0.,1200.)  for j in [1,2,5]] for x in module_names]
heg = [[TH1F('heg_'+x+'_'+str(j),';E_{#gamma} (MeV);',1200,0.,1200.)  for j in [1,2,5]] for x in module_names]
ht = [[TH1F('htdiff_'+x+'_'+str(j),';t_{hycal}-t_{tagger} (ns);',1000,-200.,200.)  for j in [1,2,5]] for x in module_names]
hrx = [[[TH1F('hrx_'+str(i)+'_'+str(j)+'_'+str(k),';E_{cluster}/E_{#gamma};',1000,0.,5.) for k in range(9)] for j in range(20)] for i in range(8)]
hegx = [[[TH1F('hegx_'+str(i)+'_'+str(j)+'_'+str(k),';E_{cluster}/E_{#gamma};',1200,0.,1200.) for k in range(9)] for j in range(20)] for i in range(8)]
tr = {1:0,2:1,5:2}
trans = [353.09,352.75]
f = TFile(indir+'/trees/tree_'+run+'.root')
t = f.hycal
n = t.GetEntries()

if 's' in var:
  base_name = name[:-3]
  mean_peak = readfile(indir+'/d_t40_xy/mean_t40_xy_s2.txt',[float]*9,cols=range(1,10))
  sigma_peak = readfile(indir+'/d_t40_xy/sigma_t40_xy_s2.txt',[float]*9,cols=range(1,10))

if 'x' in var or 'y' in var:
  trans_offset = readfile(os.environ['HOME']+'/python/transporter_offset.txt',[str,int,int,int,int],out=dict)

for (i,x) in enumerate(progress(t,n=n,show=show,modulo=10000)):

  if t.n_cl!=1: continue
  if t.trigger not in [1,2,5]: continue
  j = int((t.E[0]-200)/100.)
  if j>8: continue
  if t.Eg==0: continue
  b = False

  for v,c in zip(var,cut):
    if v=='e': b = b or (t.E[0]<c)
    if v=='n': b = b or t.nh[0]<c
    if v=='chi2': b = b or t.chi2[0]>c
    if v=='x': b = b or abs(t.xpos-t.x[0]+trans_offset[run][0])>c*trans_offset[run][2]
    if v=='y': b = b or abs(t.ypos-t.y[0]+trans_offset[run][1])>c*trans_offset[run][3]
    if v=='t': b = b or abs(t.tg-t.tcl[0])>c
    if v=='tr': b = b or t.trigger>c
    if v=='s': b = b or abs(t.E[0]/t.Eg-mean_peak[t.id[0]][j])>c*sigma_peak[t.id[0]][j]
  if b: continue

  he[t.id[0]][tr[t.trigger]].Fill(t.E[0])
  heg[t.id[0]][tr[t.trigger]].Fill(t.Eg)
  ht[t.id[0]][tr[t.trigger]].Fill(t.tg-t.tcl[0])
  hr[t.id[0]][j].Fill(t.E[0]/t.Eg)

  k1,k2 = [],[]
  if abs(t.x[0])<trans[0] and abs(t.y[0]-trans[1])<50.: 
    k1.append(0)
    k2.append(int((t.y[0]-trans[1]+50.)/5.))
  if abs(t.y[0])<trans[1] and abs(t.x[0]+trans[0])<50.:
    k1.append(1)
    k2.append(int(-(t.x[0]+trans[0]-50.)/5.))
  if abs(t.x[0])<trans[0] and abs(t.y[0]+trans[1])<50.:
    k1.append(2)
    k2.append(int(-(t.y[0]+trans[1]-50.)/5.))
  if abs(t.y[0])<trans[1] and abs(t.x[0]-trans[0])<50.: 
    k1.append(3)
    k2.append(int((t.x[0]-trans[0]+50.)/5.))

  if t.id2[0]<1000: 
    k1.append(4)
    k2.append(int((t.x[0]-xpos[t.id[0]]*10+48.)/6.))
  if t.id2[0]<1000: 
    k1.append(5)
    k2.append(int((t.y[0]-ypos[t.id[0]]*10+48.)/6.))
  if t.id2[0]>1000: 
    k1.append(6)
    k2.append(int((t.x[0]-xpos[t.id[0]]*10+24.)/3.))
  if t.id2[0]>1000: 
    k1.append(7)
    k2.append(int((t.y[0]-ypos[t.id[0]]*10+24.)/3.))

  for j1,j2 in zip(k1,k2):
    if 0<=j2<20: 
      hrx[j1][j2][j].Fill(t.E[0]/t.Eg)
      hegx[j1][j2][j].Fill(t.Eg)


# writing histograms
fout.cd()
for i,x in enumerate(module_names): 
  u = fout.mkdir(x)
  u.cd()
  for j in range(9): hr[i][j].Write()
  for j in range(3): 
    he[i][j].Write()
    heg[i][j].Write()
    ht[i][j].Write()
for i in range(8):
  u = fout.mkdir('trans_'+['top','right','bottom','left','x_lg','y_lg','x_pwo','y_pwo'][i])
  u.cd()
  for j in range(20):
    v = u.mkdir('x_'+str(j))
    v.cd()
    for k in range(9): 
      hrx[i][j][k].Write()
      hegx[i][j][k].Write()
fout.Close()
