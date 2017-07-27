#!/apps/anaconda/anaconda-2.0.1/bin/python

from params import *

run = ''
output = work_folder+'/calib'
indir = work_folder+'/calib_file'
batch = False
show = 0
gains_file = ''
final = None
method = 'island'

# input parameters
args = sys.argv
#if len(args)<3: sys.exit('no arguments')
for i,a in enumerate(args):
  if a in ['-r','-run']: run = args[i+1]
  elif a in ['-o','-out']: output = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-i','-in']: indir = args[i+1]
  elif a in ['-b','-batch']: batch = True
  elif a in ['-s','-show']: show = int(args[i+1])
  elif a in ['-f','-final']: final = args[i+1]
  elif a in ['-g','-gains']: gains_file = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-m','-method']: method = args[i+1]

# uniting gains
if final!=None:
  # raw gains
  runs = periods_prod
  l0 = [[readfile(output+'_'+str(run)+'.txt',[str,float,int]) for run in x] for x in runs]
  l2 = [[[l[0][i][0],sum([x[i][1]*x[i][2] for x in l])/sum([x[i][2] for x in l]),sum([x[i][2] for x in l])] if sum([x[i][2] for x in l])!=0 else [l[0][i][0],0.25,0] for i in range(len(l[0]))] for l in l0]
  for i,x in enumerate(l2): writelist(os.path.dirname(output)+'/gains_prod'+final+'_pr'+str(i)+'.txt',x)

  # gaussian gains
  if gains_file == '': gains0 = [0.25 for i in range(1728)]
  else: gains0 = readfile(gains_file,[float],cols=[1])
  f = [[TFile(output+'_'+str(run)+'.root') for run in x] for x in runs]
  h = [[[TH1F('hratio_'+x+'_calib_'+str(i)+'_'+str(j),';E_{cluster}/E_{#gamma};',500,0,2.5) for x in module_names] for i in range(len(runs))] for j in range(2)]
  l2 = [[[] for i in range(len(runs))] for j in range(2)]
  for i,y in enumerate(f):
    for j,x in enumerate(y):
      sys.stdout.write(str(runs[i][j])+'\r')
      sys.stdout.flush()
      for k,z in enumerate(module_names):
        h[0][i][k].Add(x.Get(z).Get('hratio_'+y+'_0'))
        h[1][i][k].Add(x.Get(z).Get('hratio_'+y+'_1'))

  for i in range(len(runs)):
    for k,z in enumerate(module_names):
      h[0][i][k].Fit('gaus','qww','',0.5,2.5)
      h[1][i][k].Fit('gaus','qww','',0.5,2.5)
      m0 = h[0][i][k].GetFunction('gaus').GetParameter(1)
      m1 = h[1][i][k].GetFunction('gaus').GetParameter(1)
      g0 = gains0[k]/m0 if m0!=0. else gains0[k]
      g1 = gains0[k]/m1 if m1!=0. else gains0[k]
      s0 = h[0][i][k].GetFunction('gaus').GetParameter(2)/m0 if m0!=0. else h[0][i][k].GetFunction('gaus').GetParameter(2)
      s1 = h[1][i][k].GetFunction('gaus').GetParameter(2)/m1 if m0!=0. else h[1][i][k].GetFunction('gaus').GetParameter(2)
      l2[0][i].append([z,g0,s0])
      l2[1][i].append([z,g1,s1])
  for i in range(len(runs)):
    writelist(os.path.dirname(output)+'/gains_calib'+final+'_gaus_'+str(i)+'_moller.txt',l2[0][i])
    writelist(os.path.dirname(output)+'/gains_calib'+final+'_gaus_'+str(i)+'_ep.txt',l2[0][i])
    
  #merging trees
  c = TChain('hycal','hycal')
  for x in runs:
    for run in x:
      c.Add(output+'_'+str(run)+'.root')
  c.Merge(os.path.dirname(output)+'/tree_calib'+final+'.root')    
  sys.stdout.write('merging trees finished\n')
  sys.exit()
  

if run=='': sys.exit('no run mentionned')

# batch mode
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+'_'+run,command=os.getcwd()+'/calibrate_prod.py',options='-o '+os.path.basename(output)+' -i . -r '+run+' -g '+gains_file+' -s '+str(show)+' -m '+method,input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run+'*',output_template=os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='10 GB',memory='1024 MB',time=5000)
  print
  sys.exit()
  

# initialisation
from readdst import *
from hycal import *
from misc_root import *

if gains_file == '': gains0 = [0.25 for i in range(1728)]
else: gains0 = readfile(gains_file,[float],cols=[1])

t = dst(run=int(run),fmt=fmt_calib_prod,lload=['','adc'],folder=indir)

ped0 = get_pedestal(t.run)
sigma_ped0 = get_pedestal_sigma(t.run)
lms_gains0 = get_lms_gains(t.run,1122,1515)

fout = TFile(output+'_'+run+'.root','recreate')

hratio = [[TH1F('hratio_'+module_names[i]+'_'+str(j),'',500,0.,2.5) for i in range(1728)] for j in range(2)]
l = [('iev','i',''),('trigger','i',''),('ebeam','f',''),('n_cl','i',''),('mid','i','[n_cl]'),('E_cl','f','[n_cl]'),('x_cl','f','[n_cl]'),('y_cl','f','[n_cl]'),('chi2_cl','f','[n_cl]'),('ep_cl','i','[n_cl]'),('nhit_cl','i','[n_cl]'),('nleak_cl','i','[n_cl]'),('nneigh_cl','i','[n_cl]')]
tree,var = tree_init('hycal',l,nmax=100)
lgain = [[0. for i in range(1728)] for j in range(2)]
ln = [[0 for i in range(1728)] for j in range(2)]

# loop on events
for ev in progress(t,modulo=10000,show=show):

  if event_is_bad(ev.iev,t.run): continue
  if ev.trigger not in [1,2]: continue


  # hycal event
  hev = hcevent(ev,1,gains=gains0,peds=ped0,sigma_peds=sigma_ped0,lms_gains=lms_gains0,method=method)
  hev.clusters = [c for c in hev.clusters if (c.E>50 and c.nhit>2)]
  if hev.clusters==[] or len(hev.clusters)>4: continue
  moller_cls = hev.moller_cls(r=0.5,ebeam=ev.Eb)
  ep_cls = hev.ep_cls(r=0.5,ebeam=ev.Eb)

  # tree filling 
  var['iev'][0] = ev.iev
  var['trigger'][0] = ev.trigger
  var['ebeam'][0] = ev.Eb
  n_cl = 0

  # moller gain filling
  for p in moller_cls:
    d,d1,d2,x,y = distance_center_cls(p[0],p[1])
    lgain[0][p[0].mid]+=p[0].E/(ev.Eb*d2/(d1+d2))
    lgain[0][p[1].mid]+=p[1].E/(ev.Eb*d1/(d1+d2))
    ln[0][p[0].mid]+=1
    ln[0][p[1].mid]+=1
    hratio[0][p[0].mid].Fill(p[0].E/(ev.Eb*d2/(d1+d2)))
    hratio[0][p[1].mid].Fill(p[1].E/(ev.Eb*d1/(d1+d2)))
    # tree
    var['mid'][n_cl] = p[0].mid
    var['E_cl'][n_cl] = p[0].E
    var['x_cl'][n_cl] = p[0].x
    var['y_cl'][n_cl] = p[0].y
    var['chi2_cl'][n_cl] = p[0].chi2
    var['nhit_cl'][n_cl] = len(p[0].hits)
    var['nleak_cl'][n_cl] = len(p[0].leaks)
    var['nneigh_cl'][n_cl] = len(p[0].neighbors)
    var['ep_cl'][n_cl] = 0
    n_cl+=1
    var['mid'][n_cl] = p[1].mid
    var['E_cl'][n_cl] = p[1].E
    var['x_cl'][n_cl] = p[1].x
    var['y_cl'][n_cl] = p[1].y
    var['chi2_cl'][n_cl] = p[1].chi2
    var['nhit_cl'][n_cl] = len(p[1].hits)
    var['nleak_cl'][n_cl] = len(p[1].leaks)
    var['nneigh_cl'][n_cl] = len(p[1].neighbors)
    var['ep_cl'][n_cl] = 0
    n_cl+=1

  # ep gain filling
  for c in ep_cls:
    lgain[1][c.mid]+=c.E/ev.Eb
    ln[1][c.mid]+=1
    hratio[1][c.mid].Fill(c.E/ev.Eb)
    # tree 
    var['mid'][n_cl] = c.mid
    var['E_cl'][n_cl] = c.E
    var['x_cl'][n_cl] = c.x
    var['y_cl'][n_cl] = c.y
    var['chi2_cl'][n_cl] = c.chi2
    var['nhit_cl'][n_cl] = len(c.hits)
    var['nleak_cl'][n_cl] = len(c.leaks)
    var['nneigh_cl'][n_cl] = len(c.neighbors)
    var['ep_cl'][n_cl] = 1
    n_cl+=1
  
  var['n_cl'][0] = n_cl
  tree.Fill()

# writing gains
gains1 = [[z/(x/y) if y!=0 else z for x,y,z in zip(lgain[i],ln[i],gains0)] for i in range(2)] 
writelist(output+'_'+run+'_moller.txt',[module_names,gains1[0],ln[0]],tr=True)
writelist(output+'_'+run+'_ep.txt',[module_names,gains1[1],ln[1]],tr=True)

# writing trees and histograms
fout.cd()
tree.Write()
for i in range(1728): 
  folder = fout.mkdir(module_names[i])
  folder.cd()
  for j in range(2): hratio[j][i].Write()
fout.Close()

sys.exit(0)
