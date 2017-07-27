#!/apps/anaconda/anaconda-2.0.1/bin/python

from params import *

run = ''
output = work_folder+'/calib'
indir = replay_folder_me
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
  elif a in ['-i','-in']: indir = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-b','-batch']: batch = True
  elif a in ['-s','-show']: show = int(args[i+1])
  elif a in ['-f','-final']: final = args[i+1]
  elif a in ['-g','-gains']: gains_file = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-m','-method']: method = args[i+1]

###################
## uniting gains ##
###################
if final!=None:
  from misc_root import *

  # raw gains
  runs = sum_dict_lists(runs_list,'calib')
  l = [readfile(output+'_'+str(run)+'.txt',[str,float,int]) for run in runs]
  l2 = [[l[0][i][0],sum([x[i][1]*x[i][2] for x in l])/sum([x[i][2] for x in l]),sum([x[i][2] for x in l])] if sum([x[i][2] for x in l])!=0 else [l[0][i][0],0.25,0] for i in range(len(l[0]))]
  writelist(os.path.dirname(output)+'/gains_calib'+final+'_raw.txt',l2)
  sys.stdout.write('raw gains computed\n')

  # gaussian gains
  if gains_file == '': gains0 = [0.25 for i in range(1728)]
  else: gains0 = readfile(gains_file,[float],cols=[1])
  f = [TFile(output+'_'+str(run)+'.root') for run in runs]
  h = [TH1F('hratio_'+x+'_calib',';E_{cluster}/E_{#gamma};',500,0,2.5) for x in module_names]
  l2 = []
  for j,x in enumerate(f):
    sys.stdout.write(str(runs[j])+'\r')
    sys.stdout.flush()
    for i,y in enumerate(module_names):
      h[i].Add(x.Get(y).Get('hratio_'+y))
  for i in range(len(module_names)):
    h[i].Fit('gaus','qww','',0.5,2.5)
    m = h[i].GetFunction('gaus').GetParameter(1)
    g = gains0[i]/m if m!=0. else gains0[i]
    s = h[i].GetFunction('gaus').GetParameter(2)/m if m!=0. else h[i].GetFunction('gaus').GetParameter(2)
    l2.append([module_names[i],g,s])
  writelist(os.path.dirname(output)+'/gains_calib'+final+'_gaus.txt',l2)
  sys.stdout.write('gaussian gains computed\n')

  # merging trees
  c = TChain('hycal','hycal')
  for run in runs: c.Add(output+'_'+str(run)+'.root')
  c.Merge(os.path.dirname(output)+'/tree_calib'+final+'.root')
  sys.stdout.write('merging trees finished\n')
  sys.exit()


if run=='': sys.exit('no run mentionned')


################
## batch mode ##
################
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+'_'+run,command=os.getcwd()+'/calibrate_calib.py',options='-o '+os.path.basename(output)+' -i . -r '+run+' -g '+gains_file+' -s '+str(show)+' -m '+method,input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run+'*',output_template=os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='10 GB',memory='1024 MB')
  print
  sys.exit()
  

##########
## Main ##
##########

# initialisation
from readdst import *
from tagger import *
from hycal import *
from misc_root import *

t = chao(run=int(run),folder=indir)
t.goto(20000)
t.lload = ['','adc','tdc']

if gains_file == '': gains0 = [0.25 for i in range(1728)]
else: gains0 = readfile(gains_file,[float],cols=[1])

ped0 = get_pedestal(t.run)
sigma_ped0 = get_pedestal_sigma(t.run)
lms_gains0 = get_lms_gains(t.run,889,979)

fout = TFile(output+'_'+run+'.root','recreate')

henergy = [TH1F('henergy_'+module_names[i],'',500,0.,1500) for i in range(1728)]
hratio = [TH1F('hratio_'+module_names[i],';E_{cluter}/E_{#gamma};',500,0.,2.5) for i in range(1728)]
hx = [TH1F('hx_'+module_names[i],'',400,-20.,20.) for i in range(1728)]
hy = [TH1F('hy_'+module_names[i],'',400,-20.,20.) for i in range(1728)]
htrigger = [[TH1F('htrigger_'+module_names[i]+'_'+str(j),'',500,0.,1500) for i in range(1728)] for j in range(6)]
l = [('iev','i',''),('trigger','i',''),('E_g','f',''),('xpos','f',''),('ypos','f',''),('n_cl','i',''),('mid','i','[n_cl]'),('E_cl','f','[n_cl]'),('x_cl','f','[n_cl]'),('y_cl','f','[n_cl]'),('chi2_cl','f','[n_cl]'),('nhit_cl','i','[n_cl]'),('nleak_cl','i','[n_cl]'),('nneigh_cl','i','[n_cl]'),('status_cl','i','[n_cl]')]
tree,var = tree_init('hycal',l,nmax=100)
lgain = [0. for i in range(1728)]
ln = [0 for i in range(1728)]

# loop on events
for ev in progress(t,modulo=10000,show=show):
  # stability and trigger cuts
  if t.isepics: continue
  if event_is_bad(ev.iev,t.run): continue
  if ev.trigger not in [1,2,5]: continue

   # tagger energy
  etc = TDCHits([x for x in get_etchannels(ev=ev) if abs(x.t)<25])
  if etc==[]: continue
  etc = etc.merge()
  eg = t.epics['MBSY2C_energy']*etc[0].E
  # hycal event
  hev = hcevent(ev,2,gains=gains0,peds=ped0,sigma_peds=sigma_ped0,lms_gains=lms_gains0,method=method)
  if hev.clusters==[]: continue

  # tree filling
  xtrans = (t.epics['hallb_ptrans_x_encoder']+652.05)/10.
  ytrans = (-(t.epics['hallb_ptrans_y1_encoder']+t.epics['hallb_ptrans_y2_encoder'])/2.-3761.5)/10.
  var['iev'][0] = ev.iev
  var['trigger'][0] = ev.trigger
  var['E_g'][0] = eg
  var['xpos'][0] = xtrans
  var['ypos'][0] = ytrans
  n_cl = 0

  # clusters
  for cl in hev.clusters:
    if cl.E/eg<0.5: continue
    htrigger[ev.trigger][cl.mid].Fill(eg)

    if cl.E<100 or cl.nhit<3: continue
    lgain[cl.mid]+=cl.E/eg
    ln[cl.mid]+=1
    # tree
    var['mid'][n_cl] = cl.mid
    var['E_cl'][n_cl] = cl.E
    var['x_cl'][n_cl] = cl.x
    var['y_cl'][n_cl] = cl.y
    var['chi2_cl'][n_cl] = cl.chi2
    var['nhit_cl'][n_cl] = len(cl.hits)
    var['nleak_cl'][n_cl] = len(cl.leaks)
    var['nneigh_cl'][n_cl] = len(cl.leaks)
    var['status_cl'][n_cl] = cl.status
    n_cl+=1
    # histos
    if ev.trigger not in [1,2]: continue
    henergy[cl.mid].Fill(cl.E,cl.E/eg)
    hratio[cl.mid].Fill(cl.E/eg)
    hx[cl.mid].Fill(cl.x-xtrans)
    hy[cl.mid].Fill(cl.y-ytrans)
    
  var['n_cl'][0] = n_cl
  tree.Fill()

# writing gains
gains1 = [z/(x/y) if y!=0 else z for x,y,z in zip(lgain,ln,gains0)] 
writelist(output+'_'+run+'.txt',[module_names,gains1,ln],tr=True)

# writing trees and histograms
fout.cd()
tree.Write()
for i in range(1728):
  folder = fout.mkdir(module_names[i])
  folder.cd()
  henergy[i].Write()
  hratio[i].Write()
  hx[i].Write()
  hy[i].Write()
  for j in range(6): htrigger[j][i].Write()
fout.Close()

sys.exit(0)
