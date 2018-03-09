#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('n','','','suffix'),
         ('o','.','','output_folder'),
         ('l',1.0,'','limit'),
         ('i',workf+'/trees/trees_42','','input_dir'),
         ('gem',1,'','use_gem_coordinate'),
         ('rdead',1.,'','radius_for_dead_module'), # cut
         ('fiducial',True,'','fiducial_cut'), # cut
         ('spacer',False,'','gem_spacer_and_dead_area'), # cut
       ]
aparser = ArgParser(largs,sub=1)
[lrun,suffix,outdir,limit,indir,gem,rdead,fiducial,spacer] = aparser.argl


# batch processing ######################################################################

if batch[0]:
  runs = get_runs_between(lrun[0],lrun[-1]) if len(lrun)<=2 else lrun
  jsub_fast(aparser,'eresolution',['r','o'],[['r'],runs],outdir,jobname='eresolution'+suffix,time=wtime[0])
  sys.exit()

print ' '.join(sys.argv)
date = subprocess.check_output(['date'])
print 'start: ',date[:-1]

# initialisation ########################################################################

from misc_root import *
from selection import *
from moddensity_utils import *

run = lrun[0]
f = TFile(indir+'/tree_island_{0}.root'.format(run))
t = f.event
t.SetBranchStatus("*",0)
for x in ['iev','Ebeam','n_cl','id','nh','E','xhycal','yhycal','zhycal','xgem','ygem','zgem']: t.SetBranchStatus(x,1)
fout = TFile(outdir+'/eresolution_{0}{1}.root'.format(run,suffix),'recreate')

load_run_params(run)
load_run_calib(run)
load_positions(run=run)
load_dead_pos(run=run)

Ebeam = beam_energy[0]/1000.
hth = listofhisto('hth',range(150),1000,-1,1,title=';E\'/E_{th}-1;')
htx = listofhisto('htx',[],1000,-1,1,title=';;')
htxe = listofhisto('htxe',[],1000,-1,1,title=';;E\'/E_{th}-1')

# post process ##########################################################################
def custom_exit():
  writelisto([hth,htx,htxe],fout)

def signal_handler(signal,frame):
  custom_exit()
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)



# loop ##################################################################################

for _ in progress(enumerate(t),show=show_progress[0],n=t.GetEntries(),precision=1,limit=limit):
  
  # variables
  [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=fiducial,exclude_dead=0,match=0,rdead=rdead,spacer=spacer)

  theta = [[ftheta(y) if y!=[] else 0 for y in x] for x in c]
  thetadeg = [[abs(degrees(y)) for y in x] for x in theta]

  # ep
  for j in range(len(idx)):
    th = thetadeg[j][gem]
    if c[j][gem]==[]: continue
    if th<0.7: continue
    if gem_fiducial(c[j][gem],rdead): continue
    E_ep = [ep_energy_el(theta[j][i],Ebeam) for i in range(2)]
    elas_ep = E[j]/E_ep[gem]/1000-1
    elasc = [0.024,0.062][lg[j]]/E_ep[gem]**0.5
    #if abs(elas_ep)>10*elasc: continue

    for k in range(120):
      if th<0.5+k*0.025 or th>0.5+(k+1)*0.025: continue
      hth[k].Fill(elas_ep)

    for k in range(30):
      if th<3.5+k*0.1 or th>3.5+(k+1)*0.1: continue
      hth[k+120].Fill(elas_ep)
      
    if abs(elas_ep)>6*elasc: continue
    cid1,tx,ty = correct_module_id(cid[j],c[j][0][0],c[j][0][1])
    htx.Fill(tx)
    htxe.Fill(tx,abs(elas_ep))
      
 
# exit ##################################################################################

custom_exit()

date = subprocess.check_output(['date'])
print 'end: ',date[:-1]
