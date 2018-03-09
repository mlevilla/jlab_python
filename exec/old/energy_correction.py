#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('i',workf+'/trees/trees_42','','input_folder'),
         ('o','.','','output folder'),
         ('n','','','suffix'),
         ('d','','','density_file'),
         ('c','','','do_energy_correction'),
         ('nbin',200,'','nbin')]

aparser = ArgParser(largs,sub=1)
[lrun,indir,outdir,suffix,density,correct,nbin] = aparser.argl

# batch mode ############################################################################

if batch[0]:
  runs = get_runs_between(lrun[0],lrun[-1],'prod')
  jsub_fast(aparser,'ecorrect',['r','o'],[['r'],runs],outdir,jobname='ecorrect'+suffix,time=wtime[0])
  sys.exit()

from misc_root import *
from moddensity_utils import *
from selection import *

# initialisation ########################################################################

load_positions(run=lrun[0]) 
load_run_calib(lrun[0])
Ebeam = beam_energy[0]/1000.
erange = erange_2GeV if Ebeam>2. else erange_1GeV

nindex = 112
 
cpar = readfile(density,[float]*5,cols=range(1,6)) 
nenergy = len(cpar)/(2*112)
  
def correct_density(t,ie,iregion,ix):
  i = iregion*nenergy+ie+len(cpar)/2*ix
  return t+cpar[i][0]*t*(t**4+cpar[i][1]*t**2+cpar[i][2])*(t**2-0.25)*(t**2-cpar[i][3])
  

f = TFile(indir+'/tree_island_'+str(lrun[0])+'.root')
t = f.event
t.SetBranchStatus('*',0)
for x in ['n_cl','id','xhycal','yhycal','zhycal','nh','E']: t.SetBranchStatus(x,1)

nindex = 353

fout = TFile(outdir+'/ecorrect_'+str(lrun[0])+suffix+'.root','recreate')
hx = listofhisto('hx',[range(nindex),range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hy = listofhisto('hy',[range(nindex),range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxy = listofhisto('hxy',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hE = listofhisto('hE',[range(nindex),range(2)],2500,0,2500,title=';E\' (MeV);')
hxE = listofhisto('hxE',[range(nindex),range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hyE = listofhisto('hyE',[range(nindex),range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxyE = listofhisto('hxyE',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hthE = listofhisto('hthE',range(2),500,0,10,1250,0,2500,title=';#theta (deg);E (MeV)')
if correct:
  hEc = listofhisto('hEc',[range(nindex),range(2)],2500,0,2500,title=';E\' (MeV);')
  hxEc = listofhisto('hxEc',[range(nindex),range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hxEc = listofhisto('hyEc',[range(nindex),range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
  hxyEc = listofhisto('hxyEc',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  hthEc = listofhisto('hthEc',range(2),500,0,10,1250,0,2500,title=';#theta (deg);E\' (MeV)')

# post process ##########################################################################

def custom_exit():
  print
  writelisto([hx,hy],fout,'dist',style=[1,0,0])
  writelisto([hxE,hyE],fout,'Edist',style=[1,0,0])
  writelisto([hxy,hxyE],fout,'occupancy')
  writelisto(hE,fout,'energy',style=[1,0,0])
  writelisto(hthE,fout,'theta')
  if correct:
   writelisto([hxEc,hyEc],fout,'Edist_correct',style=[1,0,0])
   writelisto([hxyc,hxyEc],fout,'occupancy_correct')
   writelisto(hEc,fout,'energy_correct',style=[1,0,0])
   writelisto(hthEc,fout,'theta_correct')

  fout.Close()
  f.Close()

def signal_handler(signal,frame):
  custom_exit()
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# loop ##################################################################################

for _ in progress(t,show=show_progress[0],n=t.GetEntries(),precision=1):
  for k in range(t.n_cl):
    if t.nh[k]<2: continue
    cid0,x,y = t.id[k], t.xhycal[k], t.yhycal[k]
    if cid0 in edge+[1561,1562,1595,1596]: continue
    cid1,tx,ty = correct_module_id(cid0,x,y)
    if cid1 in edge+deads: continue
    idx_density = module2index(cid1)
    idx_module = module2index_asym(cid1)
    iregion = regions_id(idx_density)
    theta = ftheta([x,y,t.zhycal[k]])
    E = correct_linearity(t.E[k],cid1)
    E_theo = [ep_energy_el(theta,Ebeam),moller_energy(theta,Ebeam)]
    elas = [E/etheo/1000.-1 for etheo in E_theo] 
    elasc = [4*[0.024,0.062][cid1<1000]/etheo**0.5 for etheo in E_theo]
    
    for ie in range(2):
      if ie==1 and t.n_cl!=1: continue
      if abs(elas[ie])>elasc[ie] : continue
      new_tx = correct_density(tx,ie+3,iregion,0)
      new_ty = correct_density(ty,ie+3,iregion,1)
      _ = hx[idx_module][ie].Fill(tx)
      _ = hy[idx_module][ie].Fill(ty)
      _ = hxy[ie].Fill(x,y)
      _ = hE[idx_module][ie].Fill(E)
      _ = hxE[idx_module][ie].Fill(tx,E/E_theo[ie])
      _ = hyE[idx_module][ie].Fill(ty,E/E_theo[ie])
      _ = hxyE[ie].Fill(x,y,E/E_theo[ie])
      _ = hxyE[ie].Fill(x,y,E/E_theo[ie])
      _ = hthE[ie].Fill(theta/degrad,E)
        

custom_exit()
