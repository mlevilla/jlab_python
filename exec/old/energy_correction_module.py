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
  jsub_fast(aparser,'ecorrect_mod',['r','o'],[['r'],runs],outdir,jobname='ecorrect_mod'+suffix,time=wtime[0])
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
  
if correct: 
  cpar_ene = readfile(correct,[float]*6,cols=range(1,7))
  lpar_ene = len(cpar_ene)/2

def correct_energy(t,ie,iregion,ix):
  i = lpar_ene*ie+iregion
  return 1/(cpar_ene[i][ix*3]+cpar_ene[i][ix*3+1]*t**2+cpar_ene[i][ix*3+2]*t**4)
  
  

f = TFile(indir+'/tree_island_'+str(lrun[0])+'.root')
t = f.event
t.SetBranchStatus('*',0)
for x in ['iev','Ebeam','n_cl','id','xhycal','yhycal','zhycal','xgem','ygem','zgem','nh','E']: t.SetBranchStatus(x,1)


fout = TFile(outdir+'/ecorrect_mod_'+str(lrun[0])+suffix+'.root','recreate')
hx = listofhisto('hx',[module_names,range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hy = listofhisto('hy',[module_names,range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxy = listofhisto('hxy',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
helas = listofhisto('helas',[module_names,range(2)],1000,-1,1,title=';E\' (MeV);')
hE = listofhisto('hE',[module_names,range(2)],2500,0,2500,title=';E\' (MeV);')
hxE = listofhisto('hxE',[module_names,range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hyE = listofhisto('hyE',[module_names,range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxyE = listofhisto('hxyE',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hthE = listofhisto('hthE',range(2),500,0,10,1250,0,2500,title=';#theta (deg);E (MeV)')
if correct:
  helasc = listofhisto('helasc',[module_names,range(2)],1000,-1,1,title=';E\' (MeV);')
  hEc = listofhisto('hEc',[module_names,range(2)],2500,0,2500,title=';E\' (MeV);')
  hxEc = listofhisto('hxEc',[module_names,range(2)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hyEc = listofhisto('hyEc',[module_names,range(2)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
  hxyEc = listofhisto('hxyEc',range(2),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  hthEc = listofhisto('hthEc',range(2),500,0,10,1250,0,2500,title=';#theta (deg);E\' (MeV)')

# post process ##########################################################################

def custom_exit():
  print
  writelisto([hx,hy],fout,'dist',style=[1,0,0])
  writelisto([hxE,hyE],fout,'Edist',style=[1,0,0])
  writelisto([hxy,hxyE],fout,'occupancy')
  writelisto(hE,fout,'energy',style=[1,0,0])
  writelisto(helas,fout,'elasticity',style=[1,0,0])
  writelisto(hthE,fout,'theta')
  if correct:
   writelisto([hxEc,hyEc],fout,'Edist_correct',style=[1,0,0])
   writelisto([hxyEc],fout,'occupancy_correct')
   writelisto(hEc,fout,'energy_correct',style=[1,0,0])
   writelisto(helasc,fout,'elasticity_correct',style=[1,0,0])
   writelisto(hthEc,fout,'theta_correct')

  fout.Close()
  f.Close()

def signal_handler(signal,frame):
  custom_exit()
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# loop ##################################################################################

for _ in progress(t,show=show_progress[0],n=t.GetEntries(),precision=1):

  [E,c,idx,lg,cid,_] = get_variables(t,lincorr=1,mgem=1,exclude_edge=1,exclude_dead=0,match=0,rdead=0,spacer=0,eloss=1)

  theta = [[ftheta(y) if y!=[] else 0 for y in x] for x in c]
  thetadeg = [[abs(degrees(y)) for y in x] for x in theta]
  tmp = [correct_module_id(cid[j],c[j][0][0],c[j][0][1]) for j in range(len(idx))]
  cid1 = [x[0] for x in tmp]
  tx = [x[1] for x in tmp]
  ty = [x[2] for x in tmp]


  # ep
  for j in range(len(idx)):
    if len(idx)!=1: break
    if c[0][1]==[]: continue
    E_ep = ep_energy_el(theta[0][1],Ebeam)
    elas_ep = E[0]/E_ep/1000-1
    elasc = 4*[0.024,0.062][lg[0]]/E_ep**0.5
    if cid1[0] in edge+dead_modules: continue
    if abs(elas_ep)>elasc : continue
    idx_density = module2index(cid1[0])
    new_tx = correct_density(tx[0],3,idx_density,0)
    new_ty = correct_density(ty[0],3,idx_density,1)
    _ = hx[cid1[0]][0].Fill(new_tx)
    _ = hy[cid1[0]][0].Fill(new_ty)
    new_x = new_tx*cell_size[int(cid1[0]<1000)][0]+module_pos[cid1[0]][0]
    new_y = new_ty*cell_size[int(cid1[0]<1000)][1]+module_pos[cid1[0]][1]
    _ = hxy[0].Fill(new_x,new_y)
    _ = hE[cid1[0]][0].Fill(E[0])
    _ = helas[cid1[0]][0].Fill(elas_ep)
    _ = hxE[cid1[0]][0].Fill(new_tx,elas_ep+1)
    _ = hyE[cid1[0]][0].Fill(new_ty,elas_ep+1)
    _ = hxyE[0].Fill(new_x,new_y,elas_ep+1)
    _ = hxyE[0].Fill(new_x,new_y,elas_ep+1)
    _ = hthE[0].Fill(thetadeg[0][1],E[0])
    if correct:
      E[0] = E[0]*correct_energy(new_tx,0,cid1[0],0)*correct_energy(new_ty,0,cid1[0],1)
      elas_ep = E[0]/E_ep/1000-1
      _ = hEc[cid1[0]][0].Fill(E[0])
      _ = helasc[cid1[0]][0].Fill(elas_ep)
      _ = hxEc[cid1[0]][0].Fill(new_tx,elas_ep+1)
      _ = hyEc[cid1[0]][0].Fill(new_ty,elas_ep+1)
      _ = hxyEc[0].Fill(new_x,new_y,elas_ep+1)
      _ = hxyEc[0].Fill(new_x,new_y,elas_ep+1)
      _ = hthEc[0].Fill(thetadeg[0][1],E[0])

  Em = [[moller_energy(y,Ebeam) if y!=0 else 0 for y in x] for x in theta]
  if len(idx)!=2: continue
  if c[0][1]==[] or c[1][1]==[]: continue
  if cid1[0] in edge+dead_modules: continue
  if cid1[1] in edge+dead_modules: continue
  dphi = degrees((fphi(c[0][1])-fphi(c[1][1]))%(2*pi)-pi)
  dE = E[0]+E[1]-Ebeam*1e3
  elas_sum = dE/Ebeam/1e3
  dEcut = 4*(([0.024,0.062][lg[0]]*E[0]**0.5)**2+([0.024,0.062][lg[1]]*E[1]**0.5)**2)**0.5*1000**0.5
  if abs(dphi)>5: continue
  if abs(dE)>dEcut: continue
  for j in range(2):
    elas_m = E[j]/Em[j][1]/1000.-1
    elasc = 4*[0.024,0.062][lg[j]]/(Em[j][1])**0.5
    if abs(elas_m)>elasc: continue
    idx_density = module2index(cid1[j])
    new_tx = correct_density(tx[j],4,idx_density,0)
    new_ty = correct_density(ty[j],4,idx_density,1)
    _ = hx[cid1[j]][1].Fill(new_tx)
    _ = hy[cid1[j]][1].Fill(new_ty)
    new_x = new_tx*cell_size[int(cid1[j]<1000)][0]+module_pos[cid1[j]][0]
    new_y = new_ty*cell_size[int(cid1[j]<1000)][1]+module_pos[cid1[j]][1]
    _ = hxy[1].Fill(new_x,new_y)
    _ = hE[cid1[j]][1].Fill(E[1])
    _ = helas[cid1[j]][1].Fill(elas_m)
    _ = hxE[cid1[j]][1].Fill(new_tx,elas_m+1)
    _ = hyE[cid1[j]][1].Fill(new_ty,elas_m+1)
    _ = hxyE[1].Fill(new_x,new_y,elas_m+1)
    _ = hxyE[1].Fill(new_x,new_y,elas_m+1)
    _ = hthE[1].Fill(thetadeg[j][1],E[j])
    if correct:
      E[j] = E[j]*correct_energy(new_tx,1,cid1[j],0)*correct_energy(new_ty,1,cid1[j],1)
      elas_m = E[j]/Em[j][1]/1000-1
      _ = hEc[cid1[j]][1].Fill(E[j])
      _ = helasc[cid1[j]][1].Fill(elas_m)
      _ = hxEc[cid1[j]][1].Fill(new_tx,elas_m+1)
      _ = hyEc[cid1[j]][1].Fill(new_ty,elas_m+1)
      _ = hxyEc[1].Fill(new_x,new_y,elas_m+1)
      _ = hxyEc[1].Fill(new_x,new_y,elas_m+1)
      _ = hthEc[1].Fill(thetadeg[j][1],E[j])
        

custom_exit()
