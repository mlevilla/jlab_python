#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('i',workf+'/trees/trees_42','','input_folder'),
         ('o','.','','output folder'),
         ('n','','','suffix'),
         ('c','','','do_correction'),
         ('nbin',200,'','nbin')]

aparser = ArgParser(largs,sub=1)
[lrun,indir,outdir,suffix,correct,nbin] = aparser.argl

# batch mode ############################################################################

if batch[0]:
  runs = get_runs_between(lrun[0],lrun[-1],'prod')
  jsub_fast(aparser,'density_ene',['r','o'],[['r'],runs],outdir,jobname='density_ene'+suffix,time=wtime[0])
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

if correct: 
  cpar = readfile(correct,[float]*5,cols=range(1,6)) 
  use_regions = 'region' in correct
  nenergy = len(cpar)/(2*nindex)
  
  def correct_density(t,ie,iregion,ix):
    i = iregion*nenergy+ie+len(cpar)/2*ix
    return t+cpar[i][0]*t*(t**4+cpar[i][1]*t**2+cpar[i][2])*(t**2-0.25)*(t**2-cpar[i][3])
  

f = TFile(indir+'/tree_island_'+str(lrun[0])+'.root')
t = f.event
t.SetBranchStatus('*',0)
for x in ['n_cl','id','xhycal','yhycal','zhycal','nh','E']: t.SetBranchStatus(x,1)

nindex = 112

fout = TFile(outdir+'/density_ene_'+str(lrun[0])+suffix+'.root','recreate')
hx = listofhisto('hx',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hy = listofhisto('hy',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxy = listofhisto('hxy',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hE = listofhisto('hE',[range(nindex),range(2)],2500,0,2500,title=';E\' (MeV);')
hxE = listofhisto('hxE',[range(nindex),range(3,5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hyE = listofhisto('hyE',[range(nindex),range(3,5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxyE = listofhisto('hxyE',range(3,5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
if correct:
  hxc = listofhisto('hxc',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hyc = listofhisto('hyc',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
  hxyc = listofhisto('hxyc',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  hxcE = listofhisto('hxcE',[range(nindex),range(3,5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hycE = listofhisto('hycE',[range(nindex),range(3,5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
  hxycE = listofhisto('hxycE',range(3,5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  h2x = listofhisto('h2x',[],2156,0,2156,2000,-2,2,title=';id;(x\'-x)/d')
  h2y = listofhisto('h2y',[],2156,0,2156,2000,-2,2,title=';id;(y\'-y)/d')

# post process ##########################################################################

def custom_exit():
  print
  hx2 = listofhisto('hx2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
  hy2 = listofhisto('hy2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
  hx2E = listofhisto('hx2E',[range(nindex),range(3,5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
  hy2E = listofhisto('hy2E',[range(nindex),range(3,5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
  if correct:
    hxc2 = listofhisto('hxc2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
    hyc2 = listofhisto('hyc2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
    hxc2E = listofhisto('hxc2E',[range(nindex),range(3,5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
    hyc2E = listofhisto('hyc2E',[range(nindex),range(3,5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
  for i in range(nindex):
    for j in range(5):
      hx2[i][j].SetEntries(hx[i][j].GetEntries())
      hy2[i][j].SetEntries(hy[i][j].GetEntries())
      if correct:
        hxc2[i][j].SetEntries(hxc[i][j].GetEntries())
        hyc2[i][j].SetEntries(hyc[i][j].GetEntries())
      for k in range(nbin/2):
        hx2[i][j][k+1] = hx[i][j][k+nbin/2+1] + hx[i][j][nbin/2+1-k]
        hy2[i][j][k+1] = hy[i][j][k+nbin/2+1] + hy[i][j][nbin/2+1-k]
        if correct:
          hxc2[i][j][k+1] = hxc[i][j][k+nbin/2+1] + hxc[i][j][nbin/2+1-k]
          hyc2[i][j][k+1] = hyc[i][j][k+nbin/2+1] + hyc[i][j][nbin/2+1-k]
    for j in range(2):
      hx2E[i][j].SetEntries(hxE[i][j].GetEntries())
      hy2E[i][j].SetEntries(hyE[i][j].GetEntries())
      if correct:
        hxc2E[i][j].SetEntries(hxcE[i][j].GetEntries())
        hyc2E[i][j].SetEntries(hycE[i][j].GetEntries())
      for k in range(nbin/2):
        hx2E[i][j][k+1] = hxE[i][j][k+nbin/2+1] + hxE[i][j][nbin/2+1-k]
        hy2E[i][j][k+1] = hyE[i][j][k+nbin/2+1] + hyE[i][j][nbin/2+1-k]
        if correct:
          hxc2E[i][j][k+1] = hxcE[i][j][k+nbin/2+1] + hxcE[i][j][nbin/2+1-k]
          hyc2E[i][j][k+1] = hycE[i][j][k+nbin/2+1] + hycE[i][j][nbin/2+1-k]
      

  writelisto([hx,hy,hxE,hyE],fout,'raw',style=[1,0,0])
  writelisto([hx2,hy2,hx2E,hy2E],fout,'final',style=[1,0,0])
  writelisto([hxy,hxyE],fout,'occupancy')
  writelisto(hE,fout,'energy',style=[1,0,0])
  if correct:
    writelisto([hxc,hyc,hxcE,hycE],fout,'raw_correct',style=[1,0,0])
    writelisto([hxc2,hyc2,hxc2E,hyc2E],fout,'final_correct',style=[1,0,0])
    writelisto([hxyc,hxycE],fout,'occupancy_correct')
    writelisto([h2x,h2y],fout,'difference')

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
    idx_module = module2index(cid1)
    if correct: 
      iregion = regions_id(idx_module) if use_regions else idx_module
    theta = ftheta([x,y,t.zhycal[k]])
    E = correct_linearity(t.E[k],cid1)
    E_theo = [ep_energy_el(theta,Ebeam),moller_energy(theta,Ebeam)]
    elas = [E/etheo/1000.-1 for etheo in E_theo] 
    elasc = [6*[0.024,0.062][cid1<1000]/etheo**0.5 for etheo in E_theo]
    for ie in range(5):
      if not ((ie<3 and erange[ie]<E<erange[ie+1]) or (ie>=3 and abs(elas[ie-3])<elasc[ie-3])): continue
      _ = hx[idx_module][ie].Fill(tx)
      _ = hy[idx_module][ie].Fill(ty)
      _ = hxy[ie].Fill(x,y)
      if ie>=3: 
        _ = hE[idx_module][ie-3].Fill(E)
        _ = hxE[idx_module][ie-3].Fill(tx,E/E_theo[ie-3])
        _ = hyE[idx_module][ie-3].Fill(ty,E/E_theo[ie-3])
        _ = hxyE[ie-3].Fill(x,y,E/E_theo[ie-3])
      if correct:
        new_tx = correct_density(tx,ie,iregion,0)
        new_ty = correct_density(ty,ie,iregion,1)
        _ = hxc[idx_module][ie].Fill(new_tx)
        _ = hyc[idx_module][ie].Fill(new_ty)
        if ie>=3:
          _ = hxcE[idx_module][ie-3].Fill(new_tx,E/E_theo[ie-3])
          _ = hycE[idx_module][ie-3].Fill(new_ty,E/E_theo[ie-3])
        new_x = new_tx*cell_size[int(cid1<1000)][0]+module_pos[cid1][0]
        new_y = new_ty*cell_size[int(cid1<1000)][1]+module_pos[cid1][1]
        _ = hxyc[ie].Fill(new_x,new_y)
        _ = h2x.Fill(cid1,(new_x-x)/cell_size[int(cid1<1000)][0])
        _ = h2y.Fill(cid1,(new_y-y)/cell_size[int(cid1<1000)][1])
        if ie>=3: _ = hxycE[ie-3].Fill(new_x,new_y,E/E_theo[ie-3])
        

custom_exit()
