#!/usr/bin/env python
from params import *

largs = [('r',-1,'several','run_number'),
         ('i',volf+'/trees/trees_sel6','','input_folder'),
         ('o','.','','output folder'),
         ('n','','','suffix'),
         ('c','','','do_correction'),
         ('nbin',200,'','nbin')]

aparser = ArgParser(largs,sub=1)
[lrun,indir,outdir,suffix,correct,nbin] = aparser.argl

# batch mode ############################################################################

if batch[0]:
  runs = get_runs_between(lrun[0],lrun[-1],'prod')
  jsub_fast(aparser,'density',['r','o'],[['r'],runs],outdir,jobname='density'+suffix,time=wtime[0])
  sys.exit()

from misc_root import *
from moddensity_utils import *
from selection import *

# initialisation ########################################################################

load_positions(run=lrun[0]) 
load_run_calib(lrun[0])
Ebeam = beam_energy[0]/1000.
erange = erange_2GeV if Ebeam>2. else erange_1GeV
nindex = 101

if correct: 
  cpar = readfile(correct,[float]*5,cols=range(1,6)) 
  use_regions = 'region' in correct
  nenergy = len(cpar)/(4*nindex)
  
  def correct_density(t,ie,iregion,ix,ipos):
    i = 2*(iregion*nenergy+ie)+ipos+len(cpar)/2*ix
    return t+cpar[i][0]*t*(t**4+cpar[i][1]*t**2+cpar[i][2])*(t**2-0.25)*(t**2-cpar[i][3])
  

f = TFile(indir+'/tree_island_'+str(lrun[0])+'.root')
t = f.event
t.SetBranchStatus('*',0)
for x in ['n_cl','id','xhycal','yhycal','zhycal','nh','E']: t.SetBranchStatus(x,1)

fout = TFile(outdir+'/density_'+str(lrun[0])+suffix+'.root','recreate')
hx = listofhisto('hx',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hy = listofhisto('hy',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')

hxy = listofhisto('hxy',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
if correct:
  hxc = listofhisto('hxc',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hyc = listofhisto('hyc',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')


  hxyc = listofhisto('hxyc',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  h2x = listofhisto('h2x',[],2156,0,2156,2000,-2,2,title=';id;(x\'-x)/d')
  h2y = listofhisto('h2y',[],2156,0,2156,2000,-2,2,title=';id;(y\'-y)/d')

# post process ##########################################################################

def custom_exit():
  print
  for i in range(nindex):
    for j in range(5):
      if sum(hx[i][j][k+1] for k in range(nbin))!=0:
        hx[i][j].Scale(float(nbin/2)/sum(hx[i][j][k+1] for k in range(nbin)))
      if sum(hy[i][j][k+1] for k in range(nbin))!=0:
        hy[i][j].Scale(float(nbin/2)/sum(hy[i][j][k+1] for k in range(nbin)))

      if correct:
        if sum(hxc[i][j][k] for k in range(nbin))!=0:
          hxc[i][j].Scale(float(nbin/2)/sum(hxc[i][j][k] for k in range(nbin)))
        if sum(hyc[i][j][k] for k in range(nbin))!=0:
          hyc[i][j].Scale(float(nbin/4)/sum(hyc[i][j][k] for k in range(nbin)))

  writelisto([hx,hy],fout,'raw',style=[1,0,0])
  writelisto(hxy,fout,'occupancy')
  if correct:
    writelisto([hxc,hyc],fout,'raw_correct',style=[1,0,0])
    writelisto(hxyc,fout,'occupancy_correct')
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
    idx_module = module2index2(cid1)
    idx_module = index_map[idx_module]
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
      if correct:
        new_tx = correct_density(tx,ie,iregion,0,int(tx>0))
        new_ty = correct_density(ty,ie,iregion,1,int(ty>0))
        _ = hxc[idx_module][ie].Fill(new_tx)
        _ = hyc[idx_module][ie].Fill(new_ty)
        new_x = new_tx*cell_size[int(cid1<1000)][0]+module_pos[cid1][0]
        new_y = new_ty*cell_size[int(cid1<1000)][1]+module_pos[cid1][1]
        _ = hxyc[ie].Fill(new_x,new_y)
        _ = h2x.Fill(cid1,(new_x-x)/cell_size[int(cid1<1000)][0])
        _ = h2y.Fill(cid1,(new_y-y)/cell_size[int(cid1<1000)][1])
        

custom_exit()
