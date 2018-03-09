#!/usr/bin/env python
from params import *

largs = [('f','','','file_name'),
         ('nf',1,'','number_of_files'),
         ('o','.','','output folder'),
         ('n','','','suffix'),
         ('c','','','do_correction'),
         ('nbin',200,'','nbin')]

aparser = ArgParser(largs,sub=1)
[filename,nf,outdir,suffix,correct,nbin] = aparser.argl

# batch mode ############################################################################

if batch[0]:
  jsub_fast(aparser,'density_sim',['nf','o'],[['n','f'],['_'+str(x) for x in range(nf)]],outdir,jobname='density_sim'+suffix,time=wtime[0])
  sys.exit()

from misc_root import *
from moddensity_utils import *
from selection import *

# initialisation ########################################################################

run = 2002 if '2GeV' in filename else 2001
print run
load_positions(run=run) 
load_run_calib(run)
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
  

f = TFile(filename+'_rec.root')
t = f.T
t.SetBranchStatus('*',0)
for x in ['HC.N','HC.CID','HC.X','HC.Y','HC.Z','HC.P']: t.SetBranchStatus(x,1)

nindex = 112

fout = TFile(outdir+'/density_sim'+suffix+'.root','recreate')
hx = listofhisto('hx',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
hy = listofhisto('hy',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
hxy = listofhisto('hxy',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
hE = listofhisto('hE',[range(nindex),range(2)],2500,0,2500,title=';E\' (MeV);')
if correct:
  hxc = listofhisto('hxc',[range(nindex),range(5)],nbin,-1.,1.,title=';(x_{hycal}-x_{center})/d;')
  hyc = listofhisto('hyc',[range(nindex),range(5)],nbin,-1.,1.,title=';(y_{hycal}-y_{center})/d;')
  hxyc = listofhisto('hxyc',range(5),1200,-600,600,1200,-600,600,title=';x (mm);y (mm)')
  h2x = listofhisto('h2x',[],2156,0,2156,2000,-2,2,title=';id;(x\'-x)/d')
  h2y = listofhisto('h2y',[],2156,0,2156,2000,-2,2,title=';id;(y\'-y)/d')

# post process ##########################################################################

def custom_exit():
  hx2 = listofhisto('hx2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
  hy2 = listofhisto('hy2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
  if correct:
    hxc2 = listofhisto('hxc2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(x_{hycal}-x_{center})/d;')
    hyc2 = listofhisto('hyc2',[range(nindex),range(5)],nbin/2,0.,1.,title=';(y_{hycal}-y_{center})/d;')
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
      
      if sum(hx2[i][j][k] for k in range(nbin/2))!=0:
        hx2[i][j].Scale(float(nbin/4)/sum(hx2[i][j][k] for k in range(nbin/2)))
      if sum(hy2[i][j][k] for k in range(nbin/2))!=0:
        hy2[i][j].Scale(float(nbin/4)/sum(hy2[i][j][k] for k in range(nbin/2)))
      if correct:
        if sum(hxc2[i][j][k] for k in range(nbin/2))!=0:
          hxc2[i][j].Scale(float(nbin/4)/sum(hxc2[i][j][k] for k in range(nbin/2)))
        if sum(hy2[i][j][k] for k in range(nbin/2))!=0:
          hyc2[i][j].Scale(float(nbin/4)/sum(hyc2[i][j][k] for k in range(nbin/2)))
              

  writelisto([hx,hy],fout,'raw',style=[1,0,0])
  writelisto([hx2,hy2],fout,'final',style=[1,0,0])
  writelisto(hxy,fout,'occupancy')
  writelisto(hE,fout,'energy',style=[1,0,0])
  if correct:
    writelisto([hxc,hyc],fout,'raw_correct',style=[1,0,0])
    writelisto([hxc2,hyc2],fout,'final_correct',style=[1,0,0])
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
  for k in range(getattr(t,'HC.N')):
    cid,x,y = getattr(t,'HC.CID')[k], getattr(t,'HC.X')[k], getattr(t,'HC.Y')[k]
    if cid in edge+[1561,1562,1595,1596]: continue
    cid,tx,ty = correct_module_id(cid,x,y)
    if cid in edge+deads: continue
    idx_module = module2index(cid)
    if correct: 
      iregion = regions_id(idx_module) if use_regions else idx_module
    theta = ftheta([x,y,getattr(t,'HC.Z')[k]-89])
    E = getattr(t,'HC.P')[k]
    E_theo = [ep_energy_el(theta,Ebeam),moller_energy(theta,Ebeam)]
    elas = [E/etheo/1000.-1 for etheo in E_theo] 
    elasc = [6*[0.024,0.062][cid<1000]/etheo**0.5 for etheo in E_theo]
    for ie in range(5):
      if not ((ie<3 and erange[ie]<E<erange[ie+1]) or (ie>=3 and abs(elas[ie-3])<elasc[ie-3])): continue
      _ = hx[idx_module][ie].Fill(tx)
      _ = hy[idx_module][ie].Fill(ty)
      _ = hxy[ie].Fill(x,y)
      if ie>=3: _ = hE[idx_module][ie-3].Fill(E)
      if correct:
        new_tx = correct_density(tx,ie,iregion,0)
        new_ty = correct_density(ty,ie,iregion,1)
        _ = hxc[idx_module][ie].Fill(new_tx)
        _ = hyc[idx_module][ie].Fill(new_ty)
        new_x = new_tx*cell_size[int(cid<1000)][0]+module_pos[cid][0]
        new_y = new_ty*cell_size[int(cid<1000)][1]+module_pos[cid][1]
        _ = hxyc[ie].Fill(new_x,new_y)
        _ = h2x.Fill(cid,(new_x-x)/cell_size[int(cid<1000)][0])
        _ = h2y.Fill(cid,(new_y-y)/cell_size[int(cid<1000)][1])
        

custom_exit()
