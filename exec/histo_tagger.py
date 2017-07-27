#!/usr/bin/env python
from params import *

# argument parsing
largs = [('r','982','','run_number'),
         ('i','','','input_folder'),
         ('m','island','','clustering_method'),
         ('o',workf,'','output_folder'),
         ('s',0,'','show_bar'),
         ('l',1.0,'','limit'),
         ('w','','','write_tagger_offset'),
         ('b',False,'','batch_mode'),
         ('f','','','path_to_offset_file'),
         ('n',[1,1,1,1],'','new_offsets')]
print_help(largs)
[run,indir,method,outdir,show,limit,write,batch,path,new_offsets] = catch_args(largs)

# batch mode
if batch:
  jsub(project='prad',track='analysis',command=os.path.abspath('histo_tagger.py'),jobname='histo_tagger_'+run,options=del_args(largs,['b'],join=1),memory='1024 MB',os='centos7')
  sys.exit()





from misc_root import *

[names,offsets,otlr,odiff,oet] = readfile(path,[str,float,float,float,float],tr=True)
if new_offsets[0]==0: deads = []
else: deads = [i for i,x in enumerate(offsets) if x==0]
[tmin,tmax] = readfile(workf+'/PRadAnalyzer/config/tagger/tagETCoin-0.9.dat',[int,int],start=1,cols=[5,6],tr=True)

for i in range(512):
  if new_offsets[0]==0: offsets[i] = 0
  if new_offsets[1]==0: otlr[i] = 0
  if new_offsets[2]==0: odiff[i] = 0
  if new_offsets[3]==0: oet[i] = 0

fout = TFile(outdir+'/histo_tagger_'+run+'.root','recreate')
if new_offsets[0]==0: hhits = listofhisto('hhits',names,1000,1000,1500)
else: hhits = listofhisto('hhits',names,1000,-300,300)
htlr = listofhisto('htlr',names[:64],1000,-100,100)
hneigh = listofhisto('hneigh',names,1000,-100,100)
het = listofhisto('het',names,1000,-300,300)
h2 = listofhisto('h2id',[],128,0,128,384,0,384)

f = TFile(indir+'/tree_'+method+'_'+run+'.root')
t = f.event
t.SetBranchStatus("*",0)
for x in ['n_tdc','idg','tg']: t.SetBranchStatus(x,1)

for _ in progress(t,show=show,n=t.GetEntries(),precision=2,limit=limit):
  idg = [x for x in t.idg]
  tg = [x-offsets[y] for x,y in zip(t.tg,t.idg)]
  idg2,tg2,idg3 = [],[],[]
  
  for i in range(t.n_tdc): 
    if idg[i] in deads: continue
    hhits[idg[i]].Fill(tg[i])
  for i in range(t.n_tdc-1):
    if idg[i] in deads: continue
    for j in range(i+1,t.n_tdc):
      if idg[j] in deads: continue
      # tlr match
      if idg[i]<64 and idg[j]<128 and idg[i]%64==idg[j]%64 and idg[i]!=idg[j]: htlr[idg[i]%64].Fill(tg[i]-tg[j]-otlr[idg[i]%64])
      if idg[j]<64 and idg[i]<128 and idg[i]%64==idg[j]%64 and idg[i]!=idg[j]: htlr[idg[i]%64].Fill(tg[j]-tg[i]-otlr[idg[i]%64])
      # neighbor match
      if idg[j]==idg[i]+1: 
        hneigh[idg[i]].Fill(tg[j]-tg[i]-odiff[idg[i]])
        if idg[i]>=128: continue
        idg2.append(idg[i])
        tg2.append((tg[i]+tg[j])/2.)
        idg3.extend([idg[i],idg[j]])
      if idg[i]==idg[j]+1: 
        hneigh[idg[j]].Fill(tg[i]-tg[j]-odiff[idg[j]])
        if idg[j]>=128: continue
        idg2.append(idg[j])
        tg2.append((tg[i]+tg[j])/2.)
        idg3.extend([idg[i],idg[j]])
  # et match
  for i in range(t.n_tdc):
    if idg[i]<128: continue
    if idg[i] in deads: continue
    for j in range(t.n_tdc):
      if idg[j]>=128: continue
      if idg[j] in deads: continue
      if idg[j] in idg3: continue
      h2.Fill(2*(idg[j]%64),idg[i]-128)
      if 2*(idg[j]%64)<tmin[idg[i]-128] or 2*(idg[j]%64)>tmax[idg[i]-128]: continue
      het[idg[i]].Fill(tg[j]-tg[i]-oet[idg[i]])
    for j in range(len(idg2)):
      if idg2[j]>=128: continue
      if idg2[j] in deads: continue
      h2.Fill(2*(idg2[j]%64)+1,idg[i]-128)
      if 2*(idg2[j]%64)+1<tmin[idg[i]-128] or 2*(idg2[j]%64)+1>tmax[idg[i]-128]: continue
      het[idg[i]].Fill(tg2[j]-tg[i]-oet[idg[i]])
       

f.Close()

s = TSpectrum(1)
ffit = TF1('ffit','gaus')

lmean = []
for x in hhits:
  s.Clear()
  s.Search(x)
  if s.GetNPeaks()==0 or x.GetEntries()==0: 
    lmean.append(0)
  else:
    peak = s.GetPositionX()[0]
    ffit.SetParameters(peak,10.,s.GetPositionY()[0])
    ffit.SetParLimits(1,peak-15,peak+15)
    ffit.SetParLimits(2,5.,15.)
    x.Fit(ffit,'rq','',peak-15,peak+15)
    mean = ffit.GetParameter(1)
    x.Fit(ffit,'rq','',mean-10,mean+10)
    lmean.append(ffit.GetParameter(1))

ltlr = []
for x in htlr:
  s.Clear()
  s.Search(x)
  if s.GetNPeaks()==0 or x.GetEntries()==0: 
    ltlr.append(0)
  else:
    peak = s.GetPositionX()[0]
    ffit.SetParameters(peak,1.,s.GetPositionY()[0])
    ffit.SetParLimits(1,peak-5,peak+5)
    ffit.SetParLimits(2,0.1,10.)
    x.Fit(ffit,'rq','',peak-5,peak+5)
    mean = ffit.GetParameter(1)
    x.Fit(ffit,'rq','',mean-5,mean+5)
    ltlr.append(ffit.GetParameter(1))
ltlr.extend([0 for i in range(448)]) 

ldiff = []    
for i,x in enumerate(hneigh):
  s.Clear()
  s.Search(x)
  if s.GetNPeaks()==0 or x.GetEntries()==0: 
    ldiff.append(0)
  else:
    peak = s.GetPositionX()[0]
    ffit.SetParameters(peak,1.,s.GetPositionY()[0])
    ffit.SetParLimits(1,peak-5,peak+5)
    ffit.SetParLimits(2,0.1,10.)
    x.Fit(ffit,'rq','',peak-5,peak+5)
    mean = ffit.GetParameter(1)
    x.Fit(ffit,'rq','',mean-5,mean+5)
    ldiff.append(ffit.GetParameter(1))

let = []    
for i,x in enumerate(het):
  s.Clear()
  s.Search(x)
  if s.GetNPeaks()==0 or x.GetEntries()==0: let.append(0)
  else: let.append(s.GetPositionX()[0])

lmean2,ltlr2,ldiff2,let2 = [],[],[],[]
for i in range(512):
  if new_offsets[0]==0: lmean2.append(lmean[i])
  elif new_offsets[0]==1: lmean2.append(lmean[i]+offsets[i])
  else: lmean2.append(offsets[i])
  
  if new_offsets[1]==0: ltlr2.append(ltlr[i])
  elif new_offsets[1]==1: ltlr2.append(ltlr[i]+otlr[i])
  else: ltlr2.append(otlr[i])

  if new_offsets[2]==0: ldiff2.append(ldiff[i])
  elif new_offsets[2]==1: ldiff2.append(ldiff[i]+odiff[i])
  else: ldiff2.append(odiff[i])

  if new_offsets[3]==0: let2.append(let[i])
  elif new_offsets[3]==1: let2.append(let[i]+oet[i])
  else: let2.append(oet[i])
                                                  
if write!='': writelist(write,[[v,int(round(w,0)),round(x,2),round(y,2),round(z,2)] for v,w,x,y,z in zip(names,lmean2,ltlr2,ldiff2,let2)],ljust=10)

writelisto([h2],fout)
writelisto(hhits,fout,'hits')
writelisto(htlr,fout,'tlr',xrg=[-20,20])
writelisto(hneigh,fout,'neigh',xrg=[-20,20])
writelisto(het,fout,'et',xrg=[-20,20])
fout.Close()
