#!/bin/env python

from params import *

largs = [('i',workf,'','input_directory'),
         ('n','','several','suffix'),
         ('s',3.,'','sigma_of_fit')]
print_help(largs)
[indir,suffix,sigma] = catch_args(largs)

from misc_root import *
gStyle.SetOptFit(1)
load_names()

f = TFile(indir+'/histo'+suffix[0]+'.root')
fout = TFile(indir+'/reso'+suffix[-1]+'.root','recreate')

nbinE = 18

hr = {k:[[f.Get(x).Get('hr_{0}_{1}_{2}'.format(x,y,z)) for z in range(nbinE)] for y in range(1,3)] if f.Get(x) else [] for k,x in module_names.items()}
hrx = [[[f.transition.Get('hrx_'+str(i)+'_'+str(j)+'_'+str(k)) for k in range(nbinE)] for j in range(20)] for i in range(8)]


le = [200+(900/nbinE/2.)+(900/nbinE)*j for j in range(nbinE)]

weizhi_params = readfile(workf+'/resolution_curve_eSigma4_ex1.txt',[str,float,float,float],out=dict)

# reso by modules
lerror,lalpha,lreso,lmean,lsigma,ldmean,ldsigma = [],[],[],[],[],[],[]
lgreso = []
for x,y in progress(module_names.items(),precision=0,show=3,insert=lambda x:x[1]):
  # print y
  lerr,lm,ls,ldm,lds = [],[],[],[],[]
  for j in range(nbinE):
    working = True
    if hr[x]==[] or x in [900,1835] or (1476<x<1775 and j*(900/nbinE)+200>650): working = False
    else:
      hr[x][0][j].Add(hr[x][1][j])
      # [m,s,dm,ds,_] = fit_gaussian(hr[x][0][j],200.,sigma,[0.5,1.5])
 
      [m,s,dm,ds,_] = fit_2gaus(hr[x][0][j],100.,[0.5,2.])
      if m<0. or dm==0. or s>0.3: working = False
      else:
        lm.append(m)
        ls.append(s)
        ldm.append(dm)
        lds.append(ds)
    if not working:
      lerr.append(j)
      lm.append(0)
      ls.append(0)
      ldm.append(0)
      lds.append(0)

  if len(lerr)==nbinE:
    print ('no stat at all for module '+y).ljust(89)
    lalpha.append([0,0]);  lreso.append([0,0,0,0,0,0,0,0]);
    lgreso.append(None)
  else:
    #sel = [j for j in range(2,8)+range(11,16) if lm[j]!=0]
    sel = range(2,10)+range(11,16) if i!=1 else range(1,9)+range(11,14)
    sel = [j for j in sel if lm[j]!=0]
    # [galpha,_,alpha,_,dalpha] = fit_list('[0]+[1]*x*0.001',le,lm,dy=ldm,sel=sel)
    alpha, dalpha = 0,0
    [greso,greso2,reso1,reso2,reso3,dreso1,dreso2,dreso3] = fit_list('sqrt(([0]/sqrt(x*0.001))^2 + ([1]/(x*0.001))^2 + [2]^2)',le,ls,dy=lds,xrg=[250,1050],sel=sel,start=weizhi_params[y])
    reso = (reso1**2+reso2**2+reso3**2)**0.5
    dreso = 0
    if reso1!=0: dreso+=(dreso1/reso1)**2
    if reso2!=0: dreso+=(dreso2/reso2)**2
    if reso3!=0: dreso+=(dreso3/reso3)**2
    if dreso>0: dreso = reso*dreso**0.5
    if abs(alpha)>0.5: alpha = 0
    if reso<0.015 or reso>0.3: reso = 0
    lalpha.append([alpha,dalpha])
    lreso.append([reso1,reso2,reso3,dreso1,dreso2,dreso3,reso,dreso])
    greso.SetName('greso_'+y)
    lgreso.append(greso)

  lmean.append(lm)
  lsigma.append(ls)
  ldmean.append(ldm)
  ldsigma.append(lds)
  lerror.append(lerr)

for x,y in zip(['mean','sigma','dmean','dsigma','error','reso','alpha'],[lmean,lsigma,ldmean,ldsigma,lerror,lreso,lalpha]):
  writelist(indir+'/'+x+suffix[-1]+'_'+str(int(sigma))+'.txt',[[a]+b for a,b in zip(module_names.values(),y)],ro=5,ljust=8)


writelisto([lgreso],fout,['modules',module_names.values()])
writelisto([[x[0] if x!=[] else [] for x in hr.values()]],fout,['modules',module_names.values()])

# transition region

# lreso,lsigma,ldsigma,lmean,ldmean,lalpha = [],[],[],[],[],[]
# lname = [x+'_'+str(j) for x in ['top','right','bottom','left'] for j in range(20)]
# for i in range(4):
#   for j in range(20):
#     lerr,lj,le,lm,ls,ldm,lds = [],[],[],[],[],[],[]
#     for k in range(9):
#       [m,s,dm,ds] = fit_gaussian(hrx[i][j][k],2.,sigma,0.8,2.)
#       if m<0. or dm==0.: 
#         lerr.append(j)
#         continue
#       lj.append(k); le.append(250+100*k); lm.append(m); ls.append(s); ldm.append(dm); lds.append(ds)
    
#     if lj==[]: lalpha.append([0,0]);  lreso.append([0,0])
#     else:
#       [galpha,_,alpha,_,dalpha] = fit_list('[0]+[1]*x*0.001',le,lm,dy=ldm)
#       [greso,reso,dreso] = fit_list('[0]/sqrt(x*0.001)',le,ls,dy=lds)
#       if abs(alpha)>0.5: alpha = -1
#       if reso<0.015 or reso>0.3: reso = 0
#       lreso.append([reso,dreso])
#       lalpha.append([alpha,dalpha])
#     lsigma.append([ls[lj.index(j)] if (j in lj and ls[lj.index(j)]<0.3) else 0 for j in range(9)])
#     ldsigma.append([lds[lj.index(j)] if j in lj else 0 for j in range(9)])
#     lmean.append([lm[lj.index(j)] if j in lj else 0 for j in range(9)])
#     ldmean.append([ldm[lj.index(j)] if j in lj else 0 for j in range(9)])

# writelist(indir+'/d'+name+'/trans_mean'+filename+'.txt',[[x]+y for x,y in zip(lname,lmean)])
# writelist(indir+'/d'+name+'/trans_sigma'+filename+'.txt',[[x]+y for x,y in zip(lname,lsigma)])
# writelist(indir+'/d'+name+'/trans_dmean'+filename+'.txt',[[x]+y for x,y in zip(lname,ldmean)])
# writelist(indir+'/d'+name+'/trans_dsigma'+filename+'.txt',[[x]+y for x,y in zip(lname,ldsigma)])
# writelist(indir+'/d'+name+'/trans_reso'+filename+'.txt',[[x]+y for x,y in zip(lname,lreso)])
# writelist(indir+'/d'+name+'/trans_alpha'+filename+'.txt',[[x]+y for x,y in zip(lname,lalpha)])

# particular region
name_regions = ['trans','center_lg','center_pwo']
regions = [get_regions(x) for x in name_regions]
hr_regions = [[TH1F('hr_'+x+'_'+str(j),';E_{cluster}/E_{#gamma};',1000,0.,5.) for j in range(nbinE)] for x in name_regions]
for x,y in module_names.items():
  for i,z in enumerate(regions):
    if x in z and hr[x]!=[]:
      for j in range(nbinE): hr_regions[i][j].Add(hr[x][0][j])

lreso,lsigma,ldsigma,lmean,ldmean,lalpha = [],[],[],[],[],[]
lgreso = []
for i in range(len(regions)):
  lerr,lm,ls,ldm,lds = [],[],[],[],[]
  for k in range(nbinE):
    # [m,s,dm,ds,_] = fit_gaussian(hr_regions[i][k],1000.,sigma,[0.5,1.5])
    [m,s,dm,ds,_] = fit_2gaus(hr_regions[i][k],1000.,[0.5,2.])
    if m<0. or dm==0. or s>0.3: 
      lerr.append(j)
      lm.append(0)
      ls.append(0)
      ldm.append(0)
      lds.append(0)
    else:
      lm.append(m)
      ls.append(s)
      ldm.append(dm);
      lds.append(ds)
    
  if len(lerr)==nbinE:
    lalpha.append([0,0])
    lreso.append([0,0])
  else:
    #sel = range(2,8)+range(11,15)
    sel = range(2,10)+range(11,16) if i!=1 else range(1,9)+range(11,14)
    [galpha,galpha2,_,alpha,_,dalpha] = fit_list('[0]+[1]*x*0.001',le,lm,dy=ldm,sel=sel)
    [greso,greso2,reso1,reso2,reso3,dreso1,dreso2,dreso3] = fit_list('sqrt(([0]/sqrt(x*0.001))^2 + ([1]/(x*0.001))^2 + [2]^2)',le,ls,dy=lds,xrg=[250,1050],sel=sel)
    # [greso,reso1,reso2,reso3,dreso1,dreso2,dreso3] = fit_list('sqrt(([0]/sqrt(x*0.001))^2+([1]/(x*0.001))^2+([2]/(x*0.001)^2)^2)',le,ls,dy=lds,xrg=[250,1050],sel=sel)
    reso = (reso1**2+reso2**2+reso3**2)**0.5
    dreso = 0
    if reso1!=0: dreso+=(dreso1/reso1)**2
    if reso2!=0: dreso+=(dreso2/reso2)**2
    if reso3!=0: dreso+=(dreso3/reso3)**2
    if dreso>0: dreso = reso*dreso**0.5
    if abs(alpha)>0.5: alpha = 0
    if reso<0.015 or reso>0.3: reso = 0
    lreso.append([reso1,reso2,reso3,dreso1,dreso2,dreso3,reso,dreso])
    greso.SetName('greso_'+name_regions[i])
    greso2.SetName('greso2_'+name_regions[i])
    lgreso.append([greso,greso2])
    lalpha.append([alpha,dalpha])

  lmean.append(lm)
  lsigma.append(ls)
  ldmean.append(ldm)
  ldsigma.append(lds)

for x,y in zip(['mean','sigma','dmean','dsigma','reso','alpha'],[lmean,lsigma,ldmean,ldsigma,lreso,lalpha]):
  writelist(indir+'/region_'+x+suffix[-1]+'_'+str(int(sigma))+'.txt',[[a]+b for a,b in zip(module_names.values(),y)],ro=5,ljust=8)

writelisto([[x,y] for x,y in zip(lgreso,hr_regions)],fout,[name_regions])

fout.Close()

# resolution versus location inside the module
# lreso,lsigma,ldsigma,lmean,ldmean,lalpha = [],[],[],[],[],[]
# lname = [x+'_'+str(j) for x in ['x_lg','y_lg','x_pwo','y_pwo'] for j in range(20)]
# for i in range(4):
#   for j in range(20):
#     lerr,lj,le,lm,ls,ldm,lds = [],[],[],[],[],[],[]
#     for k in range(9):
#       [m,s,dm,ds] = fit_gaussian(hrx[i+4][j][k],2.,sigma,0.8,2.)
#       if m<0. or dm==0.: 
#         lerr.append(j)
#         continue
#       lj.append(k); le.append(250+100*k); lm.append(m); ls.append(s); ldm.append(dm); lds.append(ds)
    
#     if lj==[]: lalpha.append([0,0]);  lreso.append([0,0])
#     else:
#       [galpha,_,alpha,_,dalpha] = fit_list('[0]+[1]*x*0.001',le,lm,dy=ldm)
#       [greso,reso,dreso] = fit_list('[0]/sqrt(x*0.001)',le,ls,dy=lds)
#       if abs(alpha)>0.5: alpha = -1
#       if reso<0.015 or reso>0.3: reso = 0
#       lreso.append([reso,dreso])
#       lalpha.append([alpha,dalpha])
#     lsigma.append([ls[lj.index(j)] if (j in lj and ls[lj.index(j)]<0.3) else 0 for j in range(9)])
#     ldsigma.append([lds[lj.index(j)] if j in lj else 0 for j in range(9)])
#     lmean.append([lm[lj.index(j)] if j in lj else 0 for j in range(9)])
#     ldmean.append([ldm[lj.index(j)] if j in lj else 0 for j in range(9)])

# writelist(indir+'/d'+name+'/inside_mean'+filename+'.txt',[[x]+y for x,y in zip(lname,lmean)])
# writelist(indir+'/d'+name+'/inside_sigma'+filename+'.txt',[[x]+y for x,y in zip(lname,lsigma)])
# writelist(indir+'/d'+name+'/inside_dmean'+filename+'.txt',[[x]+y for x,y in zip(lname,ldmean)])
# writelist(indir+'/d'+name+'/inside_dsigma'+filename+'.txt',[[x]+y for x,y in zip(lname,ldsigma)])
# writelist(indir+'/d'+name+'/inside_reso'+filename+'.txt',[[x]+y for x,y in zip(lname,lreso)])
# writelist(indir+'/d'+name+'/inside_alpha'+filename+'.txt',[[x]+y for x,y in zip(lname,lalpha)])
