#!/usr/bin/env python

from params import *
from misc_root import *

if len(sys.argv)==1:
  input_files = add_args_list(sys.argv,'i',['_'+str(x) for x in range(300)],write=jobf+'/gun_theta')
  jsub(project='prad',track='analysis',jobname='gun_theta',command=execf+'/gun_theta.py args.in',input_files=input_files,input_data='args.in',memory='4096 MB',os='centos7')
  sys.exit()

largs = [('i',-1,'','run_number')]
[number] = catch_args(largs)

fout = TFile(workf+'/gun_theta_distribution'+str(number)+'.root','recreate')
names = [['/work/hallb/prad/mlevilla/new_evgen_files/ep/1GeV/ep_1GeV_1e6_th04_75_CO1_1000_FF5','/work/hallb/prad/mlevilla/new_evgen_files/ep/2GeV/ep_2GeV_1e6_th04_75_CO1_2000_FF5'],['/work/hallb/prad/mlevilla/new_evgen_files/moller/1GeV/moller_1GeV_1e6_th04_75_CO1_1000','/work/hallb/prad/mlevilla/new_evgen_files/moller/2GeV/moller_2GeV_1e6_th04_35_CO1_2000']]

h = listofhisto('hth',[['ep','ee'],['1GeV','2GeV'],['elas','rad']],10000,0,10)
h2 = listofhisto('hEth',[['ep','ee'],['1GeV','2GeV'],['elas','rad']],1000,0,10,1000,0,2500)

for i1 in range(2):
  for i2 in range(2):
    name = names[i1][i2]+'_'+str(number)+'_e-.dat' if i1==0 else names[i1][i2]+'_'+str(number)+'.dat'
    [l_E_e1,l_th_e1,l_E_e2,l_th_e2,l_E_g,l_th_g] = readfile(name,[float]*6,cols=[0,1,3,4,6,7],tr=True)
    for j in range(1000000):
      if l_E_e2[j]!=l_E_e2[j]: continue
      for k in range(2):
        _=h[i1][i2][k].Fill(l_th_e1[j]/degrad)
        _=h2[i1][i2][k].Fill(l_th_e1[j]/degrad,l_E_e1[j])
        if i1==1: 
          _=h[i1][i2][k].Fill(l_th_e2[j]/degrad)
          _=h2[i1][i2][k].Fill(l_th_e2[j]/degrad,l_E_e2[j])
      if l_E_g[j]>0: 
        _=h[i1][i2][1].Fill(l_th_g[j]/degrad)
        _=h2[i1][i2][1].Fill(l_th_g[j]/degrad,l_E_g[j])

writelisto([h,h2],fout)
fout.Close()
