#!/apps/python/2.7.12/bin/python
from misc import *

run = catch_arg('r',1345)
suffix = catch_arg('n','')

if '-b' in sys.argv:
  jsub(project='prad',track='analysis',jobname='test_zcarbon_'+str(run),command=os.path.abspath('test_zcarbon.py -r '+str(run)),memory='2 GB')
  sys.exit()

from params import *
from misc_root import *
from geo3D import *
from selection import *

fout = TFile('{0}/test_zcarbon2_{1}.root'.format(work,run),'recreate')
f = TFile('{0}/tree/island/tree_island_{1}.root'.format(work,run))
t = f.event

load_run_params(int(run))
edge = get_regions('inner')+get_regions('outer')
zgem0 = [5175,5215]

# histograms
hzhycal = [TH1F('hzhycal_'+str(i),';z_{vertex} (mm);',1000,4000,6500) for i in range(2)]
hzgem = [TH1F('hzgem_'+str(i),';z_{vertex} (mm);',1000,4000,6500) for i in range(2)]
hzhycal2 = TH1F('hzhycal2',';z_{vertex} (mm);',1000,4000,6500)
hzgem2 = [TH1F('hzgem2_'+str(i),';z_{vertex} (mm);',1000,4000,6500) for i in range(2)]
hzvse = [TH2F('hzvse_'+str(i),';E (GeV);z_{hycal};',500,0.,2.5,500,5500,6000) for i in range(2)]
hdxhycal = [TH1F('hdxhycal_'+str(i),';#Deltax_{hycal} (mm);',1000,-50,50) for i in range(2)]
hdyhycal = [TH1F('hdyhycal_'+str(i),';#Deltay_{hycal} (mm);',1000,-50,50) for i in range(2)]

for _ in progress(t,show=2,n=t.GetEntries(),modulo=10000):
  #preselection
  [E,c,idx,Ebeam] = get_variables(t,1,['t.n_cl>1','(t.id[icl] not in '+str(edge)+' and t.E[icl]>50)','not has_br([c])'],0)

  for i1 in range(len(idx)-1):
    for i2 in range(i1+1,len(idx)):
      cj,Ej,idxj = [c[i1],c[i2]],[E[i1],E[i2]],[idx[i1],idx[i2]]
      elas1 = elasticity(cj,Ej,Ebeam,0)
      elas2 = elasticity_moller(cj,Ej,Ebeam,0)
      dphi = degrees(deltaphi(cj,0))
      if abs(elas1)>0.1 or abs(elas2)>0.1 or abs(dphi)>5.: continue
      
      # hycal
      rhycal = [(x[0][0]**2+x[0][1]**2)**0.5 for x in cj]
      dzhycal = cj[0][0][2]-cj[1][0][2]
      zhycal = ((m_e+beam_energy/1000.)*rhycal[0]*rhycal[1]/2./m_e)**0.5
      zhycal2 = [((m_e+beam_energy/1000.)*rhycal[0]*rhycal[1]/2./m_e+dzhycal**2/4.)**0.5-dzhycal/2.-cj[0][0][2]+5640,((m_e+beam_energy/1000.)*rhycal[0]*rhycal[1]/2./m_e+dzhycal**2/4.)**0.5+dzhycal/2.-cj[1][0][2]+5640]                       
      if t.id[idxj[0]]>1000 and t.id[idxj[1]]>1000: _ = hzhycal[0].Fill(zhycal)
      else: _ = hzhycal[1].Fill(zhycal)
      for j in range(2): hzhycal2.Fill(zhycal2[j])

      # gem
  
      if cj[0][1]!=[] and cj[1][1]!=[]:
        rgem1 = [(x[1][0]**2+x[1][1]**2)**0.5 for x in cj]
        zgem1 = ((m_e+beam_energy/1000.)*rgem1[0]*rgem1[1]/2./m_e)**0.5
        _ = hzgem[0].Fill(zgem1)
      if cj[0][2]!=[] and cj[1][2]!=[]:
        rgem2 = [(x[2][0]**2+x[2][1]**2)**0.5 for x in cj]
        zgem2 = ((m_e+beam_energy/1000.)*rgem2[0]*rgem2[1]/2./m_e)**0.5
        _ = hzgem[1].Fill(zgem2)
      for j1 in range(1,3):
        for j2 in range(1,3):
          if cj[0][j1]!=[] and cj[1][j2]!=[]:
            rgem = [(x[0]**2+x[1]**2)**0.5 for x in [cj[0][j1],cj[1][j2]]]
            dzgem = cj[0][j1][2]-cj[1][j2][2]
            zgem = [((m_e+beam_energy/1000.)*rgem[0]*rgem[1]/2./m_e+dzgem**2/4.)**0.5-dzgem/2,((m_e+beam_energy/1000.)*rgem[0]*rgem[1]/2./m_e+dzgem**2/4.)**0.5+dzgem/2]
            for k in range(2): _ = hzgem2[[j1,j2][k]-1].Fill(zgem[k])
              
      # z hycal
      for j in range(2):
        for k in range(2):
          if cj[j][k+1]==[]: continue
          theta = angle(cj[j][k+1])
          Emoller = moller_energy(theta,beam_energy/1000.)
          inter = clst_l2l(cj[j][k+1],[0,0,0],cj[j][0][:2]+[0.],cj[j][0][:2]+[1.])
          m = int(t.id[idxj[j]]<1000)
          _ = hzvse[m].Fill(Emoller,0.5*(inter[0][2]+inter[1][2]))
          _ = hdxhycal[m].Fill(inter[0][0]-inter[1][0])
          _ = hdyhycal[m].Fill(inter[0][1]-inter[1][1])
        

fout.cd()
for x in hzhycal: x.Write()
for x in hzgem: x.Write()
for x in hzgem2: x.Write()
hzhycal2.Write()
for x in hdxhycal: x.Write()
for x in hdyhycal: x.Write()
for x in hzvse: x.Write()
fout.Close()
