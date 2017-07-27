#!/usr/bin/env python
from misc import *

run = catch_arg('r',1345)
name = catch_arg('n','')
limit = catch_arg('l',1.)
show = catch_arg('s',0)

if '-b' in sys.argv:
  command = '{0} -r {1}{2}'.format(os.path.abspath('test_zcarbon2.py'),run,(' -n '+name)*int(name!=''))
  jsub(project='prad',track='analysis',jobname='test_zcarbon2_'+str(run),command=command,memory='2 GB')
  sys.exit()

from params import *
from misc_root import *
from geo3D import *
from selection import *

load_run_params(run)

fout = TFile('{0}/test_zcarbon2_{1}{2}.root'.format(work,run,name),'recreate')
f = TFile('{0}/tree/island/tree_island_{1}.root'.format(work,run))
t = f.event

hz = TH1F('hz','',400,-500,500)

m_e = 0.51099893e-3
for _ in progress(t,show=show,n=t.GetEntries(),precision=2,limit=limit):
  [E,c,idx,Ebeam] = get_variables(t,lincorr=1,mgem=0,exclude_edge=1,match=1)
  if len(c)!=2: continue
  if c[0][1]==c[1][1]: continue
  if any([x==[] for x in c[0]]) or any([x==[] for x in c[1]]): continue
  r = [[(x[0]**2+x[1]**2)**0.5 for x in y] for y in c]
  theta = [abs(degrees(ftheta(x[1]))) for x in c]
  elas = elasticity(c,E,Ebeam,1,theta=theta)
  dphi = (fphi(c[0][0])-fphi(c[1][0]))%(2*pi)-pi
  if abs(elas)>0.1 or degrees(abs(dphi))>5.: continue
  if any([x<1.2 for x in theta]): continue
  center = [((m_e+t.Ebeam/1000.)*r[0][i]*r[1][i]/2./m_e)**0.5-c[0][i][2] for i in range(1,3)]
  hz.Fill((center[0]+center[1])/2.)
  

fout.cd()
hz.Write()
fout.Close()
f.Close()
sys.exit()
