#!/usr/bin/env python

from params import *

largs = [('elascut',4,'','elas_cut'),
         ('phicut',10,'','phi_cut'),
         ('zcut',200,'','zvertex_cut'),
         ('o',work+'/acceptance','','outdir'),
         ('n','geo_acceptance','','file_name'),
         ('b',False,'','batch_mode'),
         ('s',0,'','progress_bar'),
         ('suff','_res','','suffix')]

print_help(largs)
[elascut,phicut,zcut,outdir,filename,batch,show,suffix] = [catch_arg(x,y,z) for x,y,z,c in largs]

if batch:
  jobname = 'acceptance_{0}_{1}_{2}'.format(elascut,phicut,zcut)
  args = '-elascut {0} -phicut {1} -zcut {2} -o {3} -n {4} -suff {5}'.format(elascut,phicut,zcut,outdir,filename,suffix)
  jsub(project='prad',track='analysis',command=os.path.abspath('acceptance.py'),jobname=jobname,options=args,memory='2 GB')
  sys.exit()

from selection import *

fout = TFile(outdir+'/'+filename+'.root','recreate')
f = [[TFile(work+'/simulation/'+x+'_'+e+suffix+'.root') for x in ['moller','ep']] for e in ['11','22']]
#f = [[TFile(work+'/pradsim-1/output/moller_elastic_11.root')],[TFile(work+'/pradsim-1/output/moller_elastic_22.root')]]
lt = [[x.T for x in y] for y in f] 
lt = [[x[0],x[0],x[1]] for x in lt]
#t = [[x[0],x[0]] for x in t]
h = listofhisto('h',[['11','22'],['2arm','1arm','ep'],['all','cut','gem','hc','both','gem_cut','hc_cut','both_cut']],200,0.,10.,fout=fout,title=';#theta (deg);acceptance')

for k1,x1 in enumerate(lt):
  Ebeam = [1.1,2.2][k1]
  for k2,t in enumerate(x1):
    count = [0 for i in range(5)]
    for _ in enumerate(progress(t,show=show,precision=1,n=t.GetEntries())):
      # gun
      n = getattr(t,'GUN.N')
      idgen = [i for i in range(n) if getattr(t,'GUN.PID')[i]==11]
      n2 = len(idgen)
      theta_gen = [getattr(t,'GUN.Theta')[i] for i in idgen]
      E_theo = [getattr(t,'GUN.E')[i] for i in idgen]
      # gem acceptance
      n_gem = getattr(t,'GEM.N')
      idx_gem = [getattr(t,'GEM.TID')[i]-1 if (getattr(t,'GEM.PTID')[i]==0 and getattr(t,'GEM.PID')[i]==11) else -1 for i in range(n_gem)] 
      c_gem = [[[getattr(t,'GEM.'+x)[j]+[0,2911][int(x=='Z')] for x in ['X','Y','Z']] for j in range(n_gem) if i==idx_gem[j]] for i in idgen]
      c_gem = [[sum([x[j][i] for j in range(len(x))])/len(x) if len(x)!=0 else -10000 for i in range(3)] for x in c_gem]
      n_gem2 = len([x for x in c_gem if x[0]!=10000])
      theta_gem = [ftheta(x) if x[2]!=-10000 else 0. for x in c_gem]
      if n_gem2!=2:
        dphi_gem = 0.
        z_gem = 0.
      else:
        dphi_gem = degrees((fphi(c_gem[0])-fphi(c_gem[1]))%(2*pi)-pi)
        r_gem = [(x[0]**2+x[1]**2)**0.5 for x in c_gem]
        dz = c_gem[0][2] - c_gem[1][2]
        z_gem = ((m_e+Ebeam)*r_gem[0]*r_gem[1]/2./m_e+dz**2/4)**0.5 - c_gem[0][2] + dz/2.
      # hc acceptance
      n_hc = getattr(t,'HC.N')
      idx_hc = [getattr(t,'HC.TID')[i]-1 if (getattr(t,'HC.PTID')[i]==0 and getattr(t,'HC.PID')[i]==11) else -1 for i in range(n_hc)] 
      c_hc = [[[getattr(t,'HC.'+x)[j]+[0,2911][int(x=='Z')] for x in ['X','Y','Z']] for j in range(n_hc) if i==idx_hc[j]] for i in idgen]
      c_hc = [[sum([x[j][i] for j in range(len(x))])/len(x) if len(x)!=0 else -10000 for i in range(3)] for x in c_hc]
      n_hc2 = len(c_hc)
      theta_hc = [ftheta(x) for x in c_hc]
      E_hc = [[getattr(t,'HC.P')[j] for j in range(n_hc) if idx_hc[j]==i] for i in idgen] 
      E_hc = [x[0] if x!=[] else -1 for x in E_hc]
      elas1 = sum(E_hc)/sum(E_theo)-1
      elas2 = [x/y-1 for x,y in zip(E_hc,E_theo)]
      cid2 = [int(abs(c_hc[i][0])>17*20.75 and abs(c_hc[i][1])>17*20.75) for i in range(n2)]
      elasc2 = [elascut*[0.025,0.06][cid2[i]]/(E_theo[i]/1000)**0.5 for i in range(n2)]
      cid1 = int(any(cid2))
      elasc1 = elascut*[0.025,0.06][cid1]/Ebeam**0.5
      fiducial = [int((abs(c_hc[i][0])>2*20.75 or abs(c_hc[i][1])>2*20.75) and abs(c_hc[i][0])<524.4 and abs(c_hc[i][1])<524.4) if c_hc[i]!=[] else 0 for i in range(n2)]
      n_hc2 = sum(fiducial)
      if n_hc2!=2:
        dphi_hc = 0.
        z_hc = 0.
      else:
        dphi_hc = degrees((fphi(c_hc[0])-fphi(c_hc[1]))%(2*pi)-pi)
        r_hc = [(x[0]**2+x[1]**2)**0.5 for x in c_hc]
        dz = c_hc[0][2] - c_hc[1][2]
        z_hc = ((m_e+Ebeam)*r_hc[0]*r_hc[1]/2./m_e+dz**2/4)**0.5 - c_hc[0][2] + dz/2.
      # histo filling
      for i in range(n2): h[k1][k2][0].Fill(degrees(theta_gen[i]))
      if k2==0:
        if abs(elas1)<elasc1 and abs(dphi_gem)<phicut and abs(z_gem)<zcut and (not any([abs(e)>c for e,c in zip(elas2,elasc2)])): 
          for i in range(n2): h[k1][k2][1].Fill(degrees(theta_gen[i]))
        if all([c[2]!=-10000 for c in c_gem]): 
          for i in range(n2): h[k1][k2][2].Fill(degrees(theta_gen[i]))
        if all([c[2]!=-10000 for c in c_hc]) and all(fiducial): 
          for i in range(n2): h[k1][k2][3].Fill(degrees(theta_gen[i]))
        if all([c[2]!=-10000 for c in c_hc]) and all([c[2]!=-10000 for c in c_gem]) and all(fiducial):
          for i in range(n2): h[k1][k2][4].Fill(degrees(theta_gen[i]))
        if abs(elas1)<elasc1 and abs(dphi_gem)<phicut and abs(z_gem)<zcut and all([c[2]!=-10000 for c in c_gem]):
          for i in range(n2): h[k1][k2][5].Fill(degrees(theta_gen[i]))
        if abs(elas1)<elasc1 and abs(dphi_hc)<phicut and abs(z_hc)<zcut and (not any([abs(e)>c for e,c in zip(elas2,elasc2)])) and all([c[2]!=-10000 for c in c_hc]) and all(fiducial): 
          for i in range(n2): h[k1][k2][6].Fill(degrees(theta_gen[i]))
        if abs(elas1)<elasc1 and abs(dphi_gem)<phicut and abs(z_gem)<zcut and (not any([abs(e)>c for e,c in zip(elas2,elasc2)])) and all([c[2]!=-10000 for c in c_hc]) and all([c[2]!=-10000 for c in c_gem]) and all(fiducial):
          for i in range(n2): h[k1][k2][7].Fill(degrees(theta_gen[i]))
      else:
        for i in range(n2): 
          if abs(elas2[i])<elasc2[i]: h[k1][k2][1].Fill(degrees(theta_gen[i]))
          if c_gem[i][2]!=-10000: h[k1][k2][2].Fill(degrees(theta_gen[i]))
          if c_hc[i][2]!=-10000 and fiducial[i]: h[k1][k2][3].Fill(degrees(theta_gen[i]))
          if c_gem[i][2]!=-10000 and c_hc[i][2]!=-10000 and fiducial[i]: h[k1][k2][4].Fill(degrees(theta_gen[i]))
          if abs(elas2[i])>elasc2[i]: continue
          if c_gem[i][2]!=-10000: h[k1][k2][5].Fill(degrees(theta_gen[i]))
          if c_hc[i][2]!=-10000 and fiducial[i]: h[k1][k2][6].Fill(degrees(theta_gen[i]))
          if c_gem[i][2]!=-10000 and c_hc[i][2]!=-10000 and fiducial[i]: h[k1][k2][7].Fill(degrees(theta_gen[i]))        
    for i in range(7): h[k1][k2][i+1].Divide(h[k1][k2][0])

writelisto(h,fout)
fout.Close()
