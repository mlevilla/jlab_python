from params import *
from moller import *
from esepp import *
from energy_loss import energy_loss
from mdict import *

edge = get_regions('inner')+get_regions('outer')

current = Mdict()

# utils #################################################################################

def get_coordinate(t,i):
  chycal = [t.xhycal[i],t.yhycal[i],t.zhycal[i]]
  cgem = [[] if t.zgem[2*i+j]<=0 else [t.xgem[2*i+j],t.ygem[2*i+j],t.zgem[2*i+j]] for j in range(2)]
  return [chycal]+cgem

def proj(c,z): 
  if c[2]==z: return c
  return [x*z/c[2] for x in c]

def clear_current():
  current['iev'] = 0
  current['idx'] = []
  current['cid'] = []
  current['c0'] = []
  current['cm'] = []
  current['E0'] = []
  current['E1'] = []
  current['E2'] = []
  current['lg'] = []

def print_event(t=None,c=None,E=None,iev=None,typ='ep'):
  if t is None:
    print 'event',iev
    for i in range(len(E)):
      print 'cluster',i
      print 'energy =',E[i]
      print 'hycal coordinate =',c[i][0]
      print 'gem coordinate =',c[i][1]
      print 'angle =',
      print 'E_theo' 
  
  
#### hycal ##############################################################################

def fiducial_cut(cid,c,exclude_edge,exclude_dead,rdead=1.):
  if c==[]: return True
  [x,y,z] = proj(c,zhycal)
  [x2,y2] = to_hycal_frame([x,y])
  if exclude_edge and ((cid in edge) or (abs(x2)<2*cell_size[0][0] and abs(y2)<2*cell_size[0][1]) or (abs(x2)>17*cell_size[0][0]+5*cell_size[1][0]) or (abs(y2)>17*cell_size[0][1]+5*cell_size[1][1])): 
    return True
  if exclude_dead:
    for i,[xdead,ydead] in zip(dead_modules,dead_pos):
      k = int(i<1000)
      if ((x-xdead)**2+(y-ydead)**2)**0.5<rdead*dead_radius[k]: return True
  return False

#### gem ################################################################################

def merge_gem(c):
  if c[1]==c[2]==[]: return [c[0],[]],0
  elif c[1]==[]: return [c[0],c[2]],2
  elif c[2]==[]: return [c[0],c[1]],1
  # else: return [c[0],[0.5*(c[1][i]+c[2][i]) for i in range(3)]]
  else:
    c_proj1 = [proj(x,zhycal) for x in c]
    dr1 = [((x[0]-c_proj1[0][0])**2+(x[1]-c_proj1[0][1])**2)**0.5 for x in c_proj1[1:]]
    c_proj2 = proj(c[1],c[2][2])
    dr2 = ((c_proj2[0]-c[2][0])**2+(c_proj2[1]-c[2][1])**2)**0.5
    if dr2>10*0.08 and dr1[0]>dr1[1]: return [c[0],c[2]],2
    elif dr2>10*0.08: return [c[0],c[1]],1
    else: return [c[0],[0.5*(c[1][i]+c[2][i]) for i in range(3)]],3

def gem_dead(x,y,igem):
  if igem!=2: return False
  return  (21<=x-gem2_center[0]<=434 and -327.05<=y-gem2_center[1]<=-337.05) or (434<=x-gem2_center[0]<=436 and -337.05<=y-gem2_center[1])

def gem_spacers(x,y,igem=0,width=20.):
  if igem==1: x,y = x-gem1_center[0],y-gem1_center[1]
  elif igem==2: x,y = x-gem2_center[0],y-gem2_center[1]
  elif igem==3: x,y = x-0.5*(gem1_center[0]+gem2_center[0]),y-0.5*(gem1_center[1]+gem2_center[1])
  x_spacer = [-344.45,-161.55,161.55,344.45]
  y_spacer = [-409.3,-204,0,204.1,409.4]
  return any([abs(x-xs)<width/2. for xs in x_spacer]) or any([abs(y-ys)<width/2. for ys in y_spacer])

def gem_beam_square(x,y):
  return (abs(x)<50 and abs(y)<50)

def gem_fiducial(c,rdead=1.):
  [x,y,z] = proj(c,zhycal)
  for i,[xdead,ydead] in zip(dead_modules,dead_pos):
    k = int(i<1000)
    if ((x-xdead)**2+(y-ydead)**2)**0.5<rdead*dead_radius[k]: return True
  for [xdead,ydead] in dead_pos_comp:
    if ((x-xdead)**2+(y-ydead)**2)**0.5<15.: return True
  return False
  

#### data ####

def correct_linearity(E,cid):
  return E/(1+linfactor[cid]*(E-ecalib[cid])/1000.)

def get_variables(t,lincorr=0,mgem=1,match=0,exclude_edge=0,exclude_dead=0,rdead=1.,spacer=0,eloss=0):
  #clear_current()
  E,c,idx,lg,cid = [],[],[],[],[]
  for i in range(t.n_cl):
    if t.nh[i]<2: continue
    # coordinate
    ctmp = get_coordinate(t,i)
    #current.c0.append(get_coordinate(t,i))
    # hycal
    if fiducial_cut(t.id[i],ctmp[0],exclude_edge,exclude_dead,rdead): continue
    # gem
    if mgem: ctmp,igem = merge_gem(ctmp)
    if match and ((len(ctmp)==2 and ctmp[1]==[]) or (len(ctmp)==3 and ctmp[1]==ctmp[2]==[])): continue
    if spacer and mgem and (ctmp[1]==[] or gem_spacers(*ctmp[1][:2],igem=igem) or gem_dead(*ctmp[1][:2],igem=igem)): continue
    # energy
    if lincorr: Etmp = correct_linearity(t.E[i],t.id[i])
    else: Etmp = t.E[i]
    if Etmp<0: continue
    # ionization energy loss
    if eloss and ctmp[1]!=[]: Etmp += energy_loss(ftheta(ctmp[1]),Etmp)
    elif eloss: Etmp += energy_loss(ftheta(ctmp[0]),Etmp)
    E.append(Etmp)
    c.append(ctmp)
    idx.append(i)
    lg.append(int(t.id[i]<1000 or is_in_region(t.id[i],'transition_pwo'))) 
    cid.append(t.id[i])
  return [E,c,idx,lg,cid,t.Ebeam/1000.]


#### simulation ####

def correct_linearity_sim(E):
  return E*(exp(-E*1.53438e-04)+1.1133e-4*E+7.17932e-2)

def efficiency_sim(cid):
  if cid==-1: return False
  eff = hycal_trgeff[cid]
  a = random.uniform(0.,1.)
  return a<eff

def match_sim(c,E,lg,sigma=6):
  c_proj = [proj(x,zhycal) for x in c]
  d = ((c_proj[1][0]-c_proj[0][0])**2+(c_proj[1][1]-c_proj[0][1])**2)**0.5
  res = 6.5*sigma/(E/1000)**0.5 if lg else 2.5*sigma/(E/1000)**0.5
  #print c,c_proj,d,res
  return d<res

def get_variables_sim(t,match=0,exclude_edge=0,exclude_dead=0,lincorr=0,hycal_eff=0,rdead=1.,spacer=0,do_matching=0):
  E,c,idx,lg,cid = [],[],[],[],[]
  for i in range(getattr(t,'HC.N')):
    # hycal
    E_hycal = correct_linearity_sim(getattr(t,'HC.P')[i]) if lincorr else getattr(t,'HC.P')[i]
    if E_hycal<30.: continue
    c_hycal = [getattr(t,'HC.'+x)[i]-[0,ztarget][int(x=='Z')] for x in ['X','Y','Z']]
    cid_hycal = getattr(t,'HC.CID')[i]
    lg_hycal = int(cid_hycal<1000 or is_in_region(cid_hycal,'transition_pwo'))
    if hycal_eff and (not efficiency_sim(c_hycal)): continue
    #gem
    c_gem = [getattr(t,'GEM.'+x)[i]-[0,ztarget][int(x=='Z')] for x in ['X','Y','Z']] if getattr(t,'GEM.Z')[i]>0 else [] # maxime
    if match and c_gem==[]: continue
    if fiducial_cut(cid_hycal,c_hycal,exclude_edge,exclude_dead,rdead): continue
    if match and do_matching and (not match_sim([c_hycal,c_gem],E_hycal,lg_hycal)): continue
    if spacer and (c_gem==[] or gem_spacers(*c_gem[:2],igem=0) or gem_dead(*c_gem[:2],igem=2)): continue

    idx.append(i)
    cid.append(cid_hycal)
    E.append(E_hycal)
    c.append([c_hycal,c_gem])
    lg.append(lg_hycal)
  return [E,c,idx,lg,cid,0]

#### selection ####
 
def has_br(c):
  if len(c)==0: return False
  elif len(c[0])==2: return any([x[1]==[] for x in c])
  else: return  any([(x[1]==[] and x[2]==[]) for x in c])


def is_ep(theta=None,E=None,Ebeam=None,c=None,E_theo=None,elas=None,elasc=None,lg=None,sigma=None,thetacut=0.):
  if theta is None and c==[]: return False
  if theta is None: theta = ftheta(c) 
  if theta<thetacut*degrad: return False
  if elas is None:
    if E_theo is None: E_theo = ep_energy_el(theta,Ebeam)
    elas = E/E_theo/1000.-1
  elasc = sigma*[0.024,0.062][lg]/E_theo**0.5 if elasc is None else elasc
  return abs(elas)<elasc

def is_ee1(theta=None,E=None,Ebeam=None,c=None,E_theo=None,elas=None,elasc=None,lg=None,sigma=None,thetacut=0.):
  if theta is None and c==[]: return False
  if theta is None: theta = ftheta(c) 
  if theta<thetacut*degrad: return False
  if elas is None:
    if E_theo is None: E_theo = moller_energy(theta,Ebeam)
    elas = E/E_theo/1000.-1
  elasc = sigma*[0.024,0.062][lg]/E_theo**0.5 if elasc is None else elasc
  return abs(elas)<elasc

def is_ee2(theta=None,E=None,Ebeam=None,c=None,E_theo=None,elas=None,elasc=None,lg=None,sigma=None,dphi=None,zvertex=None,gem=[0,1],phicut=[10,5],zcut=[500,150],dEcut=None,thetacut=[0.,0.]):
  # double arm
  if dphi is None and theta is None and any(x[gem[0]]==[] for x in c): return [0,0,0]
  if dphi is None: dphi = degrees((fphi(c[0][gem[0]])-fphi(c[1][gem[0]]))%(2*pi)-pi)
  theta_flag = False
  if theta is None: 
    theta_flag = True
    theta = [ftheta(x[gem[0]]) for x in c]
  etheo_flag = False
  if E_theo is None:
    etheo_flag = True
    E_theo = [moller_energy(th,Ebeam) for th in theta]
  dE = sum(E)-Ebeam*1000
  if dEcut is None: dEcut = sigma[0]*(([0.024,0.062][lg[0]]*E_theo[0]**0.5)**2+([0.024,0.062][lg[1]]*E_theo[1]**0.5)**2)**0.5*1000
  if zvertex is None:
    cproj = [[x[gem[0]][m]*zhycal/x[gem[0]][2] for m in range(2)] for x in c]
    rj = [(x[0]**2+x[1]**2)**0.5 for x in cproj]
    zvertex = ((m_e+Ebeam)*rj[0]*rj[1]/2./m_e)**0.5-zhycal
  if abs(dphi)>phicut[gem[0]] or abs(zvertex)>zcut[gem[0]] or abs(dE)>dEcut or any(th<thetacut[gem[0]]*degrad for th in theta): return [0,0,0]
  # single arm
  if theta_flag: theta = [ftheta(x[gem[1]]) if x[gem[1]]!=[] else -1 for x in c]
  if elas is None:
    if etheo_flag: E_theo = [moller_energy(th,Ebeam) if th!=-1 else -1 for th in theta]
    elas = [e/e_theo/1000.-1 for e,e_theo in zip(E,E_theo)]
  if elasc is None: elasc = [sigma[1]*[0.024,0.062][i]/e_theo**0.5 if e_theo!=-1 else 0 for i,e_theo in zip(lg,E_theo)]
  return [1]+[abs(x)<y or th>thetacut[gem[1]]*degrad  if th!=-1 else False for x,y,th in zip(elas,elasc,theta)]
