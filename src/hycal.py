from params import *
from math import log,pi,atan2,atan

time_start = time.time()

# tables for cell_hyc
llg = readfile(dataf+'/prof_lg.dat',[int,int,float,float])
lpwo = readfile(dataf+'/prof_pwo.dat',[int,int,float,float])
acell = [[[0. for i in range(501)] for j in range(501)] for k in range(2)]
ad2c = [[[0. for i in range(501)] for j in range(501)] for k in range(2)]
for x in lpwo: 
  acell[0][x[0]][x[1]] = x[2]
  acell[0][x[1]][x[0]] = x[2]
  ad2c[0][x[0]][x[1]] = x[3]
  ad2c[0][x[1]][x[0]] = x[3]
for x in llg: 
  acell[1][x[0]][x[1]] = x[2]
  acell[1][x[1]][x[0]] = x[2]
  ad2c[1][x[0]][x[1]] = x[3]
  ad2c[1][x[1]][x[0]] = x[3]

bad_modules = [148,710,720,799]
offsets = [0.,0.]
cell_size = [[2.077,2.075],[3.815,3.815]]
trans = [35.309,35.275]

########################
### Useful functions ###
#######################

# def g_f(x,y,b): return atan(x/b)+atan(y/b)+atan(x*y/(b*(b*b+x*x+y*y)**0.5))
  
# def cell_hyc(dx,dy,k):
#  if k==0: a,b = 0.98, 1.042
#  elif k==1: a,b = 1.00, 1.206
#  else: a,b = 0.,0.1
#  return a/pi/2.*(g_f(dx+0.5,dy+0.5,b)-g_f(dx+0.5,dy-0.5,b)-g_f(dx-0.5,dy+0.5,b)+g_f(dx-0.5,dy-0.5,b))

def cell_hyc(dx,dy,k):
  ax,ay = 100.*abs(dx),100.*abs(dy)
  i,j = int(ax),int(ay)
  if i>499 or j>499: return 0.
  wx,wy = ax-i,ay-j
  return acell[k][i][j]*(1-wx)*(1-wy)+acell[k][i+1][j]*wx*(1-wy)+acell[k][i][j+1]*(1-wx)*wy+acell[k][i+1][j+1]*wx*wy 

def d2c(dx,dy,k):
  ax,ay = 100.*abs(dx),100.*abs(dy)
  i,j = int(ax),int(ay)
  if i>499 or j>499: return 0.
  wx,wy = ax-i,ay-j
  return ad2c[k][i][j]*(1-wx)*(1-wy)+ad2c[k][i+1][j]*wx*(1-wy)+ad2c[k][i][j+1]*(1-wx)*wy+ad2c[k][i+1][j+1]*wx*wy

def sigma2(dx,dy,fc,e,k):
  if dx**2+dy**2>25: return 1.
  return 0.816*fc+(32.1+1.72*(e/10.)**0.5)*d2c(dx,dy,k)+0.2/e*10.

def distance_center_cls(c1=None,c2=None,x=None,y=None):
  if c1==None and c2==None: x1,x2,y1,y2 = x[0],x[1],y[0],y[1]
  else: x1,x2,y1,y2 = c1.x, c2.x, c1.y, c2.y
  a = (y2-y1)
  b = x1-x2
  c = -a*x1-b*y1
  d = abs(c)/(a**2+b**2)**0.5
  x = -a*c/(a**2+b**2)
  y = -b*c/(a**2+b**2)
  d1 = ((x1-x)**2+(y1-y)**2)**0.5
  d2 = ((x2-x)**2+(y2-y)**2)**0.5
  return d,d1,d2,x,y
  
  
#################
### HyCal Hit ###
#################
  
class hchit:
  def __init__(self,mid=-1,adc=0.,offsetx=0.,offsety=0.,gains=None,lms_gains=None,peds=None,sigma_peds=None,E=None):
    self.mid = mid
    self.adc = adc
    self.E = E if E!=None else 0.
    self.x = offsetx
    self.y = offsety
    self.lg = -1
    if self.mid!=-1:
      self.x += xpos[self.mid]
      self.y += ypos[self.mid]
      self.lg = 0 if (abs(self.x)<trans[0] and abs(self.y)<trans[1]) else 1
      if self.mid<1728 and (self.mid not in bad_modules) and E==None:
        lms = lms_gains[self.mid] if lms_gains!=None else 1.
        ped = peds[self.mid] if peds!=None else 0.
        sigma_ped = sigma_peds[self.mid] if sigma_peds!=None else 0.
        gain = gains[self.mid] if gains!=None else 1.
        self.E = round(adc-ped)*gain*lms if (adc>(ped+5*sigma_ped)) else 0.
      
  def __repr__(self): return 'Hit(mod={0},value={1},x={2},y={3},E={4})'.format(module_names_more[self.mid],self.adc,round(self.x,2),round(self.y,2),round(self.E,2))

  def distance(self,hit): return ((self.x-hit.x)**2+(self.y-hit.y)**2)**0.5

  def isNeighbor(self,cluster): return (self.mid in cluster.neighbors)


#####################
### HyCal Cluster ###
#####################
  
class hccluster:
  def __init__(self,hit):
    self.mid = hit.mid
    self.x = hit.x
    self.y = hit.y
    self.lg = 0 if (abs(self.x)<trans[0] and abs(self.y)<trans[1]) else 1
    self.nhit = 1
    self.hits = [hit]
    self.E = hit.E
    self.chi2 = 1000000.
    self.neighbors = set(neighbors[hit.mid])
    self.leaks = []
    self.status = 1
      
  def __repr__(self): return 'Cluster(center={0},nhit={1},x={2},y={3},E={4},chi2={5})'.format(module_names_more[self.mid],self.nhit,round(self.x,2),round(self.y,2),round(self.E,2),round(self.chi2,2))
      
  def add_hit(self,hit,r=1.,reset=False):
    hit_in_self = -1
    for i,h in enumerate(self.hits):
      if h.mid == hit.mid:
        hit_in_self = i
        break
    if hit_in_self!=-1:
      self.hits[i].E = hit.E*r if reset else self.hits[i].E+hit.E*r
      self.E = hit.E*r if reset else self.E+hit.E*r
    else:
      self.hits.append(hit)
      self.E += hit.E*r 
      self.nhit += 1
      self.neighbors = self.neighbors | set(neighbors[hit.mid]) - set([h.mid for h in self.hits])

  def update(self):
    t0 = time.time()
    self.lg = 0 if (abs(self.x)<trans[0] and abs(self.y)<trans[1]) else 1
    for h in self.hits:
      size_i = cell_size[h.lg]
      if abs(h.x-self.x)<cell_size[h.lg][0] and abs(h.y-self.y)<cell_size[h.lg][1]:
        self.mid = h.mid
        break

  def get_info(self,leak=0):
    i,ecorr0 = 0,0.
    t0 = time.time()
    self.coord()
    if leak:
      ecorr = self.correct_leakage()
      while self.status and abs(ecorr-ecorr0)>0.001 and i<50:
        self.coord()
        ecorr0 = ecorr
        ecorr = self.correct_leakage()
        i+=1
    if i==50: self.status=0
    else: self.chi2 = self.chi2_xy(fast=leak)
    if leak: self.update()
            
  def correct_leakage(self):
    self.leaks = []
    for i in self.neighbors: 
      if (i not in bad_modules) and i<1728: continue
      h = hchit(i)
      h.E = self.cell_hyc(h)
      self.leaks.append(h)
    ecorr = sum([h.E for h in self.leaks])
    if ecorr<0.2: 
      self.E=sum([h.E for h in self.hits])/(1-ecorr)
      for h in self.leaks: h.E*=self.E 
    else: self.status = 0
    return ecorr 
    
  def coord(self):
    if self.E==0.: return
    sum_weight = 0.
    self.x, self.y = 0.,0.
    for h in self.hits+self.leaks:
      weight =  4.6 + log(h.E/self.E) if (self.E*h.E!=0.) else 0.
      if weight>0:
        self.x += weight*h.x 
        self.y += weight*h.y
        sum_weight += weight
    self.x/=sum_weight
    self.y/=sum_weight

  def cell_hyc(self,h,x=None,y=None):
    x = self.x if x==None else x
    y = self.y if y==None else y
    if (self.lg==0 and h.lg==0) or (self.lg==1 and h.lg==1): return cell_hyc((h.x-x)/cell_size[self.lg][0],(h.y-y)/cell_size[self.lg][1],self.lg)
    else:
      sides = [[(x-trans[0])*(h.x-trans[0])<0,(x+trans[0])*(h.x+trans[0])<0],[(y-trans[1])*(h.y-trans[1])<0,(y+trans[1])*(h.y+trans[1])<0]]
      ax = ((-1)**(sides[0].index(1))*trans[0]-x)/(h.x-x) if sum(sides[0])==1 else -1
      ay = ((-1)**(sides[1].index(1))*trans[1]-y)/(h.y-y) if sum(sides[1])==1 else -1
      if ax!=-1 and ay!=-1:
        x0 = x+ay*(h.x-x)
        y0 = y+ax*(h.y-y)
        if abs(x0)<trans[0]: a = ax
        elif abs(y0)<trans[1]: a = ay
        else: return 0.
      elif ax!=-1: a = ax
      else: a = ay
      return a*cell_hyc((a/cell_size[self.lg][0]+(1-a)/cell_size[h.lg][0])*(h.x-x),(a/cell_size[self.lg][1]+(1-a)/cell_size[h.lg][1])*(h.y-y),self.lg)+(1-a)*cell_hyc((a/cell_size[self.lg][0]+(1-a)/cell_size[h.lg][0])*(h.x-x),(a/cell_size[self.lg][1]+(1-a)/cell_size[h.lg][1])*(h.y-y),h.lg)

  def chi2_xy(self,x=None,y=None,fast=0):
    if self.E<=0: return 0.
    if x==None: x = self.x
    if y==None: y = self.y
    s = 0.
    l = self.hits if fast else self.hits+[hchit(neigh) for neigh in self.neighbors]
    for h in l:
      dx,dy = (x-h.x)/cell_size[self.lg][0], (y-h.y)/cell_size[self.lg][1]
      fc = self.cell_hyc(h,x,y)
      if (fc>0. or h.E>0): #s += (fc-h.E/self.E)**2
        s += self.E/10.*(fc-h.E/self.E)**2/sigma2(dx,dy,fc,self.E,self.lg)
    return s

  def best_chi2(self,stepx=0.05,stepy=0.05,cut=0.02):
    i = 0
    while (stepx>cut or stepy>cut) and i<50:
      chil = self.chi2_xy(self.x-stepx,self.y)
      chir = self.chi2_xy(self.x+stepx,self.y)
      chiu = self.chi2_xy(self.x,self.y+stepy)
      chid = self.chi2_xy(self.x,self.y-stepy)
      dx = stepx*(chil-chir)/2.
      dy = stepy*(chid-chiu)/2.
      stepx = min([stepx,abs(dx)])
      stepy = min([stepy,abs(dy)])
      self.x += dx
      self.y += dy
      i+=1
    self.chi2 = self.chi2_xy()
    self.update()
      
  def distance(self,hit): return ((self.x-hit.x)**2+(self.y-hit.y)**2)**0.5

  def isNeighbor(self,hit): return (hit.mid in self.neighbors)

  def merge(self,cluster):
    for h in cluster.hits: self.add_hit(h) 
    for h in cluster.leaks: self.leaks.append(h)

  def copy(self):
    c = hccluster(self.hits[0])
    for h in self.hits[1:]: c.add_hit(h)
    return c


###################
### HyCal Event ###
###################
  
class hcevent:
  def __init__(self,ev,clustering=0,gains=None,peds=None,sigma_peds=None,lms_gains=None,method='island'):
    t0 = time.time()
    self.hits = [hchit(ev.adc_id[j],ev.adc_val[j],gains=gains,lms_gains=lms_gains,peds=peds) for j in range(ev.n_adc) if ev.adc_id[j]<1728]
    self.hits = [h for h in self.hits if h.E>1.]
    self.hits.sort(key=lambda x:x.E,reverse=True)
    self.clusters = []
    #print 1,time.time()-t0
    if clustering>0: 
      if method=='island': self.island()
      elif method=='5by5': self.fiveby5()
    #print 2,time.time()-t0
    if clustering>1 and method=='island': self.merge_clusters()
    #print 3,time.time()-t0
    if clustering>0: 
      for c in self.clusters: c.get_info(1)
      #print 4,time.time()-t0

  def __repr__(self): return 'Event(nhits={0},nclus={1},E={2})'.format(len(self.hits),len(self.clusters),round(self.Etot,2))

  def island(self):
    rh = []
    for i,h in enumerate(self.hits):
      rc = []
      for j,c in enumerate(self.clusters):
        if h.isNeighbor(c): rc.append(j)
      if rc==[]: self.clusters+=[hccluster(h)]
      elif len(rc)==1: self.clusters[rc[0]].add_hit(h)
      else: rh.append((i,rc))
    for i,rc in rh:
      r = [self.clusters[j].E*self.clusters[j].cell_hyc(self.hits[i]) for j in rc]
      s = sum(r)
      for j,f in zip(rc,r): 
        if f>0: self.clusters[j].add_hit(self.hits[i],f/s)
    self.clusters.sort(key=lambda x:x.E,reverse=True)

  def fiveby5(self):
    left_out = []
    for h in self.hits:
      newcluster = 1
      for c in self.clusters:
        diag = 2*1.02*(cell_size[c.lg][0]**2+cell_size[c.lg][1]**2)**0.5
        if h.distance(c.hits[0])<diag:
          c.add_hit(h)
          newcluster = 2
          break
        elif h.distance(c.hits[0])<2*diag:
          newcluster = 0
      if newcluster==1: self.clusters+=[hccluster(h)]
      elif newcluster==0: left_out.append(h)
    for h in left_out: 
      for c in self.clusters:
        if h.distance(c.hits[0])<2*1.02*(cell_size[c.lg][0]**2+cell_size[c.lg][1]**2)**0.5:
          c.add_hit(h)
          break
    self.clusters.sort(key=lambda x:x.E,reverse=True)

  def merge_clusters(self):
    i = 0
    while i<len(self.clusters)-1:
      x = self.clusters[i]
      x.get_info(0)
      dof_x = max(1,len(x.hits)+len(x.neighbors)-2)
      j = i+1 
      while j<len(self.clusters):
        y = self.clusters[j]
        y.get_info(0)
        dof_y = max(1,len(y.hits)+len(y.neighbors)-2)
        if (set([h.mid for h in x.hits]) & y.neighbors)!=set():
          xy = x.copy()
          xy.merge(y)
          xy.get_info(0)
          dof_xy = max(1,len(xy.hits)+len(xy.neighbors)-2)
          if xy.chi2/dof_xy<(x.chi2+y.chi2)/(dof_x+dof_y): 
            self.clusters[i] = xy
            x = self.clusters[i]
            del self.clusters[j]
          else: j+=1
        else: j+=1
      i+=1

  def moller_cls(self,r=0.1,sym=1.,dist=5.,ebeam=1100.):
    l = []
    for i1,c1 in enumerate(self.clusters):
      if c1.nhit<3: continue
      for c2 in self.clusters[i1+1:]:
        if c2.nhit<3: continue
        if abs(c1.E+c2.E-ebeam)/ebeam>r:continue
        d = distance_center_cls(c1,c2)[0]
        if d>dist: continue
        if abs(c1.E-c2.E)/(c1.E+c2.E)/2.>sym: continue
        l.append([c1,c2])
    return l

  def ep_cls(self,r=0.1,ebeam=1100.):
    l = []
    for c in self.clusters:
      if c.nhit<3: continue
      if abs(c.E-ebeam)/ebeam>r:continue
      l.append(c)
    return l

sys.stdout.write('hycal.py loaded in '+str(round(time.time()-time_start,2))+' seconds\n')
