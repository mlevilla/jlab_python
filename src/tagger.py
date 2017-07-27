from misc import *
from params import *
from ROOT import TH1I,TH1F,TH2F,TCanvas,TTree,TFile
from array import array

config_folder = home_folder+'/config_tagger'
[tag_names,means,tlr,diff,etdiff] = readfile(config_folder+'/tagger_offsets.txt',[str,float,float,float,float],cols=range(5),tr=True)
kavg = readfile(config_folder+'/tagE-boundaries-1.1.dat',[float],start=1,cols=[3])
[kminfirst,kmin] = readfile(config_folder+'/tagT-boundaries-0.9.dat',[float,float],start=1,cols=[2,3],tr=True)
[tchmin,tchmax] = readfile(config_folder+'/tagETCoin-0.9.dat',[int,int],start=1,cols=[5,6],tr=True)



class TDCHit:
  def __init__(self,channel_id=-1,time=0.,hit_type=0,energy=0.):
    self.cid = channel_id
    self.t = time
    self.htype = hit_type
    self.E = energy
  def __eq__(self,tdc): return (self.cid==tdc.cid and self.t==tdc.t)
  def __ne__(self,tdc): return (not self==tdc)
  def __lt__(self,tdc): return (self.cid<tdc.cid or self.t<tdc.t)
  def __le__(self,tdc): return (self.cid<=tdc.cid)
  def __gt__(self,tdc): return (self.cid>tdc.cid or self.t>tdc.t)
  def __ge__(self,tdc): return (self.cid>=tdc.cid)
  def __repr__(self): return 'Hit(id={0},type={1},t={2},E={3})'.format(self.cid,self.htype,self.t,self.E)
  
class TDCHits(list):
  def __init__(self,l=[]):
    super(TDCHits,self).__init__(l)
  def items(self):
    for x in self: yield x.cid,x.t
  def ids(self): return [x.cid for x in self]
  def merge(self):
    l = []
    for x in self:
      matched = False
      for i,y in enumerate(l):
        if abs(x.cid-y[-1].cid)<4 and abs(x.t-y[-1].t)<5.:
          l[i].append(x)
          matched = True
          break
      if not matched: l.append([x])
    l2 = TDCHits()
    for x in l:
      l2.append(TDCHit(int(sum([y.cid for y in x])/len(x)),sum([y.t for y in x])/len(x),5,sum([y.E for y in x])/len(x)))
    return l2


  
class HistTDC():
  def __init__(self,name='',nt=12000,tmin=-3000,tmax=3000,ne=120,emin=0.,emax=1.2):
    self.name = name
    self.event = 0
    self.nhit = 0
    self.cid = TH1I('h'+name+'_id','',900,0,900)
    self.t = TH1F('h'+name+'_time','',nt,tmin,tmax)
    self.tdiff = TH1F('h'+name+'_timediff','',nt,tmin,tmax)
    self.E = TH1F('h'+name+'_E','',ne,emin,emax)
    self.EvsEh = TH2F('h'+name+'_EvsEh','',200,0.,1.,200,0.,2.)
    self.tvsEh = TH2F('h'+name+'_tvsEh','',100,0.,3.,600,-300,300)
    self.n = TH1I('h'+name+'_n','',20,0,20)

  def fill(self,tdc,event=0,tref=None,energy=None):
    if type(tdc)==list:
      for x in tdc: self.fill(x)  
    else:
      self.cid.Fill(tdc.cid)
      self.t.Fill(tdc.t)   
      self.E.Fill(tdc.E)
      if tref!=None: self.tdiff.Fill(tdc.t-tref)
      if energy!=None:
        self.tvsEh.Fill(energy,tdc.t-tref)
        self.EvsEh.Fill(energy,tdc.E)
      if event==self.event: self.nhit+=1
      else:
        self.n.Fill(self.nhit)
        self.event = event
        self.nhit = 1
      

  def draw(self,trange=50,idrange=900,c=None):
    if c==None: c = TCanvas()
    c.Divide(2,2)
    c.cd(1)
    self.cid.GetXaxis().SetRangeUser(0,idrange)
    self.cid.Draw()
    c.cd(2)
    center = self.tdiff.GetBinCenter(self.tdiff.GetMaximumBin())
    self.tdiff.GetXaxis().SetRangeUser(center-trange,center+trange)
    self.tdiff.Draw()
    c.cd(3)
    self.tvsEh.GetYaxis().SetRangeUser(center-trange,center+trange)
    self.tvsEh.Draw('colz')
    c.cd(4)
    self.EvsEh.Draw('colz')
    return c

  def write(f):
    dump = f.cd()
    self.cid.Write()
    self.t.Write()
    self.tdiff.Write()
    self.E.Write()
    self.tvsEh.Write()
    self.EvsEh.Write()
    self.n.Write()


    

class TaggerTree(TTree):
  def __init__(self,name='',title='',tree=None):
    if tree!=None: self=tree
    else:
      super(TaggerTree,self).__init__(name,title)
      self.SetDirectory(0)
      dinit = [('iev','i',''),('latch','i',''),('trigger','i',''),('ttrg','f',''),('etot','f',''),('n_ch','i',''),('ch_id','i','[n_ch]'),('ch_type','i','[n_ch]'),('ch_t','f','[n_ch]'),('ch_e','f','[n_ch]')]
      self.nmax = 150
      dtmp = {'i':[0],'f':[0.],'i[':[0]*self.nmax,'f[':[0.]*self.nmax}
      self.dvar = {}
      for x,y,z in dinit:
        self.dvar[x] = array(y,dtmp[y+z[:1]])
        self.Branch(x,self.dvar[x],x+z+'/'+y.upper())

  def setvar(self,iev=-1,latch=-1,trigger=-1,ttrg=0.,etot=0.,ch=TDCHits()):
    self.dvar['iev'][0] = iev
    self.dvar['latch'][0] = latch
    self.dvar['trigger'][0] = trigger
    self.dvar['ttrg'][0] = ttrg
    self.dvar['etot'][0] = etot
    self.dvar['n_ch'][0] = len(ch)
    for j,x in enumerate(ch):
      if j>self.nmax: break
      self.dvar['ch_id'][j] = x.cid
      self.dvar['ch_type'][j] = x.htype
      self.dvar['ch_t'][j] = x.t
      self.dvar['ch_e'][j] = x.E

  def getch(self): return TDCHits([TDCHit(self.ch_id[j],self.ch_t[j],self.ch_type[j],self.ch_e[j]) for j in range(self.n_ch)])

  def __len__(self): return self.GetEntries()



  
#get all hits with raw channel id
def get_hits(ev,raw=False):
  lh = TDCHits()
  for j in range(ev.n_tdc):
    if ev.tdc_id[j]<1000: continue
    if raw: lh.append(TDCHit(ev.tdc_id[j]-1000,ev.tdc_val[j]/10.,0))
    else: lh.append(TDCHit(ev.tdc_id[j]-1000,ev.tdc_val[j]/10.-means[ev.tdc_id[j]-1000],0))
  return lh


def dispatch_hits(lh=None,ev=None,raw=False):
  if lh==None: lh = get_hits(ev,raw)
  tlh,trh,eh = TDCHits(),TDCHits(),TDCHits()
  for j,v in lh.items():
    if j<64: tlh+=[TDCHit(j+1,v,0)]
    elif j<128: trh+=[TDCHit(j-64+1,v,1)]
    else: eh+=[TDCHit(j-128+1,v,5)]
  eh.sort()
  return tlh,trh,eh

def get_thits(tlh=None,trh=None,lh=None,ev=None,raw=False):
  if tlh==None: tlh,trh,_ = dispatch_hits(lh,ev,raw)
  th = TDCHits()
  for j,v in tlh.items():
    if (j in trh.ids()) and (abs(v-trh[trh.ids().index(j)].t-tlr[j-1])<3.): th.append(TDCHit(j,(v+trh[trh.ids().index(j)].t)/2.,2))
    elif j in [17,18,19,55]: th.append(TDCHit(j,v,0))
  for j,v in trh.items():
    if j in [20,32,33]: th.append(TDCHit(j,v,1))
  th.sort()
  return th

def get_tchannels(th=None,tlh=None,trh=None,lh=None,ev=None,raw=False):
  if th==None: th = get_thits(tlh,trh,lh,ev,raw)
  tc = TDCHits()
  while(th!=[]):
    if len(th)>1 and th[1].cid==th[0].cid+1:
      if (abs(th[0].t-th[1].t-diff[th[0].cid-1])<3.25):
        tc.append(TDCHit(2*th[0].cid,(th[0].t+th[1].t)/2.,4,(kminfirst[th[1].cid-1]+kmin[th[0].cid-1])/2.))
        del th[0],th[0]
      else:
        tc.append(TDCHit(2*th[0].cid-1,th[0].t,3,(kminfirst[th[0].cid-1]+kmin[th[0].cid-1])/2.))
        del th[0]
    else:
      tc.append(TDCHit(2*th[0].cid-1,th[0].t,3,(kminfirst[th[0].cid-1]+kmin[th[0].cid-1])/2.))
      del th[0]
  return tc

def get_echannels(eh=None,lh=None,ev=None,raw=False):
  if eh==None: _,_,eh = dispatch_hits(lh,ev,raw)
  ec = TDCHits()
  while(eh!=[]):
    if len(eh)>1 and eh[1].cid==eh[0].cid+1:
      if (abs(eh[0].t-eh[1].t-diff[eh[0].cid+128-1])<8.):
        ec.append(TDCHit(2*eh[0].cid+128,(eh[0].t+eh[1].t)/2.,6,kavg[2*eh[0].cid-1]))
        del eh[0],eh[0]
      else:
        ec.append(TDCHit(2*eh[0].cid-1+128,eh[0].t,5,kavg[2*eh[0].cid-2]))
        del eh[0]
    else:
      ec.append(TDCHit(2*eh[0].cid-1+128,eh[0].t,5,kavg[2*eh[0].cid-2]))
      del eh[0]
  return ec

def get_etchannels(tc=None,ec=None,th=None,tlh=None,trh=None,eh=None,lh=None,ev=None,raw=False):
  if tc==None and ec==None:
    tlh,trh,eh = dispatch_hits(lh,ev,raw)
    th = get_thits(tlh,trh)
    tc = get_tchannels(th)
    ec = get_echannels(eh)
  elif tc==None: tc = get_tchannels(th,tlh,trh,lh,ev,raw)
  elif ec==None: ec = get_echannels(eh,lh,ev,raw)
  etc = TDCHits()
  for he in ec:
    eid = (he.cid-128+1)/2
    if he.cid%2: tchan_max = tchmax[eid-1]
    else: tchan_max = tchmax[eid]
    for ht in tc:
      if tchmin[eid-1]<=ht.cid<=tchan_max:
        if (abs(ht.t-he.t-etdiff[eid])<7.):
          etc.append(TDCHit(he.cid,ht.t,7,he.E))
  etc.sort()
  return etc

      
# functions to get times of trigger or group tdcs
def get_grouptdcs(ev):
  tdcs = [TDCHit(ev.tdc_id[j],ev.tdc_val[j]/10.) for j in range(ev.n_tdc) if ev.tdc_id[j]<52]
  return tdcs

def get_sumtdcs(ev):
  stdcs = split_trgtdcs(ev)
  stdcs2 = TDCHits()
  for x in stdcs:
    if x.cid==55:
      x.t = x.t/10.-833.
      stdcs2+=[x]
  return stdcs2

def trigger_time(ev):
  tdcs = [TDCHit(ev.tdc_id[j],ev.tdc_val[j]) for j in range(ev.n_tdc) if 53<ev.tdc_id[j]<60]
  if tdcs==[]: return None
  stdcs,ntdc = TDCHits(),[]
  for x in tdcs:
    matched = False
    for j,y in enumerate(stdcs):
      if abs(x.t-y.t)<10:
        y.t = (ntdc[j]*y.t+x.t)/(ntdc[j]+1)
        ntdc[j]+=1
        if ntdc[j]==5: return y.t/10.-1000.
        matched = True
        break
    if not matched:
      stdcs.append(TDCHit(x.cid,x.t))
      ntdc.append(1)
  print 'trigger not found',[(x.cid,x.t) for x in tdcs],ntdc

def split_trgtdcs(ev):
  tdcs = [TDCHit(ev.tdc_id[j],ev.tdc_val[j]) for j in range(ev.n_tdc) if 53<ev.tdc_id[j]<60]
  stdcs,ntdc = TDCHits(),[]
  for x in tdcs:
    matched = False
    for j,y in enumerate(stdcs):
      if abs(x.t-y.t)<4:
        y.t = (ntdc[j]*y.t+x.t)/(ntdc[j]+1)
        ntdc[j]+=1
        matched = True
        break
    if not matched:
      stdcs+=[TDCHit(x.cid,x.t)]
      ntdc+=[1]
  for i,x in enumerate(stdcs):
    if ntdc[i]==6: x.cid = 58
  return stdcs

def match_tdcs(tdc1,tdc2,cut=20):
  tdcs = TDCHits()
  for x in tdc1:
    miny = TDCHit(0,100000)
    for y in tdc2:
      if abs(x.t-y.t)<abs(x.t-miny.t): miny = y
    if miny!=TDCHit(0,100000) and abs(x.t-miny.t)<cut:
      miny.t = miny.t-x.t
      tdcs+=[miny]
  return tdcs
