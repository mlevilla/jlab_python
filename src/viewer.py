from params import *
from ROOT import TCanvas,TPad,TH2F,TGraph,gStyle,TGaxis,TPaletteAxis

gStyle.SetPalette(1)
gStyle.SetNumberContours(99) 
gStyle.SetFrameLineColor(gStyle.GetCanvasColor());
TGaxis.SetMaxDigits(3)

cell_size = [[2.077,2.075],[3.815,3.815]]
offset = 0.05
size2_x = (1-4*offset)/(34*cell_size[0][0]+2*6*cell_size[1][0])*6*cell_size[1][0]
size2_y = (1-2*offset)/(34*cell_size[0][1]+2*6*cell_size[1][1])*6*cell_size[1][1]
size3_x = (1-4*offset)/(34*cell_size[0][0]+2*6*cell_size[1][0])*24*cell_size[1][0]
size3_y = (1-2*offset)/(34*cell_size[0][1]+2*6*cell_size[1][1])*24*cell_size[1][1]

load_names()

class hcviewer:
  def __init__(self,name='',log=False,zmin=10.,zmax=1100.):
    self.name = name
    self.log = log
    self.c = None
    self.p = [None]*5
    self.g_hits = [[TGraph(1) for j in range(5)] for i in range(3)]
    self.g_cls = [[TGraph(1) for j in range(5)] for i in range(3)]
    self.zmin = zmin
    self.zmax = zmax
    self.hp = TH2F('hpalette_'+self.name,'',1,0,1,1,0,1)
    self.hp.SetEntries(1)
    self.hp.SetStats(0)
    self.hp.GetXaxis().SetAxisColor(0)
    self.hp.GetXaxis().SetLabelSize(0)
    self.hp.GetXaxis().SetTickLength(0)
    self.hp.GetYaxis().SetAxisColor(0)
    self.hp.GetYaxis().SetLabelSize(0)
    self.hp.GetYaxis().SetTickLength(0)
    self.palette = TPaletteAxis()
    for i in range(3):
      for j in range(5):
        self.g_hits[i][j].SetMarkerStyle([24,26,27][i])
        self.g_hits[i][j].SetPoint(0,-1,-1)
        self.g_cls[i][j].SetMarkerStyle([24,26,27][i])
        self.g_cls[i][j].SetPoint(0,-1,-1)
        self.g_cls[i][j].SetMarkerColor(632)
    self.init_histos()

  def init_histos(self):
    hcrys = TH2F('hcrys_'+self.name,'',34,0,34,34,0,34)
    hbot = TH2F('hbot_'+self.name,'',24,0,24,6,0,6)
    htop = TH2F('htop_'+self.name,'',24,0,24,6,0,6)
    hleft = TH2F('hleft_'+self.name,'',6,0,6,24,0,24)
    hright = TH2F('hright_'+self.name,'',6,0,6,24,0,24)
    self.h = [hcrys,hbot,htop,hleft,hright]
    self.h[0].SetNdivisions(34,'X')
    self.h[0].SetNdivisions(34,'Y')
    for i in range(1,3):
      self.h[i].SetNdivisions(24,'X')
      self.h[i].SetNdivisions(6,'Y')
    for i in range(3,5):
      self.h[i].SetNdivisions(24,'Y')
      self.h[i].SetNdivisions(6,'X')
    for x in self.h:
      x.SetStats(0)
      x.SetTickLength(0,'X')
      x.SetTickLength(0,'Y')

  def init_canvas(self,c=None):
    if c==None: self.c = TCanvas('c_'+self.name,'',770,700)
    else: self.c = c
    if self.log: self.c.SetLogz()
    self.set_range(self.zmin,self.zmax)
    pcrys = TPad('pcrys_'+self.name,'',offset+size2_x,offset+size2_y,1-3*offset-size2_x,1-offset-size2_y)
    pbot = TPad('pbot_'+self.name,'',offset+size2_x,offset,offset+size2_x+size3_x,offset+size2_y)
    ptop = TPad('ptop_'+self.name,'',1-3*offset-size2_x-size3_x,1-offset-size2_y,1-3*offset-size2_x,1-offset)
    pleft = TPad('pleft_'+self.name,'',offset,1-offset-size2_y-size3_y,offset+size2_x,1-offset-size2_y)
    pright = TPad('pright_'+self.name,'',1-3*offset-size2_x,offset+size2_y,1-3*offset,offset+size2_y+size3_y)
    self.p = [pcrys,pbot,ptop,pleft,pright]
    for x in self.p:
      x.SetMargin(0,0,0,0)
      x.SetGrid()
      if self.log: x.SetLogz()

  def draw(self,c=None,cluster=True):
    self.init_canvas(c)
    for x,y in zip(self.p,self.h):
      x.cd()
      y.Draw('col')
      self.c.cd()
      x.Draw()
    if cluster:
      for g in self.g_hits:
        for x,y in zip(self.p,g):
          x.cd()
          y.Draw('p')
          self.c.cd()
          x.Draw()
      for g in self.g_cls:
        for x,y in zip(self.p,g):
          x.cd()
          y.Draw('p')
          self.c.cd()
          x.Draw()
    
  def update(self): 
    if self.c!=None: self.c.Update()

  def reset(self):
    for x in self.h: x.Reset()
    for x in self.g_hits+self.g_cls:
      for y in x:
        y.Set(1)
        y.SetPoint(0,-1,-1)
    
  def __getitem__(self,index):
    ih,xbin,ybin,ibin = self.binning(index)
    return self.h[ih][ibin]

  def __setitem__(self,index,val):
    ih,xbin,ybin,ibin = self.binning(index)
    self.h[ih].Fill(xbin,ybin,val)

  def __delitem__(self,index):
    ih,xbin,ybin,ibin = self.binning(index)
    self.h[ih][ibin] = 0.
    self.h[ih].SetEntries(self.h[ih].GetEntries()-1)

  def binning(self,index):
    if type(index)==int: index=module_names[index]
    mid = int(index[1:])-1
    if index[0]=='W': return (0,mid%34,33-mid/34,mid%34+1+(34-mid/34)*36)
    elif mid>719 and mid%30>5: return (1,mid%30-6,29-mid/30,mid%30-6+1+(30-mid/30)*26)
    elif mid<180 and mid%30<24: return (2,mid%30,5-mid/30,mid%30+1+(6-mid/30)*26)
    elif mid%30<6: return (3,mid%30,23-mid/30+6,mid%30+1+(24-mid/30+6)*8)
    elif mid%30>23: return (4,mid%30-24,23-mid/30,mid%30-24+1+(24-mid/30)*8)
    else:
      print 'hit outside of hycal',index,mid
      return (0,0,0,0)

  def position(self,x,y):
    if abs(x)<trans[0] and abs(y)<trans[1]:
      return (0,(trans[0]+x)/cell_size[0][0],(trans[1]+y)/cell_size[0][1])
    elif x>-trans[0] and y<-trans[1]:
      return (1,(trans[0]+x)/cell_size[1][0],(trans[1]+y)/cell_size[1][1]+6.)
    elif x<trans[0] and y>trans[1]:
      return (2,(-trans[0]+x)/cell_size[1][0]+24.,(-trans[1]+y)/cell_size[1][1]-6.)
    elif x<-trans[0] and y<trans[1]:
      return (3,(trans[0]+x)/cell_size[1][0]+6.,(-trans[1]+y)/cell_size[1][1]+24)
    elif x>trans[0] and y>-trans[1]:
      return (4,(-trans[0]+x)/cell_size[1][0],(trans[1]+y)/cell_size[1][1])
    else:
      print 'point is nowher on hycal'
      return (0,0,0)

  def addcluster(self,l,li,i,cl=None):
    j = [0]*5
    for k in range(5): self.g_hits[i][k].Set(len(l))
    for x,y in zip(l,li):
      ih,xbin,ybin,ibin = self.binning(y)
      self.g_hits[i][ih].SetPoint(j[ih],xbin+0.5,ybin+0.5)
      j[ih]+=1
    for k in range(5):
      if j[k]>0: self.g_hits[i][k].Set(j[k])
      else:
        self.g_hits[i][k].Set(1)
        self.g_hits[i][k].SetPoint(0,-1,-1)
    if cl!=None: 
      ih,x,y = self.position(cl[0],cl[1])
      self.g_cls[i][ih].SetPoint(0,x,y)

  def fill(self,l,li=[range(1728)],reset=False,cl=None):
    if reset: self.reset()
    for i,(x,y) in enumerate(zip(l,li)):
      if i<3 and cl!=None and i<len(cl): self.addcluster(x,y,i,cl[i])
      elif i<3 and cl!=None: self.addcluster(x,y,i)
      for a,b in zip(x,y): self[b] = a
    for x in self.p: 
      if x!=None: x.Modified()

  def set_range(self,zmin,zmax):
    self.zmin = zmin
    self.zmax = zmax
    for x in self.h: x.GetZaxis().SetRangeUser(self.zmin,self.zmax)
    self.hp.GetZaxis().SetRangeUser(self.zmin,self.zmax)
    self.hp[4] = 0. if self.log else self.zmin-1
    self.hp.Draw('colz')
    self.update()
    self.palette = self.hp.GetListOfFunctions().FindObject("palette");
    self.palette.SetX1NDC(0.88)
    self.palette.SetX2NDC(0.93)
    self.hp.Draw('colz')

  def set_log(self,log=True):
    self.log = log
    if log==True:
      for x in self.p: 
        if x!=None: x.SetLogz()
    self.hp[4] = 0. if self.log else self.zmin-1

    
