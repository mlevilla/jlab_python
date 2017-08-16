from misc import *

time_start = time.time()

import ROOT
from ROOT import gStyle, gROOT, TGraph, TGraphErrors, TH1F, TH2F, TTree, TFile, TF1, TCanvas, TPaveText, TSpectrum, TLorentzVector, TH1D, TH2D, TChain, TH1, TH2
from array import array

gStyle.SetOptStat(0)
gStyle.SetOptFit(0111)
gStyle.SetPalette(1)
gStyle.SetNumberContours(99)

def list_op(y,dy,ops):
  for i in range(len(dy)):
    if dy[i]==[]: dy[i] = [0 for _ in y[0]]
  a,da = y[0], dy[0]
  for b,db,op in zip(y[1:],dy[1:],ops):
    for i in range(len(a)):
      if op=='+':
        a[i],da[i] = a[i]+b[i],(da[i]**2+db[i]**2)**0.5
      elif op=='-':
        a[i],da[i] = a[i]-b[i],(a[i]**2+b[i]**2)**0.5
      elif op=='*' and a[i]!=0 and b[i]!=0:
        a[i],da[i] = a[i]*b[i],a[i]*b[i]*((da[i]/a[i])**2+(db[i]/b[i])**2)**0.5
      elif op=='/' and a[i]!=0 and b[i]!=0:
        a[i],da[i] = a[i]/b[i],a[i]/b[i]*((da[i]/a[i])**2+(db[i]/b[i])**2)**0.5
  return a,da

def root_to_list(x=[],y=[],dx=[],dy=[]):
  if y==dx==dy==[]: x,y = [],x
  if isinstance(x,TH1): x = [x.GetBinCenter(i+1) for i in range(x.GetNbinsX())]
  if isinstance(x,TGraph): x = [a for a in x.GetX()]
  if isinstance(y,TH1):
    if not x: x = [y.GetBinCenter(i+1) for i in range(y.GetNbinsX())]
    if not dx: dx = [y.GetBinWidth(i+1)/2. for i in range(y.GetNbinsX())]
    if not dy: dy = [y.GetBinError(i+1) for i in range(y.GetNbinsX())]
    y = [y[i+1] for i in range(y.GetNbinsX())]
  elif isinstance(y,TGraph): 
    if not x: x = [a for a in y.GetX()]
    if isinstance(y,TGraphErrors) and not dy: dy = [a for a in y.GetEY()]
    if isinstance(y,TGraphErrors) and not dx: dx = [a for a in y.GetEX()]
    y = [y.Eval(a) for a in x]
  return x,y,dx,dy

def root_op(y,dy,ops):
  if not dy: dy = [[] for _ in y]
  for i in range(len(y)):
    _,y[i],_,dy[i] = root_to_list([],y[i],[],dy[i])
  return list_op(y,dy,ops)
    
def tgraph(x,y=[],dx=[],dy=[],ylim=[],op=[],**args):
  # conversion from TH1F or TGraph
  x,y,dx,dy = root_to_list(x,y,dx,dy)
  # length
  n = len(x)
  if n!=len(y): return TGraph()
  # operation on y axis
  if op:
    _,y2,_,dy2 = root_to_list(x,op[1])
    #print y,y2,dy,dy2
    #print
    y,dy = list_op([y,y2],[dy,dy2],[op[0]])
  # remove points out of range
  if ylim!=[]:
    i0 = 0
    for i in range(n):
      if y[i-i0]<ylim[0] or y[i-i0]>ylim[1]:
        del x[i-i0],y[i-i0]
        if dy: del dy[i-i0]
        if dx: del dx[i-i0]
        i0+=1     
  # build graphs
  n = len(x)
  if n!=len(y): g = TGraph()
  elif not dy: g = TGraph(n,array('d',x),array('d',y))
  elif n==len(dy) and not dx: g = TGraphErrors(n,array('d',x),array('d',y),array('d',[0. for i in range(n)]),array('d',dy))
  elif n==len(dx): g = TGraphErrors(n,array('d',x),array('d',y),array('d',dx),array('d',dy))
  # custom
  custom_obj(g,**args)
  return g


def fitslices(h2,ratio=1,draw=0,axis=0,nmin=1,nmax=-1):
  axes = [h2.GetXaxis(),h2.GetYaxis()]
  nbin = [x.GetNbins() for x in axes]
  low = [x.GetBinLowEdge(1) for x in axes]
  up = [x.GetBinUpEdge(nbin[i]) for i,x in enumerate(axes)]
  title = [x.GetTitle() for x in axes]

  # creating and filling slice histograms
  n = nbin[axis]/ratio
  print n
  hs = [TH1D('hs_'+h2.GetName()+str(i),';'+title[not axis]+';',nbin[not axis],low[not axis],up[not axis]) for i in range(n+2)]
  res = [None for i in range(n+2)]
  means = [0. for i in range(n+2)]
  for i in range(nmin,nmax%(n+2)):
    ns = 0
    for j in range(1,nbin[not axis]+1):
      if axis==0:
        hs[i][j] = sum([h2[k+(nbin[0]+2)*j] for k in range(i*ratio+1,(i+1)*ratio+1)])
      else:
        hs[i][j] = sum([h2[j+(nbin[0]+2)*k] for k in range(i*ratio+1,(i+1)*ratio+1)])
      ns += hs[i][j]
    hs[i].SetEntries(ns)
    if ns!=0 and axis==0:
      means[i] = sum([sum([h2[k+(nbin[0]+2)*j] for j in range(1,nbin[1]+1)])*axes[0].GetBinCenter(k) for k in range(i*ratio+1,(i+1)*ratio+1)])/ns
    elif ns!=0:
      means[i] = sum([sum([h2[j+(nbin[0]+2)*k] for j in range(1,nbin[0]+1)])*axes[1].GetBinCenter(k) for k in range(i*ratio+1,(i+1)*ratio+1)])/ns

  # fitting slice histograms
  i0 = 0
  for i in range(n+2):
    okfit = False
    if hs[i-i0].GetEntries()>h2.GetEntries()/n*0.1 and means[i-i0]!=0:
      blabla = TF1('blabla','gaus+gaus(3)+pol0(6)')
      mean = hs[i-i0].GetBinCenter(hs[i-i0].GetMaximumBin())
      maxi = hs[i-i0].GetMaximum()
      rms = hs[i-i0].GetRMS()
      blabla.SetParLimits(0,0,2*maxi)
      blabla.SetParLimits(1,mean-rms,mean+rms)
      blabla.SetParLimits(2,0.02*rms,2*rms)
      blabla.SetParLimits(3,0,0.25*maxi)
      blabla.SetParLimits(4,mean-2*rms,mean+2*rms)
      blabla.SetParLimits(5,0.05*rms,3*rms)
      blabla.SetParLimits(6,0.,0.1*maxi)
      blabla.SetParameter(1,mean)
      blabla.SetParameter(2,0.1*rms)
      if draw<2: x = hs[i-i0].Fit('blabla','0QR','',mean-5*rms,mean+5*rms)
      else: x = hs[i-i0].Fit('blabla','QR','',mean-5*rms,mean+5*rms)
      if int(x)==0:
        gaus = hs[i-i0].GetFunction('blabla')
        res[i-i0]=[means[i-i0],gaus.GetParameter(1),gaus.GetParameter(2),gaus.GetParError(2)]
      else: 
        del hs[i-i0],means[i-i0],res[i-i0]
        i0+=1
    else: 
      del hs[i-i0],means[i-i0],res[i-i0]
      i0+=1
      
  # returning slice histograms and graph
  if not res: return [hs,None,[None,None]]
  if axis==0:
    g = TGraphErrors(len(res), array('d',[x[0] for x in res]), array('d',[x[1] for x in res]),array('d',[0 for x in res]),array('d',[x[2]*x[3] for x in res]))
  else:
    g = TGraphErrors(len(res), array('d',[x[1] for x in res]), array('d',[x[0] for x in res]),array('d',[x[2]*x[3] for x in res]),array('d',[0 for x in res]))
  if draw>0: g.Fit('pol1','Q')
  else: g.Fit('pol1','0Q')
  par = [g.GetFunction('pol1').GetParameter(i) for i in range(2)]
  return [hs,g,par]


def tree_init(name,l,title='',nmax=100,fout=None):
  t = TTree(name,title)
  if fout!=None: fout.SetDirectory(fout)
  tmp = {'i':[0],'f':[0.],'i[':[0]*nmax,'f[':[0.]*nmax}
  var = {}
  for x,y in l:
    var[x] = array(y[0],tmp[y[:2]])
    t.Branch(x,var[x],x+y[1:]+'/'+y[0].upper())
  del tmp
  return t,var

def add_histo(h,hplus):
  n = h.GetNbinsX()
  if hplus.GetNbinsX()!=n: 
    print "not the same number of bins"
    return
  for i in range(n+2): 
    if hplus[i]!=hplus[i]: continue
    h[i] += hplus[i]
  h.SetEntries(h.GetEntries()+hplus.GetEntries())
  

def custom_obj(h,**args):
  # default
  if 'titlesize' not in args: args['titlesize']=[0.07,0.07]
  if 'labelize' not in args: args['labelsize']=[0.05,0.05]
  if 'offset' not in args: args['offset']=[0.85,0.85]
  if 'style' not in args: args['style']=[1,0,20]
  
  # color 
  color = None
  if 'color' in args: 
    color = args['color'] 
    if isinstance(color,int): color = [color,color,color]
  elif h.GetMarkerColor()==ROOT.kBlack: color = [ROOT.kBlack,ROOT.kBlack,ROOT.kRed]
  if color is not None:
    h.SetLineColor(color[0])
    h.SetFillColor(color[1])
    h.SetMarkerColor(color[2])

  # style
  if 'style' in args:
    style = args['style']
    h.SetLineStyle(style[0])
    h.SetFillStyle(style[1])
    h.SetMarkerStyle(style[2])
  elif isinstance(h,TH1): h.SetFillStyle(0)
    
  if 'width' in args: h.SetLineWidth(args['width'])
  if 'stat' in args: h.SetStats(args['stat'])

  if 'title' in args: h.SetTitle(args['title'])
  if 'name' in args: h.SetName(args['name'])

  # axis
  a = [xa,ya] = [h.GetXaxis(),h.GetYaxis()]
  for i,x in enumerate(a):
    if 'titlesize' in args: x.SetTitleSize(args['titlesize'][i])
    if 'labelsize' in args: x.SetLabelSize(args['labelsize'][i])
    x.SetLabelOffset(0.01)
    x.SetNdivisions(505)
    if 'offset' in args: x.SetTitleOffset(args['offset'][i])
  if 'xrg' in args: xa.SetRangeUser(args['xrg'][0],args['xrg'][1])
  if 'yrg' in args: ya.SetRangeUser(args['yrg'][0],args['yrg'][1])
 

def custom(h,first=True,c=None,log=0,grid=False,legend=[],insert=[],margin=[0.13,0.13],opt='',**objargs):
  if c == None: c = TCanvas('c')
  custom_obj(h,**objargs)

  if first:
    c.SetBottomMargin(margin[0])
    c.SetLeftMargin(margin[1])
    c.SetTopMargin(0.2-margin[0])
    c.SetRightMargin(0.2-margin[1])
    # grid and log
    if grid:
      c.SetGridx()
      c.SetGridy()
    if log&1: c.SetLogx()
    if log&2: c.SetLogy()
    if log&4: c.SetLogz()

  # draw graph or histogram
  c.cd()
  if isinstance(h,TH1):  
    if first: h.Draw(opt)
    else: h.Draw('same '+opt)
  elif isinstance(h,TH2): h.Draw('colz')
  else: 
    if opt=='': opt = 'p'
    if first: h.Draw('a'+opt)
    else: h.Draw(opt)
      

  # draw legend
  if legend: c.BuildLegend(legend[0],legend[1],legend[2],legend[3],'')
  
  # draw insert
  t = []
  for x in insert:
    y = TPaveText(x[0],x[1],x[2],x[3],'NDC')
    y.InsertText(x[4])
    y.SetBorderSize(0)
    y.SetFillStyle(0)
    y.SetTextFont(42)
    y.SetTextSize(0.05)
    y.SetTextColor(x[5])
    y.Draw()
    t.append(y)
  c.Update()
  return c,t

def fit_rgaus(h,threshold=2.,sigma=3.,xrg=None):
  if xrg!=None: h.GetXaxis().SetRangeUser(xrg[0],xrg[1])
  peak = h.GetBinCenter(h.GetMaximumBin())
  std = h.GetStdDev()
  integ = h.Integral()
  height = h.GetMaximum()
  if integ<=threshold: return [height,peak,std,0]
  h.Fit('gaus','q0r','',peak-sigma*std,peak+sigma*std)
  f = h.GetFunction('gaus')
  if f==None: return [height,peak,std,0]
  return [f.GetParameter(i) for i in range(3)]+[1]

  
def fit_gaus(h,threshold=2.,sigma=3.,xrg=None):
  [_,mean,sig,flag] = fit_rgaus(h,threshold,sigma,xrg)
  if not flag: return [mean,sig,0.,0.,0]
  h.Fit('gaus','qr','',mean-sigma*sig,mean+sigma*sig)
  f = h.GetFunction('gaus')
  mean = f.GetParameter(1)
  sig = f.GetParameter(2)
  h.GetXaxis().SetRangeUser(mean-sigma*sig,mean+sigma*sig)
  integ = h.Integral()
  if integ<=threshold: return [mean,sig,0.,0.,0]
  if xrg!=None: h.GetXaxis().SetRangeUser(xrg[0],xrg[1])
  return [f.GetParameter(i) for i in range(1,3)]+[f.GetParError(i) for i in range(1,3)]+[1]

def fit_2gaus(h,threshold=2.,xrg=None):
  [height,mean,sig,flag] = fit_rgaus(h,threshold,3,xrg)
  if not flag: return [mean,sig,0.,0.,0]
  ffit = TF1('2gaus','gaus(0)+gaus(3)+[6]',h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
  ffit.SetParameters(height,mean,sig,height*0.1,mean,4*sig,0.,0.)
  ffit.SetParLimits(0,height*0.8,height*1.2)
  ffit.SetParLimits(1,mean-sig,mean+sig)
  ffit.SetParLimits(2,0.8*sig,1.2*sig)
  ffit.SetParLimits(3,0,height*0.3)
  ffit.SetParLimits(4,mean-3*sig,mean+3*sig)
  ffit.SetParLimits(5,sig,6*sig)
  ffit.SetParLimits(6,0,0.1*height)
  if xrg is None: h.Fit(ffit,'q')
  else: h.Fit(ffit,'rq','',xrg[0],xrg[1])
  return [ffit.GetParameter(1),ffit.GetParameter(2),ffit.GetParError(1),ffit.GetParError(2),1]
  
  
def fit_list(f,x,y,dx=[],dy=[],xrg=[],sel=[],start=None):
  if not sel: g,g2 = tgraph(x,y,dx,dy),TGraph()
  else: 
    g = tgraph([x[i] for i in sel],[y[i] for i in sel],[dx[i] for i in sel] if dx!=[] else [],[dy[i] for i in sel] if dy!=[] else [])
    notsel = list(set(range(len(x)))-set(sel))
    g2 = tgraph([x[i] for i in notsel],[y[i] for i in notsel],[dx[i] for i in notsel] if dx!=[] else [],[dy[i] for i in notsel] if dy!=[] else [],color=[ROOT.kBlack,ROOT.kBlack,ROOT.kBlue])
  ffit = TF1('ffit',f,min(x),max(x))
  if start is not None:
    for i in range(len(start)): ffit.SetParameter(i,start[i])
  if not xrg: g.Fit(ffit,'q')
  else: g.Fit(ffit,'rq','',xrg[0],xrg[1])
  l = [ffit.GetParameter(i) for i in range(f.count('['))]+[ffit.GetParError(i) for i in range(f.count('['))]
  return [g,g2]+l  


def writelisto(l,f=None,folder=[],**custom_arg):
  if l==None: return 
  if not isinstance(l,list) and not isinstance(l,dict):
    custom_obj(l,**custom_arg)
    l.Write() 
  else:
    lnames = []
    if isinstance(folder,str): folder = [folder]
    if isinstance(l,dict): l = l.values()
    for i,x in enumerate(l): 
      lnames = [] if f is None else [y.GetName() for y in f.GetListOfKeys()]
      if not folder: u = f
      elif f is not None:
        if not isinstance(folder[0],list): fold = str(folder[0])
        elif len(folder[0])!=len(l): fold = str(folder[0][0])
        else: fold = str(folder[0][i])
        if fold in lnames: u = f.Get(fold)
        else: u = f.mkdir(fold)
      if u is not None: u.cd()
      writelisto(x,u,folder[1:],**custom_arg)


def listofhisto(name,suffix,nbinx,xmin=None,xmax=None,nbiny=None,ymin=None,ymax=None,xbin=None,ybin=None,title='',fout=None):
  if fout!=None: fout.cd()
  if len(suffix)==0:
    if nbiny is None:
      if xbin is None: h = TH1D(name,title,nbinx,xmin,xmax)
      else: h = TH1D(name,title,nbinx,array('d',xbin))
    else:
      if xbin is None: h = TH2D(name,title,nbinx,xmin,xmax,nbiny,ymin,ymax)
      else: h = TH2D(name,title,nbinx,xbin,nbiny,ybin)
      h.SetOption('colz')
    return h
  else:
    if not isinstance(suffix[0],list): suffix = [suffix]
    suf = suffix[0]
    if not isinstance(nbinx,list): nbinx = [nbinx for y in suf]
    if not isinstance(xmin,list): xmin = [xmin for y in suf]
    if not isinstance(xmax,list): xmax = [xmax for y in suf]
    if xbin is None or not isinstance(xbin[0],list): xbin = [xbin for y in suf]
    if not isinstance(nbiny,list): nbiny = [nbiny for y in suf]
    if not isinstance(ymin,list): ymin = [ymin for y in suf]
    if not isinstance(ymax,list): ymax = [ymax for y in suf]
    if ybin is None or not isinstance(ybin[0],list): ybin = [ybin for y in suf]
    if not isinstance(title,list): title = [title for y in suf]
    return [listofhisto(name+'_'+str(suf[i]),suffix[1:],nbinx[i],xmin[i],xmax[i],nbiny[i],ymin[i],ymax[i],xbin[i],ybin[i],title[i]) for i in range(len(suf))]


def fill_pulls(p,h,hi,const=True,factor=[1.,1.],w=False,xrg=[]):
  n = h.GetNbinsX()
  if not isinstance(p,list): p = [p]
  m = len(p)
  l = []
  for i in range(1,n+1):
    if xrg and (xrg[0]>h.GetBinLowEdge(i) or xrg[1]<h.GetBinLowEdge(i+1)): continue
    s2 = (hi.GetBinError(i)*factor[1])**2-(h.GetBinError(i)*factor[0])**2 if const else (hi.GetBinError(i)*fact[1])**2
    if s2<=0: continue
    p[(i-1)/(n/m)].Fill((hi[i]*factor[1]-h[i]*factor[0])/s2**0.5)
  

def filetopdfs(f,folder,opt={}):
  if isinstance(f,str): f = TFile(f)
  l = f.GetListOfKeys()
  for i,x in enumerate(l):
    y = x.ReadObj()
    if x.IsFolder(): 
      opt2 = opt[i] if isinstance(opt,list) else opt
      if not os.access(folder+'/'+x.GetName(),os.R_OK): os.mkdir(folder+'/'+x.GetName())
      filetopdfs(y,folder+'/'+x.GetName(),opt2)
    else:
      (c,_) = custom(y,**opt)
      c.Print(folder+'/'+x.GetName()+'.pdf')

def getlisto(f):
  if f.IsFolder(): return [getlisto(x.ReadObj()) for x in  f.GetListOfKeys()]
  else: return f

def get_binning(h):
  return [h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)]

def get_chain(name,filename,n,ext='.root',start=0):
  c = TChain(name)
  for i in range(start,n+start): c.AddFile(filename+'_'+str(i)+ext)
  return c

def divide_binwidth(h):
  for i in range(h.GetNbinsX()):
    err = h.GetBinError(i)
    h[i] = h[i]/h.GetBinWidth(i)
    h.SetBinError(i,err/h.GetBinWidth(i))

def multiply_binwidth(h):
  for i in range(h.GetNbinsX()):
    err = h.GetBinError(i)
    h[i] = h[i]*h.GetBinWidth(i)
    h.SetBinError(i,err*h.GetBinWidth(i))

def hintegral(h,a,b):
  binning = get_binning(h)
  ia = max(i for i,x in enumerate(binning) if a>=x)
  ib = min(i for i,x in enumerate(binning) if b<=x)
  return sum(h[i+1] for i in range(ia,ib))
