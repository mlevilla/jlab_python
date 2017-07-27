#!/apps/python/python-2.7.1/bin/python

from params import *
from misc_root import *
from math import exp,erf
from ROOT import TF1Convolution

largs = [('i','','','input_directory'),('p','','','pass_number'),('m','fsnake','','clustering_method'),('b',False,'','batch_mode'),('s',0,'','show_bar')]
print_help(largs)
[indir,ipass,method,batch,show] = [catch_arg(x,y,z) for x,y,z,c in largs]

if ipass=='' and '/pass' in indir: 
  ind = indir.index('/pass')
  ipass = indir[ind+5:]
if indir=='': indir = os.environ['WORK']+'/calib_results/'+method+'/pass'+ipass

if batch: 
  jsub(project='prad',track='analysis',jobname='spatial_'+method,command=os.path.abspath('spatial_reso.py'),options=' '.join(sys.argv[1:]).replace('-b',''),memory='2 GB')
  sys.exit()

fout = TFile(indir+'/spatial_resolution_42.root','recreate')

#######################
## data distribution ##
#######################

runs =  [955,956,957,958,960,961,962,965,966,967,968,969,970,971,972]
ypos_ref = [[-113.2],[-133.6],[-154.0],[-174.7],[-195.6],[-216.4],[-235.0],[-255.6],[-276.4],[-297.1],[-317.9],[-338.6],[-368.1,-406.2],[-444.9,-483.1],[-521.3,-559.5]]
px = [[14.296, 5.748], [14.135, 5.688], [14.032, 5.617], [13.783, 5.623], [13.686, 5.586], [13.683, 5.504], [13.542, 5.523], [13.174, 5.543]]
py = [[14.303, 5.749], [14.127, 5.672], [13.988, 5.609], [13.809, 5.65], [13.725, 5.562], [13.6, 5.523], [13.538, 5.457], [13.393, 5.487]]

nbin0,r0 = 200,40
h2x = [[TH1F('h2x'+str(k)+'_'+str(j),'',nbin0,-r0,r0) for j in range(8)] for k in range(90)]
ftest = TF1('ftest','[0]*exp(-0.5*(x-[1])^2/[2]^2)*(1+erf([3]*(x-[1])))',-40,40)
ftest.SetParameters(100,0,10,0)
ffit = TF1('ffit','[0]+[1]/sqrt(x)+[2]/x')

for run,yref in zip(runs,ypos_ref):
  for i in range(len(yref)):
    print run,i
    f1 = TFile(indir+'/trees/tree_'+str(run)+'.root')
    t1 = f1.hycal
    lcor = readfile('/work/hallb/prad/mlevilla/transporter_time_position/transporter_time_position_'+str(run)+'.txt',[long,float,float])
    fout.cd()
    h1x = [TH1F('h1x'+'_'+str(run)+str(j),'',nbin0,-r0,r0) for j in range(8)]
    h1y = [TH1F('h1y'+'_'+str(run)+str(j),'',nbin0,-r0,r0) for j in range(8)]
    for ev in progress(t1,n=t1.GetEntries(),show=show):
      if ev.n_cl!=1: continue
      if ev.n_tag!=1: continue
      if ev.trigger not in [1,2]: continue
      if {'calib':ev.status[0]%2==0,'fsnake':False}[method]: continue
      if abs(ev.tg[0]-ev.tcl[0])>40: continue
      if ev.nh[0]<4: continue
      if ev.Eg[0]>=1000 or ev.Eg[0]<200: continue
      if abs(1-ev.E[0]/ev.Eg[0])>0.1: continue
      if abs(ev.x[0])>330: continue
      if ev.y[0]>-103.75 or ev.y[0]<=-581.65: continue
      [xpos,ypos] = correct_w_time(lcor,ev.time)
      if abs(ypos-yref[i])>0.1: continue
      _ = h1x[int((ev.Eg[0]-200)/100)].Fill(ev.x[0]-xpos)
      _ = h1y[int((ev.Eg[0]-200)/100)].Fill(ev.y[0]-ypos)
      if ev.y[0]>-352.75: j = int(-(ev.y[0]+103.75)/20.75*5.)
      else: j = int(-(ev.y[0]+352.75)/38.15*5.)+60
      _ = h2x[j][int((ev.Eg[0]-200)/100)].Fill(ev.x[0]-xpos)

    sigma_x, sigma_y, err_x, err_y = [],[],[],[]
    
    for j in range(8):
      ftest.SetParameter(3,h1x[j].GetSkewness())
      h1x[j].Fit(ftest,'','',-15,15)
      fconvx = TF1Convolution('(1+erf('+str(ftest.GetParameter(3))+'*x))*exp(-0.5*x^2/'+str(px[j][0])+'^2)/(1+x^2/'+str(px[j][1])+'^2)','[0]*exp(-0.5*(x-[1])^2/[2]^2)',-40,40)
      ffitx = TF1('ffitx',fconvx,-40,40,3)
      ffitx.SetParLimits(1,-10,10)
      ffitx.SetParLimits(2,1,10)
      ffitx.SetParameters(100,ftest.GetParameter(1),3.)
      h1x[j].Fit(ffitx,'LR','',-15,15)
      sigma_x.append(ffitx.GetParameter(2))
      err_x.append(ffitx.GetParError(2))

    # for j in range(8):
    #   ftest.SetParameter(3,h1y[j].GetSkewness())
    #   h1y[j].Fit(ftest,'','',-15,15)
    #   fconvy = TF1Convolution('(1+erf('+str(ftest.GetParameter(3))+'*x))*exp(-0.5*x^2/'+str(py[j][0])+'^2)/(1+x^2/'+str(py[j][1])+'^2)','[0]*exp(-0.5*(x-[1])^2/[2]^2)',-40,40)
    #   ffity = TF1('ffity',fconvy,-40,40,3)
    #   ffity.SetParameters(100,ftest.GetParameter(1),3.)
    #   h1y[j].Fit(ffity,'LR','',-15,15)
    #   sigma_y.append(ffity.GetParameter(2))
    #   err_y.append(ffity.GetParError(2))
      
    gx = tgraph([0.25+0.1*j for j in range(8)],sigma_x,dy=err_x)
    gx.Fit(ffit,'Q')

    # gy = tgraph([0.25+0.1*j for j in range(8)],sigma_y,dy=err_y)
    # gy.Fit(ffit,'Q')
    
    # sigma_r = [(0.5*(x**2+y**2))**0.5 for x,y in zip(sigma_x,sigma_y)]
    # err_r = [(0.5*(x**2+y**2))**0.5 for x,y in zip(err_x,err_y)]
    # gr = tgraph([0.25+0.1*j for j in range(8)],sigma_r,dy=err_r)
    # gr.Fit(ffit,'Q')

    u = fout.mkdir(str(run)+'_'+str(i))
    u.cd()
    for j in range(8):
      c,t = custom(h1x[j],first=True,xtitle='x_{cluster}-x_{transporter} (mm)')
      h1x[j].Write()
    custom(gx,c=c,first=True,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{x} (mm)')
    gx.Write()

    # for j in range(8): 
    #   c,t = custom(h1y[j],c=c,first=True,xtitle='y_{cluster}-y_{transporter} (mm)')
    #   h1y[j].Write()
    # custom(gy,c=c,first=True,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{y} (mm)')
    # gy.Write()

    # custom(gr,c=c,first=True,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{r} (mm)')
    # gr.SetName('gr_'+str(run)+'_'+str(i))
    # gr.Write()       
    

    f1.Close()

sigma_x, err_x = [],[]
for k in range(90):
  sigma_x2, err_x2, = [],[]
  for j in range(8):
    ftest.SetParameter(3,h2x[k][j].GetSkewness())
    h1x[j].Fit(ftest,'','',-15,15)
    fconvx = TF1Convolution('(1+erf('+str(ftest.GetParameter(3))+'*x))*exp(-0.5*x^2/'+str(px[j][0])+'^2)/(1+x^2/'+str(px[j][1])+'^2)','[0]*exp(-0.5*(x-[1])^2/[2]^2)',-40,40)
    ffitx = TF1('ffitx',fconvx,-40,40,3)
    ffitx.SetParLimits(1,-10,10)
    ffitx.SetParLimits(2,1,10)
    ffitx.SetParameters(100,ftest.GetParameter(1),3.)
    h2x[k][j].Fit(ffitx,'LR','',-15,15)
    sigma_x2.append(ffitx.GetParameter(2))
    err_x2.append(ffitx.GetParError(2))
  sigma_x.append(sigma_x2)
  err_x.append(err_x2)

gx = [tgraph([0.25+0.1*j for j in range(8)],sigma_x[k],dy=err_x[k]) for k in range(90)]
lx = []
for k in range(90):
  gx[k].Fit(ffit,'Q')
  lx.append(ffit.Eval(1.))

gtrans = tgraph([103.75+20.75/5.*(j+0.5) for j in range(60)]+[352.75+38.15/5.*(j+0.5) for j in range(30)],lx)

u = fout.mkdir('all')
u.cd()
for k in range(90):
  v = u.mkdir(str(k))
  v.cd()
  custom(gx[k],c=c,first=True,xtitle='E_{#gamma} (GeV)',ytitle='#sigma_{x} (mm)')
  gx[k].SetName('gx_'+str(i))
  gx[k].Write()
  for j in range(8):
    custom(h2x[k][j],c=c,first=True,xtitle='x_{cluster}-x_{transporter} (mm)')
    h2x[k][j].Write()
u.cd()
custom(gtrans,c=c,first=True,xtitle='y_{hycal} (mm)',ytitle='#sigma_{x} (mm)')
gtrans.SetName('gtrans')
gtrans.Write()

fout.Close()
    
