#!/apps/anaconda/anaconda-2.0.1/bin/python

from params import *

indir = catch_arg('i','')
gaindir_in = catch_arg('g','')
ipass = catch_arg('p','')
method = catch_arg('m','fsnake')
write_gain = catch_arg('w',False)
show = catch_arg('s',0)
name = catch_arg('n','')
run = catch_arg('r','')
batch = catch_arg('b',False)
histo = catch_arg('se',False)
hadd = catch_arg('hadd',False)
final = catch_arg('f',False)
[var,cut] = catch_arg('c',['',0.],'several')

if ipass=='' and indir=='': sys.exit('not enough arguments: '+sys.argv)

###################################################################################################
################ batch jobs #######################################################################
###################################################################################################
if batch: 
  jsub(project='prad',track='analysis',jobname='histo'+name+'_'+run,command=os.path.abspath('snake.py'),options=' '.join(sys.argv[1:]).replace('-b',''),memory='2 GB')
  sys.exit()

if ipass=='' and '/pass' in indir: 
  ind = indir.index('/pass')
  ipass = indir[ind+5:indir.index('/',ind+2)]
if indir=='': indir = os.environ['WORK']+'/calib_results/'+method+'/pass'+ipass
if gaindir_in=='': gaindir_in = indir.replace('/calib_results/','/calib_files/')

###################################################################################################
################ merging histograms ###############################################################
###################################################################################################

if hadd:
  jsub(project='prad',track='analysis',jobname='hadd'+name,command='hadd -f '+indir+'/histos'+name+'.root '+indir+'/histos'+name+'_*.root',memory='4 GB')
  sys.exit()

from misc_root import *
from ROOT import *

print ' '.join(sys.argv)
print
 
mlist = module_names

###################################################################################################
################ creating histograms from trees ###################################################
###################################################################################################
if histo:
  print run
  #mean_gain = readfile(indir+'/mean_0.txt',[float]*9,cols=range(1,10))
  #sigma_gain = readfile(indir+'/sigma_0.txt',[float]*9,cols=range(1,10))
  fout = TFile(indir+'/histos'+name+'_'+run+'.root','recreate')
  hr = [[TH1F('hr_'+x+'_'+str(j),';E_{cluster}/E_{#gamma};',1000,0.,5.) for j in range(9)] for x in mlist]
  he = [[TH1F('he_'+x+'_'+str(j),';E_{cluster} (MeV);',1200,0.,1200.)  for j in [1,2,5]] for x in mlist]
  heg = [[TH1F('heg_'+x+'_'+str(j),';E_{#gamma} (MeV);',1200,0.,1200.)  for j in [1,2,5]] for x in mlist]
  ht = [[TH1F('htdiff_'+x+'_'+str(j),';t_{hycal}-t_{tagger} (ns);',1000,-200.,200.)  for j in [1,2,5]] for x in mlist]
  tr = {1:0,2:1,5:2}
  f = TFile(indir+'/tree_'+run+'.root')
  t = f.hycal
  n = t.GetEntries()
  for (i,x) in enumerate(progress(t,n=n,show=show,modulo=10000)):
    if t.n_cl!=1: continue
    if t.trigger not in [1,2,5]: continue
    j = int((t.E[0]-200)/100.)
    if j>8: continue
    if t.Eg==0: continue
    b = False
    for v,c in zip(var,cut):
      if v=='e': b = b or (t.E[0]<c)
      #if v=='sigma': b = b or t.Eg==0. or abs(t.E[0]/t.Eg-mean_gain[t.id[0]][j])>c*sigma_gain[t.id[0]][j]
      if v=='n': b = b or t.nh[0]<c
      if v=='chi2': b = b or t.chi2[0]>c
      if v=='x': b = b or abs(t.xpos-t.x[0])>c
      if v=='y': b = b or abs(t.ypos-t.y[0])>c
      if v=='t': b = b or abs(t.tg-t.tcl[0])>c
      if v=='tr': b = b or t.trigger>c
    if b: continue
    he[t.id[0]][tr[t.trigger]].Fill(t.E[0])
    heg[t.id[0]][tr[t.trigger]].Fill(t.Eg)
    ht[t.id[0]][tr[t.trigger]].Fill(t.tg-t.tcl[0])
    hr[t.id[0]][j].Fill(t.E[0]/t.Eg)

    k = 1728 if t.id2[0]<1001 else 1729
    he[k][tr[t.trigger]].Fill(t.E[0])
    heg[k][tr[t.trigger]].Fill(t.Eg)
    ht[k][tr[t.trigger]].Fill(t.tg-t.tcl[0])
    hr[k][j].Fill(t.E[0]/t.Eg)

    if module_names[t.id[0]] in names_trans_lg:
      he[1730][tr[t.trigger]].Fill(t.E[0])
      heg[1730][tr[t.trigger]].Fill(t.Eg)
      ht[1730][tr[t.trigger]].Fill(t.tg-t.tcl[0])
      hr[1730][j].Fill(t.E[0]/t.Eg)

    if module_names[t.id[0]] in names_trans_pwo:
      he[1731][tr[t.trigger]].Fill(t.E[0])
      heg[1731][tr[t.trigger]].Fill(t.Eg)
      ht[1731][tr[t.trigger]].Fill(t.tg-t.tcl[0])
      hr[1731][j].Fill(t.E[0]/t.Eg)

    if abs(t.x[0])<2*20.77 and abs(t.y[0])<2*20.75:
      he[1732][tr[t.trigger]].Fill(t.E[0])
      heg[1732][tr[t.trigger]].Fill(t.Eg)
      ht[1732][tr[t.trigger]].Fill(t.tg-t.tcl[0])
      hr[1732][j].Fill(t.E[0]/t.Eg)

  # writing histograms
  fout.cd()
  for i in range(1733): 
    u = fout.mkdir(mlist[i])
    u.cd()
    for j in range(9): hr[i][j].Write()
    for j in range(3): 
      he[i][j].Write()
      heg[i][j].Write()
      ht[i][j].Write()
  fout.Close()

###################################################################################################
################ fit of histograms and output #####################################################
###################################################################################################
if final:
  gains0 = readfile(gaindir_in+'/calib_889-979.txt',[float],cols=[1])
  gains0.extend([1.,1.,1.,1.,1.])
  # mean_gain = readfile(indir+'/sigma_0.txt',[float]*8,cols=range(1,9))
  # sigma_gain = readfile(indir+'/sigma_0.txt',[float]*8,cols=range(9,17))

  f = TFile(indir+'/histos'+name+'.root')
  fout = TFile(os.environ['SCRATCH']+'/tmp.root','recreate')
  hr = [[f.Get(x).Get('hr_'+x+'_'+str(j)) for j in range(9)] for x in mlist]
  he = [[f.Get(x).Get('he_'+x+'_'+str(j)) for j in [1,2,5]] for x in mlist]
  heg = [[f.Get(x).Get('heg_'+x+'_'+str(j)) for j in [1,2,5]] for x in mlist]
  ht = [[f.Get(x).Get('htdiff_'+x+'_'+str(j)) for j in [1,2,5]] for x in mlist]

  # computation of gains and cie
  lerror,lindex,lalpha,lreso,lmean,lsigma,ldmean,ldsigma = [],[],[],[],[],[],[],[]
  for i,x in enumerate(mlist):
    print i, mlist[i]
    lerr,lj,le,lm,ls,ldm,lds = [],[],[],[],[],[],[]
    for j in range(9):
      if x=='G900' or x=='W835' or (x[0]=='W' and 476<int(x[1:])<775 and j>4):
        continue
      hr[i][j].GetXaxis().SetRangeUser(0.8,2.)
      if x=='G101' or x=='G131': hr[i][j].GetXaxis().SetRangeUser(0.5,2.)
      peak = hr[i][j].GetBinCenter(hr[i][j].GetMaximumBin())
      std = hr[i][j].GetStdDev()
      integ = hr[i][j].Integral()
      hr[i][j].Fit('gaus','qwwr','',peak-3*std,peak+3*std)
      ffit = hr[i][j].GetFunction('gaus')
      if integ<2 or ffit==None:
        print 'no stat for {0} ({1}<E<{2})'.format(x,200+100*j,300+100*j)
        lerr.append(j)
        continue
      m = ffit.GetParameter(1)
      if m<=0.:
        print 'negative fit for {0} ({1}<E<{2})'.format(x,200+100*j,300+100*j)
        lerr.append(j)
        continue
      s = ffit.GetParameter(2)
      dm = s/(integ-1)**0.5
      ds = s/(2*(integ-1))**0.5
      lj.append(j); le.append(250+100*j); lm.append(m); ls.append(s); ldm.append(dm); lds.append(ds)

    if lj==[]:
      print 'no stat at all for module '+x
      lalpha.append(0);  lreso.append(0)
    else:
      gm = tgraph(le,lm,[0. for x in lj],ldm)
      ffit = TF1('ffit_gain_'+mlist[i],'[0]+[1]*x*0.001',250.,1000.)
      gm.Fit(ffit,'q')
      if -1<ffit.GetParameter(1)<1: lalpha.append(ffit.GetParameter(1))
      else: 
        print 'non linearity not found for', mlist[i], ffit.GetParameter(1)
        lalpha.append(0)
      gs = tgraph(le,ls,[0. for x in lj],lds)
      ffit = TF1('ffit_reso_'+mlist[i],'[0]/sqrt(x*0.001)',250.,1000.)
      gs.Fit(ffit,'q')
      if 0.015<ffit.GetParameter(0)<0.2: lreso.append(ffit.GetParameter(0))
      else: 
        print 'resolution not found for', mlist[i], ffit.GetParameter(0)
        lreso.append(0)
    lindex.append(lj)
    lerror.append([j if j in lerr else -1 for j in range(9)])
    lmean.append([lm[lj.index(j)] if j in lj else 0 for j in range(9)])
    lsigma.append([ls[lj.index(j)] if (j in lj and ls[lj.index(j)]<0.3) else 0 for j in range(9)])
    ldmean.append([ldm[lj.index(j)] if j in lj else 0 for j in range(9)])
    ldsigma.append([lds[lj.index(j)] if j in lj else 0 for j in range(9)])
    
  writelist(indir+'/mean'+name+'.txt',[[x]+y for x,y in zip(mlist,lmean)])
  writelist(indir+'/sigma'+name+'.txt',[[x]+y for x,y in zip(mlist,lsigma)])
  writelist(indir+'/dmean'+name+'.txt',[[x]+y for x,y in zip(mlist,ldmean)])
  writelist(indir+'/dsigma'+name+'.txt',[[x]+y for x,y in zip(mlist,ldsigma)])
  writelist(indir+'/error'+name+'.txt',[[x]+y for x,y in zip(mlist,lerror)])
  writelist(indir+'/results'+name+'.txt',[['module','gain_550_600','alpha','reso']]+[[a,gains0[i]/b[3],c,d] if b[3]!=0 else [a,0.25,c,d] for i,(a,b,c,d) in enumerate(zip(mlist,lmean,lalpha,lreso))])

  # creating files with trigger stats for efficiency calculation
  ln1 = [[int(heg[i][0].Integral(201+50*j,201+50*(j+1))) for j in range(17)] for i in range(1733)]
  ln2 = [[int(heg[i][1].Integral(201+50*j,201+50*(j+1))) for j in range(17)] for i in range(1733)]
  ln5 = [[int(heg[i][2].Integral(201+50*j,201+50*(j+1))) for j in range(17)] for i in range(1733)]
  writelist(indir+'/stats'+name+'.txt',[[w]+x+y+z for w,x,y,z in zip(mlist,ln1,ln2,ln5)])
  ln1 = [[int(heg[i][0].Integral(201+50*j,201+50*(j+1))) if j/2 in lindex[i] else 0 for j in range(17)] for i in range(1733)]
  ln2 = [[int(heg[i][1].Integral(201+50*j,201+50*(j+1))) if j/2 in lindex[i] else 0 for j in range(17)] for i in range(1733)]
  ln5 = [[int(heg[i][2].Integral(201+50*j,201+50*(j+1))) if j/2 in lindex[i] else 0 for j in range(17)] for i in range(1733)]
  writelist(indir+'/stats_cleaned'+name+'.txt',[[w]+x+y+z for w,x,y,z in zip(mlist,ln1,ln2,ln5)])

  # tdiff calculation
  ldiff = [[f.Get(x).Get('htdiff_'+x+'_'+str(i)).GetMean() for i in [1,2]] if x[0]=='G' else [0.,f.Get(x).Get('htdiff_'+x+'_2').GetMean()] for x in module_names]
  ldiff = [[a]+[int(round(y,1)) for y in x] for a,x in zip(module_names,ldiff)]
  writelist(indir+'/tdiff'+name+'.txt',ldiff)

  # creating calibration files
  if write_gain!='':
    lms = [get_lms_gains(889,1),get_lms_gains(889,2),get_lms_gains(889,3)]
    lms = [[lms[i][j] for i in range(3)] for j in range(1728)]
    if metjod=='fsnake':
      lf = [[x,'%.4f'%(gains0[i]/y[3])]+z+[0.0,'%.4f'%b] if y[3]!=0 else [x,'%.4f'%(0.25)]+z+[0.0,'%.4f'%b] for i,(x,y,b,z) in enumerate(zip(mlist,lmean,lalpha,lms))]
    else:
      lf = [[x,'%.4f'%(gains0[i]/y[3]),'%.4f'%b]+z if y[3]!=0 else [x,'%.4f'%(0.25),'%.4f'%b]+z for i,(x,y,b,z) in enumerate(zip(mlist,lmean,lalpha,lms))]
    writelist(gaindir_in.replace('pass'+ipass,'pass'+str(int(ipass)+1)+'/calib_889-979.txt',lf,ljust=10)
