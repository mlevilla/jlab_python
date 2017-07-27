#!/apps/anaconda/anaconda-2.0.1/bin/python

from params import *
from misc_root import *

largs = [('r','','','run_number'),
         ('o','.','','output_folder'),
         ('i','.','','input_folder'),
         ('b',False,'','batch_mode'),
         ('f',False,'','final_calculation'),
         ('s',0,'','show_progress_bar')]

print_help(largs)
[run,outdir,indir,batch,final,show] = catch_args(largs)

load_names()
for i in range(1,3): module_names[-i] = ['LMS'+str(i)]

# group all runs
if final:
  runs = sum_dict_lists(runs_list)
  f = TFile(output+'.root','recreate')
  f.cd()
  gped = [TGraphErrors(len(runs)) for i in range(1731)]
  gped_sigma = [TGraphErrors(len(runs)) for i in range(1731)]
  gled = [TGraphErrors(len(runs)) for i in range(1731)]
  galpha = [TGraphErrors(len(runs)) for i in range(3)]
  gref = [TGraphErrors(len(runs)) for i in range(3)]
  ggains = [[TGraphErrors(len(runs)) for j in range(3)] for i in range(1728)]
  for i in range(1731):
    gped[i].SetName('gped_'+module_names[i])
    gped_sigma[i].SetName('gped_sigma_'+module_names[i])
    gled[i].SetName('gled_'+module_names[i])
    gped[i].SetTitle('PED '+module_names[i]+';#run;ped')
    gped_sigma[i].SetTitle('PED SIGMA'+module_names[i]+';#run;#sigma_{ped}')
    gled[i].SetTitle('LED '+module_names[i]+';#run;led')
  for i in range(1728):
    for j in range(3):
      ggains[i][j].SetName('ggains_'+module_names[i]+'_LMS'+str(j+1))
      ggains[i][j].SetTitle('Gain '+module_names[i]+' wrt ref '+str(j+1)+';#run;gain')
  for i in range(3):
    galpha[i].SetName('galpha_LMS'+str(i+1))
    gref[i].SetName('gref_LMS'+str(i+1))
    galpha[i].SetTitle('alpha '+str(i+1)+';#run;alpha')
    gref[i].SetTitle('ref '+str(i+1)+';#run;alpha/lms')

  for irun,run in enumerate(runs):
    print run
    frun = TFile(output+'_'+str(run)+'.root')
    for i in range(1728):
      u = frun.hycal.Get(module_names[i])
      fped = u.Get('hped_'+module_names[i]+'_'+str(run)).GetFunction('gaus')
      gped[i].SetPoint(irun,run,fped.GetParameter(1))
      gped[i].SetPointError(irun,0.,fped.GetParError(1))
      gped_sigma[i].SetPoint(irun,run,fped.GetParameter(2))
      gped_sigma[i].SetPointError(irun,0.,fped.GetParError(2))
      fled = u.Get('hled_'+module_names[i]+'_'+str(run)).GetFunction('gaus')
      gled[i].SetPoint(irun,run,fled.GetParameter(1))
      gled[i].SetPointError(irun,0.,fled.GetParError(1))
    for i in range(3):
      u = frun.hycal.Get(module_names[i+1728])
      fled = u.Get('hled_'+module_names[i+1728]+'_'+str(run)).GetFunction('gaus')
      halpha = u.Get('halpha_'+module_names[i+1728]+'_'+str(run))
      halpha.Fit('gaus','rwwq0','',0,800)
      if run in [982,983,987,988]: fped = halpha.GetFunction('gaus')
      else: fped = u.Get('hped_'+module_names[i+1728]+'_'+str(run)).GetFunction('gaus')
      ped = fped.GetParameter(1)
      dped = fped.GetParameter(2)/100.
      led = fled.GetParameter(1)
      dled = fled.GetParameter(2)/100.
      falpha = TF1('falpha','gaus',1000.,8000.)
      halpha.Fit(falpha,'rwwq0')
      alpha = falpha.GetParameter(1) if falpha!=None else 0.
      nalpha = sum([halpha[k] for k in range(halpha.GetXaxis().FindBin(ped+200),halpha.GetXaxis().FindBin(16000))])
      if nalpha<=1: alpha,ref,dalpha,dref = 0.,0.,0.,0.
      else:
        dalpha = falpha.GetParameter(2)/(nalpha-1)**0.5 if falpha!=None else 0.
        ref = (alpha-ped)/(led-ped)
        dref = ref*((dalpha**2+dped**2)/(alpha-ped)**2+(dled**2+dped**2)/(led-ped)**2)**0.5
      galpha[i].SetPoint(irun,run,alpha)
      galpha[i].SetPointError(irun,0.,dalpha)
      gref[i].SetPoint(irun,run,ref)
      gref[i].SetPointError(irun,0.,dref)
      for j in range(1728):
        u = frun.hycal.Get(module_names[j])
        fped = u.Get('hped_'+module_names[j]+'_'+str(run)).GetFunction('gaus')
        fled = u.Get('hled_'+module_names[j]+'_'+str(run)).GetFunction('gaus')
        ped = fped.GetParameter(1)
        dped = fped.GetParError(1)
        led = fled.GetParameter(1)
        dled = fled.GetParError(1)
        if nalpha<=1: gain,dgain = 0.,0.
        else:
          gain = ref*(led-ped)
          dgain = gain*(dref**2/ref**2+(dled**2+dped**2)/(led-ped)**2)**0.5
        ggains[j][i].SetPoint(irun,run,gain)
        ggains[j][i].SetPointError(irun,0.,dgain)

  for i in range(1728):
    folder = f.mkdir(module_names[i])
    folder.cd()
    gped[i].Write()
    gped_sigma[i].Write()
    gled[i].Write()
    for j in range(3):
      ggains[i][j].Write()
  for i in range(3):
    folder = f.mkdir('LMS'+str(i+1))
    folder.cd()
    gped[1728+i].Write()
    gped_sigma[1728+i].Write()
    gled[1728+i].Write()
    galpha[i].Write()
    gref[i].Write()

  f.Close()
  sys.exit('runs grouped')


if run=='': sys.exit('no run mentionned')

# batch mode
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+'_'+run,command=os.getcwd()+'/stability.py',options='-o '+os.path.basename(output)+' -i . -r '+run,input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run+'*',output_template=os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='25 GB',memory='1024 MB')
  print
  sys.exit()

# initialisation
from readdst import *

f = TFile(output+'_'+run+'.root','recreate')
t = chao(run=int(run),folder=indir,lload=['','adc','dsc'])
f.cd()
hped = [TH1F('hped_'+module_names[i]+'_'+run,';ped;#counts',1024,0,1024) for i in range(1731)]
hled = [TH1F('hled_'+module_names[i]+'_'+run,';LMS;#counts',4096,0,16384) for i in range(1731)]
halpha = [TH1F('halpha_'+module_names[i]+'_'+run,';alpha;#counts',2048,0,8192) for i in range(1728,1731)]
l_epics,l_dsc = [],[]

# loop on events
for ev in progress(t,show=show,modulo=10000):
  if t.isepics: l_epics.append([t.epics['iev']]+[t.epics[x] for x in t.epics_name])
  elif t.e.typ==0xc1:
    time = t.e.ungated[7]/500000.
    live_time = 1-t.e.gated[7]/t.e.ungated[7]
    beam_charge =  (t.e.ungated[6]-100.)/906.2
    beam_current = beam_charge/time
    l_dsc.append([t.e.iev,live_time,beam_charge,beam_current])
  elif t.i<30000:
    if ev.trigger==3: # led
      for i in range(ev.n_adc):
        if ev.adc_id[i]<1731: hled[ev.adc_id[i]].Fill(ev.adc_val[i])
    elif ev.trigger==4: # alpha 
      for i in range(ev.n_adc):
        if ev.adc_id[i]<1728: hped[ev.adc_id[i]].Fill(ev.adc_val[i])
        elif 1728<=ev.adc_id[i]<1731: halpha[ev.adc_id[i]-1728].Fill(ev.adc_val[i])
    elif ev.trigger in [1,2,5,6]: #physics
      for i in range(ev.n_adc):
        if 1728<=ev.adc_id[i]<1731: hped[ev.adc_id[i]].Fill(ev.adc_val[i])

# fit
ped_mean,ped_sigma,led_mean,led_sigma,alpha_mean,alpha_sigma = [],[],[],[],[],[]
for i in range(1731):
  cmin,cmax = -10,10
  if module_names[i]=='W6': cmin,cmax=-60,60
  elif module_names[i] in ['W605','W604','W575','W269']: cmin,cmax=-30,30
  elif module_names[i] in ['G632','G238','G171','G55','G41']: cmin,cmax=-20,10

  peak = hped[i].GetBinCenter(hped[i].GetMaximumBin())
  hped[i].Fit('gaus','rqww','',peak+cmin,peak+cmax)
  fped = hped[i].GetFunction('gaus')
  if fped!=None:
    ped_mean.append(fped.GetParameter(0))
    ped_sigma.append(fped.GetParError(0))
  else: 
    print 'no stat for module',module_names[i]
    ped_mean.append(1)
    ped_sigma.append(0)
  hled[i].Fit('gaus','qww')
  fled = hled[i].GetFunction('gaus')
  if fled!=None:
    led_mean.append(fled.GetParameter(1))
    led_sigma.append(fled.GetParError(0))
  else: 
    print 'no stat for module',module_names[i]
    led_mean.append(1)
    led_sigma.append(0)
for i in range(3):
  halpha[i].Fit('gaus','rqww','',1000.,8192.)
  falpha = halpha[i].GetFunction('gaus')
  if falpha!=None:
    alpha_mean.append(falpha.GetParameter(1))
    alpha_sigma.append(falpha.GetParError(1))
  else: 
    print 'no stat for module',module_names[i]
    alpha_mean.append(1)
    alpha_sigma.append(0)

# writing histograms
folder = f.mkdir('hycal')
folder.cd()
for i in range(1728):
  folder2 = folder.mkdir(module_names[i])
  folder2.cd()
  hped[i].Write()
  hled[i].Write()
for i in range(3):
  folder2 = folder.mkdir('LMS'+str(i+1))
  folder2.cd()
  halpha[i].Write()
  hped[1728+i].Write()
  hled[1728+i].Write()

# writing epics and dsc graphs
g_epics = [TGraph(len(l_epics),array('f',[x[0] for x in l_epics]),array('f',[x[i] for x in l_epics])) for i in range(1,len(l_epics[0]))]
print l_dsc
if l_dsc!=[]:
  g_dsc = [TGraph(len(l_dsc),array('f',[x[0] for x in l_dsc]),array('f',[x[i] for x in l_dsc])) for i in range(1,len(l_dsc[0]))]
else: g_dsc = []
for x,y in zip(g_epics,t.epics_name): 
  x.SetName('g'+y)
  x.SetTitle(';#event;'+y)
for x,y in zip(g_dsc,['live_time','beam_charge','beam_current']): 
  x.SetName('g'+y)
  x.SetTitle(';#event;'+y)
folder = f.mkdir('epics')
folder.cd()
for x in g_epics: x.Write()
folder = f.mkdir('dsc')
folder.cd()
for x in g_dsc: x.Write()
f.Close()

t.close()

# calculation of gains
ref_mean = [(alpha_mean[i]-ped_mean[1728+i])/(led_mean[1728+i]-ped_mean[1728+i]) for i in range(3)]
ref_sigma = [ref_mean[i]*((alpha_sigma[i]**2+ped_sigma[1728+i]**2)/(alpha_mean[i]-ped_mean[1728+i])**2+(led_sigma[1728+i]**2+ped_sigma[1728+i]**2)/(led_mean[1728+i]-ped_mean[1728+i])**2)**0.5 for i in range(3)]
gains_mean = [[ref_mean[i]*(led_mean[j]-ped_mean[j]) for i in range(3)] for j in range(1728)]
gains_sigma = [[gains_mean[j][i]*(ref_sigma[i]**2/ref_mean[i]**2+(led_sigma[j]**2+ped_sigma[j]**2)/(led_mean[j]-ped_mean[j])**2)**0.5 for i in range(3)] for j in range(1728)]

# writing info for graphs
f_txt = open(output+'_'+run+'.txt','w')
f_txt.write('run '+run+'\n')
for i in range(1731):
  f_txt.write('ped '+str(i)+' '+str(ped_mean[i])+' '+str(ped_sigma[i])+'\n')
  f_txt.write('led '+str(i)+' '+str(led_mean[i])+' '+str(led_sigma[i])+'\n')
for i in range(3):
  f_txt.write('alpha '+str(i)+' '+str(alpha_mean[i])+' '+str(alpha_sigma[i])+'\n')
  f_txt.write('ref '+str(i)+' '+str(ref_mean[i])+' '+str(ref_sigma[i])+'\n')
for i in range(1728):
  for j in range(3):
    f_txt.write('gains '+str(i)+' '+str(j)+' '+str(gains_mean[i][j])+' '+str(gains_sigma[i][j])+'\n')
f_txt.close()
    
sys.exit(0)
  
