#!/apps/python/PRO/bin/python

from params import *
from readdst import *
from tagger import *

run = ''
output = work_folder+'/prad'
indir = replay_folder_me
batch = False
show = 0
calib = 1

# input parameters
args = sys.argv
if len(args)<3: sys.exit('no arguments')
for i,a in enumerate(args):
  if a in ['-r','-run']: run = args[i+1]
  elif a in ['-o','-out']: output = os.path.abspath(args[i+1])
  elif a in ['-i','-in']: indir = os.path.abspath(args[i+1])
  elif a in ['-b','-batch']: batch = True
  elif a in ['-s','-show']: show = int(args[i+1])
  elif a in ['-c','-calib']: calib = int(args[i+1])
  
if run=='': sys.exit('no run mentionned')

# batch mode
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+run,command=os.getcwd()+'/toIlya.py',options='-r '+run+' -i . -o '+os.path.basename(output)+' -c '+str(calib),input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run+'_ilya.dst',output_template=os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='25 GB',memory='1024 MB')
  print
  sys.exit()
  
# opening file and skip LMS part
t = chao(run=int(run),folder=indir)
t.goto(20000)
t.lload = ['','adc','tdc'] if calib else ['','adc']
# output file
f = open(output+'_'+run+'_ilya.dst','wb')

# loop on events
for ev in progress(t,modulo=10000,show=show):
  if t.isepics: continue
  if ev.trigger not in [1,2,5]: continue

  # get ET matching photons (or not)
  if calib:
    etc = TDCHits([x for x in get_etchannels(ev=ev) if abs(x.t)<25])
    if etc==[]: 
      hastagger = 0
      eg = t.epics['MBSY2C_energy']
    else:
      hastagger = 1
      etc = etc.merge()
      eg = t.epics['MBSY2C_energy']*etc[0].E
  else: 
    hastagger = 0
    eg = t.epics['MBSY2C_energy']

  # get hycal hits
  adc,old_id = [],[]
  for i in range(ev.n_adc):
    if ev.adc_id[i]>1727: continue
    name = module_names[ev.adc_id[i]]
    old_id.append(int(name[1:])+1000*int(name[0]=='W'))
    adc.append(ev.adc_val[i])
    
  f.write(pack('I',t.run))
  f.write(pack('I',ev.iev))
  f.write(pack('I',ev.time/25000))
  f.write(pack('f',eg))
  f.write(pack('I',ev.trigger))
  f.write(pack('I',hastagger))
  xpos = t.epics['hallb_ptrans_x_encoder']+652.05
  ypos = -(t.epics['hallb_ptrans_y1_encoder']+t.epics['hallb_ptrans_y2_encoder'])/2.-3761.5
  f.write(pack('i',int(1000*xpos)))
  f.write(pack('i',int(1000*ypos)))
  f.write(pack('I',len(adc)))
  for i in range(len(adc)): f.write(pack('I',old_id[i]<<16 | adc[i]))
   
f.close()
t.close()
sys.exit(0)
