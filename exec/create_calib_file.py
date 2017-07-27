#!/apps/python/PRO/bin/python

from readdst import *

run = ''
output = work_folder+'/prad'
indir = replay_folder
batch = False
show = 0

# input parameters
args = sys.argv
if len(args)<2: sys.exit('no arguments')
for i,a in enumerate(args):
  if a in ['-r','-run']: run = args[i+1]
  elif a in ['-o','-out']: output = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-i','-in']: indir = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-b','-batch']: batch = True
  elif a in ['-s','-show']: show = int(args[i+1])
  
if run=='': sys.exit('no run mentionned')

# batch mode
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+'_'+run,command=os.getcwd()+'/create_calib_file.py',options='-r '+run+' -i . -o '+os.path.basename(output)+'_'+run,input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run.zfill(6)+'.dst',output_template=os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='25 GB',memory='1024 MB')
  print
  sys.exit()

# opening file and skip LMS part
t = chao(run=int(run),folder=indir)
t.goto(20000)
t.lload = ['','adc']

# output file
f = open(output+'.dst','wb')
f.write(pack('I',int(run)))

# loop on events
n = 0
for ev in progress(t,modulo=10000,show=show):
  if t.isepics: continue
  if ev.trigger not in [1,2]: continue
  n+=1
  f.write(pack('I',ev.iev))
  f.write(pack('B',ev.trigger))
  f.write(pack('f',t.epics['MBSY2C_energy']))
  f.write(pack('H',ev.n_adc))
  for i in range(ev.n_adc):
    f.write(pack('H',ev.adc_id[i]))
    f.write(pack('H',ev.adc_val[i]))

f.write(pack('I',n))

f.close()
t.close()
sys.exit()
