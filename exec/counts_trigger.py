#!/apps/python/PRO/bin/python

from params import *

run = ''
output =  work_folder+'/counts_trigger'
indir = replay_folder_me
batch = False
show = 0

args = sys.argv
if len(args)<3: sys.exit('no arguments')
for i,a in enumerate(args):
  if a in ['-r','-run']: run = args[i+1]
  elif a in ['-o','-out']: output = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-i','-in']: indir = args[i+1] if (args[i+1][0]=='/') else os.getcwd()+'/'+args[i+1]
  elif a in ['-b','-batch']: batch = True
  elif a in ['-s','-show']: show = int(args[i+1])

if run=='': sys.exit('no run mentionned')

# batch mode
if batch:
  jsub(project='prad',track='analysis',jobname=os.path.basename(output)+'_'+run,command=os.getcwd()+'/counts_trigger.py',options='-o '+os.path.basename(output)+' -i . -r '+run,input_files=indir+'/prad_'+run.zfill(6)+'.dst',output_data=os.path.basename(output)+'_'+run+'*',output_template==os.path.dirname(output)+'/@OUTPUT_DATA@',disk_space='25 GB',memory='1024 MB')
  print
  sys.exit()

#initialisation
from readdst import *
from tagger import *
from hycal import *

t = chao(run=int(run),folder=indir,step=100000,lload=['','tdc'])
triggers0 = [0 for i in range(8)]
triggers1 = [0 for i in range(8)]
for ev in progress(t,modulo=10000,show=show):
  if t.isepics: continue
  if event_is_bad(t.e.iev,t.run): continue
  triggers0[ev.trigger]+=1
  #etc = [x for x in get_etchannels(ev=ev) if abs(x.t)<25]
  #if etc==[]: continue
  #triggers1[ev.trigger]+=1
  
  
f = open(output+'_'+run+'.txt','w')
f.write(run+'\n')
f.write(str(t.i)+'\n')
f.write(' '.join([str(t.ancre[i]) for i in range(len(t.ancre))])+'\n')
f.write(' '.join([str(x) for x in triggers0])+'\n')
f.write(' '.join([str(x) for x in triggers1])+'\n')
f.close()

sys.exit(0)
