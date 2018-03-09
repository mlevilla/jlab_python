import os
from arg_parser import *

def jsub(**args):
  necessary = ['project','track','command','jobname']
  for x in necessary:
    if x not in args:
      sys.stderr.write('necessary argument "'+x+'" not found\n')
      return
  args['jobname'] = args['jobname'][:50]
  print os.environ['JOBS']+'/'+args['jobname']+'.txt'
  f = open(os.environ['JOBS']+'/'+args['jobname']+'.txt','w')
  for x,y in args.items():
    if isinstance(y,str): f.write(x.upper()+': '+y+'\n')
    elif isinstance(y,list): f.write(x.upper()+': '+' '.join(y)+'\n')
  f.close()
  os.system('jsub '+os.environ['JOBS']+'/'+args['jobname']+'.txt')
  sys.stdout.write('job '+args['jobname']+' submitted!\n')

def jsub_fast(aparser,name='job',rem=[],add=[[],[]],outdir='/work/hallb/prad/mlevilla',extend=[],output_suff='_*.*',jobname=None,**oargs):
  rem.extend([x[0] for x in largs0[0]])
  aparser.delete(rem)
  aparser.argv.extend(extend+['-'+x[0] for x in largs0[1]])
  if jobname is None: jobname = name
  input_files = aparser.add_list(add[0],add[1],write=os.environ['JOBS']+'/'+jobname)
  if 'memory' not in oargs: oargs['memory']='2048 MB'
  if 'os' not in oargs: oargs['os']='centos7'
  if 'time' not in oargs: oargs['time']='300'
  jsub(project='prad',track='analysis',command=os.path.abspath(sys.argv[0])+' args.in',jobname=jobname,input_files=input_files,input_data='args.in',output_data=name+output_suff,output_template=outdir+'/@OUTPUT_DATA@',**oargs)
