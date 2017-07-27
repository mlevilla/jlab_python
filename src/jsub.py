import os,sys

def jsub(**args):
  necessary = ['project','track','command','jobname']
  for x in necessary:
    if x not in args:
      sys.stderr.write('necessary argument "'+x+'" not found\n')
      return
  f = open(os.environ['JOBS']+'/'+args['jobname']+'.txt','w')
  for x,y in args.items():
    if isinstance(y,str): f.write(x.upper()+': '+y+'\n')
    elif isinstance(y,list): f.write(x.upper()+': '+' '.join(y)+'\n')
  f.close()
  os.system('jsub '+os.environ['JOBS']+'/'+args['jobname']+'.txt')
  sys.stdout.write('job '+args['jobname']+' submitted!\n')

def jsub_fast(largs,name='job',rem=[],add=[[],[]],outdir='/work/hallb/prad/mlevilla',extend=[],output_suff='_*.*',**oargs):
  args = del_args(largs,rem)
  args.extend(extend)
  input_files = add_args_list(args,add[0],add[1],write=os.environ['JOBS']+'/'+name)
  if 'memory' not in oargs: oargs['memory']='2048 MB'
  if 'os' not in oargs: oargs['os']='centos7'
  jsub(project='prad',track='analysis',command=os.path.abspath(sys.argv[0])+' args.in',jobname=name,input_files=input_files,input_data='args.in',output_data=name+output_suff,output_template=outdir+'/@OUTPUT_DATA@',**oargs)
