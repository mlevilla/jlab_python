import os,sys,time,subprocess
from math import log,exp,cos,sin,tan,atan2,asin,acos,pi,degrees,radians
from random import randint,gauss
degrad = pi/180.

#----------------------------------------------------------------------------------------
# progress bar 
#----------------------------------------------------------------------------------------

def write_and_flush(st):
  sys.stdout.write(st)
  sys.stdout.flush()
  sys.stdout.write('\r')

def parse_time(t): return '%.2i'%((int(t)/60)/60)+':'+'%.2i'%((int(t)/60)%60)+':'+'%.2i'%(int(t)%60)

def progress_bar(show,count,precision,modulo,i=None,n=None,width=None,m=None,t0=None):
  if show==0: return ''
  if n==-1: # no info on length
    st=str(i)
  elif show==1: # just percentage
    st = ('%0.'+str(precision)+'f')%(count*100/float(modulo))+'%'
  elif show==2: # progress bar + percentage
    st = m[0]+m[1]*(((i+1)*width)/n)+' '*(width-((i+1)*width)/n)+m[2]+' '+('%0.'+str(precision)+'f')%(count*100/float(modulo))+'%'
  elif show==3: # progress bar + percentage + time
    ti = time.time()
    st = m[0]+m[1]*(((i+1)*width)/n)+' '*(width-((i+1)*width)/n)+m[2]+' '+('%0.'+str(precision)+'f')%(count*100/float(modulo))+'%  - time: '+parse_time(ti-t0)+' - eta: '+parse_time((ti-t0)*(n-i-1)/(i+1))
  write_and_flush(st)
  return st

def progress(l,modulo=100,precision=0,show=2,width=30,m='<=>',n=-1,limit=1.,insert=None):
  # find length according to type of object
  if isinstance(l,int) or isinstance(l,long): l=range(l)
  if n==-1:
    try: n = len(l)
    except: n = -1
  if n!=-1 and n<width: width=n
  # options
  if modulo!=100: precision = int(log(modulo)/log(10))-2
  if precision!=0: modulo = int(10**(2+precision))
  if n!=-1: modulo = min(modulo,n)
  if limit<=1. and n!=-1: limit *= n
  if limit==1. and n==-1: limit = sys.maxint
  # loop
  t0 = time.time()
  count = 0
  for i,x in enumerate(l):
    if i>=limit: break
    yield x
    if show==0: continue
    if (n==-1 and i%modulo==0) or (n!=-1 and i>=count*n/float(modulo)):
      st = progress_bar(show,count,precision,modulo,i,n,width,m,t0)
      count+=1
    # insert
    if insert!=None:
      ins = str(insert(x))
      n_ins = len(ins)+2
      st = st[:width/2-n_ins/2]+ins.center(n_ins)+st[width/2-n_ins/2+n_ins:]
      write_and_flush(st)
  # end of loop
  if n==-1: pass
  progress_bar(show,count,precision,modulo,i,n,width,m,t0)
  sys.stdout.write('\n')


#----------------------------------------------------------------------------------------
# input/output utils 
#----------------------------------------------------------------------------------------
 
def readfile(filename,types=[],out=list,start=0,end=0,cols=[],rows=[],tr=False,space=' ',header=[],simple=True,raw=False):

  # test access to file
  if not os.access(filename,os.R_OK):
    sys.stderr.write(filename+' file not readable\n')
    return []

  # read all the lines
  f = open(filename,'r')
  l = f.readlines()
  f.close()

  # taking care of rows
  if not l or end>len(l) or (rows and max(rows)>len(l)-1):
    sys.stderr.write('not enough rows in the file '+filename+'\n')
    return l
  end = (end-1)%len(l)+1
  if not rows: rows = range(start,end)

  # creating list of lists
  l2 = []
  for i in rows:
    x = l[i].split(space)
    while '' in x: del x[x.index('')]
    if x[-1][-1]=='\n': x[-1] = x[-1][:-1]
    l2.append(x)
  l = l2

  if raw: return l
  
  # handling columns
  if not cols and not types:
    cols = range(len(l[0]))
    types = [str]*len(l[0])
  elif not cols: cols = range(len(types))
  elif not types: types = [str]*len(cols)
  elif len(types)>len(cols): types=types[:len(cols)]
  elif len(types)<len(cols): types.extend([str]*(len(cols)-len(types)))

  # formatting header
  if header: 
    keys = [x(y) for x,y in zip(header,[l[0][j] for j in cols if j<len(l[0])])]
    del l[0]
  # formatting for list of list
  for i in range(len(l)):
    l[i] = [x(y) for x,y in zip(types,[l[i][j] for j in cols if j<len(l[i])])]
    if simple and len(l[i])==1: l[i] = l[i][0]
  # transposing
  if tr: l = [[x[j] for x in l] for j in range(len(l[0]))]
  # changing to list of dicts
  elif header:
    for i in range(len(l)): l[i] = {x:y for x,y in zip(keys,l[i])}

  if simple and len(l)==1: l = l[0]
  
  # formatting for dict ouput 
  if out==dict:
    if header: l = {keys[i]:x[0] if (simple and len(x)==1) else x for i,x in enumerate(l)}
    else: l = {x[0]:x[1] if (simple and len(x)==2) else x[1:] for x in l}
  return l


def writelist(path,l,tr=False,space=' ',ljust=None,ro=-1):
  f = open(path,'w')
  if not isinstance(l[0],list): l = [[x] for x in l]
  if ro!=-1: l = [[round(y,ro) if isinstance(y,float) else y for y in x] for x in l]
  if tr: l = [[x[j] for x in l] for j in range(len(l[0]))]
  for x in l:
    if ljust is None: f.write(space.join([str(y) for y in x])+'\n')
    elif isinstance(ljust,int): f.write(space.join([str(y).ljust(ljust) for y in x])+'\n')
    elif isintance(ljust,list) and len(ljust)==len(l[0]): f.write(space.join([str(y).ljust(ljust[i]) for i,y in enumerate(x)])+'\n')
    else: sys.stdout.write('ljsut with wrong length\n')
  f.close()


#----------------------------------------------------------------------------------------
# argument parser 
#----------------------------------------------------------------------------------------

# largs0 = [[('sp',0,'','show_progress'),
#           ('b',False,'','batch_mode'),
#           ('wt','120','','wall_time')],
#           [('iv',False,'','interactive_mode')]]
# show_progress, batch, interactive, wtime = [0], [False], [False], ['120']

# def catch_arg(mark,default,opt='',args=[]):
#   if args==[]: args = sys.argv
#   if not isinstance(default,list): default = [default]
#   if opt=='several': r = [default]
#   else: r = default
#   if '-'+mark in args:
#     i = args.index('-'+mark)
#     if opt=='':
#       if isinstance(default[0],bool): r = [not default[0]]
#       else: r = [type(default[j])(args[i+j+1]) for j in range(len(default)) if i+j+1<len(args)]
#     elif opt=='several':
#       k = i+1
#       while k<len(args) and (args[k]=='' or args[k][0]!='-'): k+=1
#       r = [[type(default[j])(args[i+1+len(default)*n+j]) for n in range((k-i-1)/len(default))] for j in range(len(default))]
#     elif opt=='end':
#       r = [' '.join(args[i+1:])]
#   if len(default)==1: r = r[0]
#   return r

# def catch_args(largs):
#   if '.in' in sys.argv[1]: args = readfile(sys.argv[1])
#   else: args = sys.argv
#   linput = []
#   for x in largs:
#     if x[0] in [y[0] for y in largs0[0]+largs0[1]]:
#       print 'warning: standard argument -'+x[0]+' rewritten'
#     linput.append(catch_arg(*x[:-1],args=args))
#   global show_progress, batch, interactive, wtime          
#   for x,y in zip(largs0[0]+largs0[1],[show_progress,batch,wtime,interactive]):
#     if x[0] in [a[0] for a in largs]: continue
#     largs.append(x)
#     y[0] = catch_arg(*x[:-1],args=args)
#   return linput

# def print_help(largs,empty=1,lmax=15):
#   if '-h' in sys.argv or '--help' in sys.argv or (empty and len(sys.argv)==1):
#     print 'command:',sys.argv[0]
#     print 'inputs:'
#     max_length0 = max(len(x[0]) for x in largs)
#     max_length1 = max(min(len(str(x[1])),lmax+3)+len(x[3]) for x in largs)
#     for x,y,z,c in largs:
#       mark = '\'' if isinstance(y,str) else ''
#       if len(str(y))>lmax: y = '...'+str(y)[-15:]
#       sys.stdout.write('  '+('-'+x).ljust(max_length0+1)+'\t'+(' < '+c+' = '+mark+str(y)+mark+' > ').ljust(max_length1+9)+('\t mode: '+{'several':'several','end':'end','':''}[z])*int(z!='')+'\n')
#     sys.exit()

# def del_args(largs,ldel,join=0):
#   args = sys.argv[1:]
#   for x in largs:
#     if x[0] not in ldel: continue
#     if '-'+x[0] not in args: continue
#     i = args.index('-'+x[0])
#     del args[i]
#     if x[2]=='several':
#       while args[i][0]!='-': del args[i]
#     elif x[2]=='end': args = args[:i]
#     elif isinstance(x[1],list):
#       for x in x[1]: del args[i]
#     elif not isinstance(x[1],bool): del args[i]
#   if join: args = ' '.join(args)
#   return args

# def add_args_list(args,mark,l,write=''):
#   if not isinstance(mark,list): mark = [mark]
#   idx,arg0 = [],[]
#   for x in mark:
#     if '-'+x in args:
#       i = args.index('-'+x)+1
#       while i<len(args) and args[i][0]!='-':
#         idx.append(i)
#         arg0.append(args[idx[-1]])  
#         i+=1
#     else:
#       args.extend(['-'+x,''])
#       idx.append(len(args)-1)
#       arg0.append(args[idx[-1]])  
#   largs = []
#   for x in l:
#     for i,a in zip(idx,arg0): args[i] = a+str(x)
#     if write!='':
#       f = open(write+str(x)+'.in','w')
#       f.write(' '.join(args))
#       f.close()
#     largs.append(args)
#   if write=='': return largs
#   else: return ' '.join([write+str(x)+'.in' for x in l])


#----------------------------------------------------------------------------------------
# other utils 
#----------------------------------------------------------------------------------------
 
def cut(b,name,lv,d,binning):
  if not b: return True
  if not isinstance(lv,list): lv = [lv]
  for x in lv:
    if v<binning[0]: d[name][0]+=1
    elif v>=binning[-1]: d[name][-1]+=1
    else:
     for i in range(len(binning)-1):
       if binning[i]<=v<binning[i+1]:
         d[name][i+1]+=1
         break
  return False

def unnest(l):
  l2 = []
  for x in l:
    if isinstance(x,list): l2.extend(unnest(x))
    else: l2.append(x)
  return l2  


#----------------------------------------------------------------------------------------
# jlab farm submission 
#----------------------------------------------------------------------------------------

# def jsub(**args):
#   necessary = ['project','track','command','jobname']
#   for x in necessary:
#     if x not in args:
#       sys.stderr.write('necessary argument "'+x+'" not found\n')
#       return
#   f = open(os.environ['JOBS']+'/'+args['jobname']+'.txt','w')
#   for x,y in args.items():
#     if isinstance(y,str): f.write(x.upper()+': '+y+'\n')
#     elif isinstance(y,list): f.write(x.upper()+': '+' '.join(y)+'\n')
#   f.close()
#   os.system('jsub '+os.environ['JOBS']+'/'+args['jobname']+'.txt')
#   sys.stdout.write('job '+args['jobname']+' submitted!\n')

# def jsub_fast(largs,name='job',rem=[],add=[[],[]],outdir='/work/hallb/prad/mlevilla',extend=['-'+x[0] for x in largs0[1]],output_suff='_*.*',jobname=None,**oargs):
#   rem.extend([x[0] for x in largs0[0]])
#   args = del_args(largs,rem)
#   args.extend(extend)
#   if jobname is None: jobname = name
#   input_files = add_args_list(args,add[0],add[1],write=os.environ['JOBS']+'/'+jobname)
#   if 'memory' not in oargs: oargs['memory']='2048 MB'
#   if 'os' not in oargs: oargs['os']='centos7'
#   if 'time' not in oargs: oargs['time']='300'
#   jsub(project='prad',track='analysis',command=os.path.abspath(sys.argv[0])+' args.in',jobname=jobname,input_files=input_files,input_data='args.in',output_data=name+output_suff,output_template=outdir+'/@OUTPUT_DATA@',**oargs)
