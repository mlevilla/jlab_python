import os,sys,time,subprocess,re,signal
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

  filename = catch_env_var(filename)
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


def writelist(path,l,tr=False,space=' ',ljust=None,ro=[]):
  f = open(catch_env_var(path),'w')
  if not isinstance(l[0],list): l = [[x] for x in l]
  if isinstance(ro,int): ro = [ro for x in l[0]]
  if ro!=[]: l = [[round(y,r) if isinstance(y,float) else y for y,r in zip(x,ro)] for x in l]
  if tr: l = [[x[j] for x in l] for j in range(len(l[0]))]
  for x in l:
    if ljust is None: f.write(space.join([str(y) for y in x])+'\n')
    elif isinstance(ljust,int): f.write(space.join([str(y).ljust(ljust) for y in x])+'\n')
    elif isinstance(ljust,list) and len(ljust)==len(l[0]): f.write(space.join([str(y).ljust(ljust[i]) for i,y in enumerate(x)])+'\n')
    else: sys.stdout.write('ljsut with wrong length\n')
  f.close()





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

def catch_env_var(x):
  i0 = 0
  while '$' in x[i0:]:
    i = x.index('$',i0)
    if i!=0 and x[i-1]=='\\': 
      i0 = i+1
      continue
    if i==len(x)-1: break
    if x[i+1]=='{':
      j1 = i+2
      j2 = x.index('}',i+2)
      j3 = j2+1
    else:
      j1 = i+1
      j2 = x.index('/',i+1) if '/' in x[i+1:] else len(x)
      j3 = j2
    i0 = j3
    x = x[:i]+os.environ[x[j1:j2]]+x[j3:]
  return x
      
def median(l):
  if not l: 
    sys.stderr('empty list for median calculation')
    return 0
  l2 = sorted(l)
  n = len(l2)
  if n%2==1: return l2[n/2]
  else: return 0.5*(l2[n/2]+l2[n/2-1])
