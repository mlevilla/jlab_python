import sys
from mdict import Mdict

largs0 = [[('sp',0,'','show_progress'),
           ('b',False,'','batch_mode'),
           ('wt','120','','wall_time')],
          [('iv',False,'','interactive_mode')]]
show_progress, batch, interactive, wtime = [0], [False], [False], ['120']


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

class ArgParser:
  
  def __init__(self,largs=[],argv=[],empty=1,lmax=15,sub=0):
    self.largs = largs
    self.argv = []
    self.argd = Mdict()
    self.argl = []
    self.imax = 0
    self.rest = []
    
    self.read_argv(argv)
    self.catch_args()
    if sub: self.catch_globs()
    self.print_help(empty,lmax)

  def read_argv(self,argv):
    self.argv = argv if argv else [catch_env_var(x) for x in sys.argv]
    if len(self.argv)>1 and '.in' in self.argv[1]:
      f = open(self.argv[1])
      self.argv = [catch_env_var(x) for x in f.readline().split(' ')] 
      f.close()

  def catch(self,mark,default,opt=''):
    if not isinstance(default,list): default = [default]
    r = [default] if opt=='several' else default
    if '-'+mark in self.argv:
      i = self.argv.index('-'+mark)
      if opt=='':
        if isinstance(default[0],bool):
          r = [not default[0]]
          self.imax = max(self.imax,i)
        else:
          r = [type(default[j])(self.argv[i+j+1]) for j in range(len(default)) if i+j+1<len(self.argv)]
          self.imax = max(self.imax,min(i+len(default),len(self.argv)))
      elif opt=='several':
        k = i+1
        while k<len(self.argv) and (self.argv[k]=='' or self.argv[k][0]!='-'): k+=1
        r = [[type(default[j])(self.argv[i+1+len(default)*n+j]) for n in range((k-i-1)/len(default))] for j in range(len(default))]
        self.imax = max(self.imax,k-1)
      elif opt=='end':
        r = [' '.join(self.argv[i+1:])]
        self.imax = len(self.argv)
    if len(default)==1: r = r[0]
    return r

  def catch_args(self):
    for x in self.largs:
      if x[0] in [y[0] for y in largs0[0]+largs0[1]]:
        print 'warning: standard argument -'+x[0]+' rewritten'
      y = self.catch(*x[:-1])
      self.argd[x[0]] = y
      self.argl.append(y)
    self.rest = self.argv[self.imax+1:]
  
  def catch_globs(self):
    global show_progress, batch, interactive, wtime          
    for x,y in zip(largs0[0]+largs0[1],[show_progress,batch,wtime,interactive]):
      if x[0] in [a[0] for a in self.largs]: continue
      self.largs.append(x)
      y[0] = self.catch(*x[:-1])
      self.argd[x[0]] = y[0]

  def print_help(self,empty=1,lmax=15):
    if '-h' in self.argv or '--help' in self.argv or '--h' in self.argv or (empty and len(self.argv)==1):
      print 'command:',self.argv[0]
      print 'inputs:'
      max_length0 = max(len(x[0]) for x in self.largs)
      max_length1 = max(min(len(str(x[1])),lmax+3)+len(x[3]) for x in self.largs)
      for x,y,z,c in self.largs:
        mark = '\'' if isinstance(y,str) else ''
        if len(str(y))>lmax: y = '...'+str(y)[-15:]
        sys.stdout.write('  '+('-'+x).ljust(max_length0+1)+'\t'+(' < '+c+' = '+mark+str(y)+mark+' > ').ljust(max_length1+9)+('\t mode: '+{'several':'several','end':'end','':''}[z])*int(z!='')+'\n')
      sys.exit()
    
  def delete(self,ldel,join=0):
    for x in self.largs:
      if x[0] not in ldel: continue
      if '-'+x[0] not in self.argv: continue
      i = self.argv.index('-'+x[0])
      del self.argv[i]
      if x[2]=='several':
        while i<len(self.argv) and self.argv[i][0]!='-': 
          del self.argv[i]
      elif x[2]=='end': 
        self.argv = self.argv[:i]
      elif isinstance(x[1],list):
        for x in x[1]: del self.argv[i]
      elif not isinstance(x[1],bool): 
        del self.argv[i]
    if join: return ' '.join(self.argv)
    return self.argv

  def add_list(self,mark,l,write=''):
    if not isinstance(mark,list): mark = [mark]
    idx,arg0 = [],[]
    for x in mark:
      if '-'+x in self.argv:
        i = self.argv.index('-'+x)+1
        while i<len(self.argv) and self.argv[i][0]!='-':
          idx.append(i)
          arg0.append(self.argv[idx[-1]])  
          i+=1
      else:
        self.argv.extend(['-'+x,''])
        idx.append(len(self.argv)-1)
        arg0.append(self.argv[idx[-1]])  
    largs = []
    for x in l:
      for i,a in zip(idx,arg0): self.argv[i] = a+str(x)
      if write!='':
        f = open(write+str(x)+'.in','w')
        f.write(' '.join(self.argv[1:]))
        f.close()
      largs.append(self.argv)
    if write=='': return largs
    return ' '.join([write+str(x)+'.in' for x in l])
