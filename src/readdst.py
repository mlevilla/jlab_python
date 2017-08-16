from struct import *
from params import *


fmt0 = [['=IIIfIIii','',['irun','iev','time','Eg','trigger','hastagger','xpos','ypos']],['i','n',['n_adc']],['HH','adc',['adc_val','adc_id']]]

fmt1 = [['=iBQfff','',['iev','trigger','time','xpos','ypos','Eg']],['HHH','n',['n_adc','n_tdc','n_tag']],['HH','adc',['adc_id','adc_val']],['HH','tdc',['tdc_id','tdc_val']],['HH','tag',['tag_id','tag_val']]]

fmt2 = [['=iBQfff','',['iev','trigger','time','xpos','ypos','ebeam']],['HHHH','n',['n_adc','n_tdc','n_tag','n_et']],['HH','adc',['adc_id','adc_val']],['HH','tdc',['tdc_id','tdc_val']],['HH','tag',['tag_id','tag_val']],['iff','et',['et_id','et_t','et_e']]]

fmt3 = [['=IBQf','',['iev','trigger','time','Eb']],['HH','n',['n_adc','n_tdc']],['HH','adc',['adc_id','adc_val']],['HH','tdc',['tdc_id','tdc_val']]]

fmt_calib_prod = [['=IBf','',['iev','trigger','Eb']],['H','n',['n_adc']],['HH','adc',['adc_id','adc_val']]]

class dst:
  def __init__(self,path='',fmt=fmt2,lload=[],step=100000,n=-1,run=0,folder=replayf):
    if path=='' and run==0: 
      print 'no path or run given'
      return
    elif path=='': path = folder+'/prad_'+str(run).zfill(6)+'.dst'
    self.f = open(path,'rb')
    self.e = event()
    self.i = -1
    self.fmt = fmt
    self.lload = lload
    self.step = step
    self.n = n
    self.start = 0
    self.d = {}
    self.cs = {}
    for x in fmt:
      for y in x[2]: self.d[y] = []
    self.read_header()
    self.read_footer()

  def __repr__(self): return 'DST(run={0}, entries={1})'.format(self.run,self.n)
    
  def read_header(self):
    if self.fmt!=fmt0: self.run = unpack('I',self.f.read(4))[0]
    self.ancre = {0:self.f.tell()}

  def read_footer(self):
    if self.fmt!=fmt0:
      self.f.seek(-4,2)
      self.n = unpack('I',self.f.read())[0]
      self.f.seek(self.ancre[0])
     
  def __iter__(self):
    #self.goto(self.start-1)
    return self

  def itr(self,start=0,lload=[]):
    self.start = start
    self.lload = lload
    return self.__iter__()

  def next(self):
    if self.end(): raise StopIteration
    else:
      self.new_event()
      return self.e

  def next_cond(self,lload=[],**cond):
    self.read([''])
    b = all([self.e[x]==y for x,y in cond.items()])
    while not b:
      self.read([''])
      b = all([self.e[x]==y for x,y in cond.items()])
    for x in lload: self.load(x)
    return self.e

  def end(self): return self.i==self.n-1

  def goto(self,n):
    min_ancre = min(self.ancre.keys(),key=lambda x: n/self.step-x if n/self.step>x else float('inf'))
    if n<self.i or min_ancre>self.i/self.step:
      self.f.seek(self.ancre[min_ancre])
      self.i = min_ancre*self.step-1
    else: self.f.seek(self.cs['start'],0)
    for i in range(n-self.i-1): self.read()
    if n!=-1: self.new_event()
    else: self.e = event()

  def __len__(self):
    if self.n!=-1: return self.n-self.i-1
    start = self.f.tell()
    istart = self.i
    self.f.seek(0,0)
    self.i = -1
    while not self.end(): self.skip()
    self.n = self.i+1
    self.f.seek(start,0)
    self.i = istart
    return self.n-self.i-1

  def close(self): self.f.close()

  def read(self,lload=[]):
    self.cs['start'] = self.f.tell()
    for x in self.fmt:
      if (x[1]=='n') or (x[1]=='' and ('' in lload)):
        buff = unpack(x[0],self.f.read(calcsize(x[0])))
        for i,y in enumerate(x[2]): self.e[y] = buff[i]
      elif x[1] in lload:
        for y in x[2]: self.e[y] = []
        buff = unpack(self.e['n_'+x[1]]*x[0],self.f.read(self.e['n_'+x[1]]*calcsize(x[0])))
        for j in range(self.e['n_'+x[1]]):
          for i,y in enumerate(x[2]): self.e[y].append(buff[i+j*len(x[2])])
      else:
        self.cs[x[1]] = self.f.tell()
        size = 1 if x[1]=='' else self.e['n_'+x[1]]
        self.f.seek(calcsize(x[0])*size,1)
    self.cs['end'] = self.f.tell()
    self.i += 1
    if self.i%self.step==0: self.ancre[self.i/self.step] = self.cs['start']

  def load(self,name):
    self.f.seek(self.cs[name],0)
    x = self.fmt[[y[1] for y in self.fmt].index(name)]
    if name=='':
      buff = unpack(x[0],self.f.read(calcsize(x[0])))
      for i,y in enumerate(x[2]): self.e[y] = buff[i]
    else:
      for y in x[2]: self.e[y] = []
      buff = unpack(self.e['n_'+x[1]]*x[0],self.f.read(self.e['n_'+x[1]]*calcsize(x[0])))
      for j in range(self.e['n_'+x[1]]):
        for i,y in enumerate(x[2]): self.e[y].append(buff[i+j*len(x[2])])
    self.f.seek(self.cs['end'],0)

  def eclear(self):
    del self.e
    self.e = event()

  def new_event(self):
    self.eclear()
    self.read(self.lload)

  def store(self,lstore=[]):
    for x,y in self.e.items():
      if lstore==[] or (x in lstore): self.d[x].append(y)

  #def draw(self,name):
  #  return hdraw(l=self.d[name])
    


    
class event(dict):
  def __init__(self): super(event,self).__init__()   
  def __getattr__(self,name):
    if name in self: return self[name]
    else: raise AttributeError("No such attribute: " + name)
  def __setattr__(self,name,value): self[name] = value
  def __delattr__(self,name):
    if name in self: del self[name]
    else: raise AttributeError("No such attribute: " + name)


runs_ancre = {}  


fmt_chao = [['=iBBQ','',['iev','typ','trigger','time']],['I','n',['n_adc']],['HH','adc',['adc_id','adc_val']],['I','n',['n_tdc']],['HH','tdc',['tdc_id','tdc_val']],['I','n',['n_gem']],['','gem',[['BBB','',['addr0','addr1','addr2']],['I','n',['n_val']],['f','val',['gem_val']]]],['I','n',['n_dsc']],['II','dsc',['gated','ungated']]]  
    
class chao:
  def __init__(self,path='',run=0,lload=[],step=100000,folder=replayf):
    if path=='' and run==0: 
      print 'no path or run given'
      return
    elif path=='': path = folder+'/prad_'+str(run).zfill(6)+'.dst'
    self.f = open(path,'rb')
    self.read_header()
    self.e = event()
    self.epics = event()
    self.i = -1
    self.lload = lload
    self.isepics = 0
    self.step = step
    self.cs = {}
    self.n = -1
    if self.run in runs_length: self.n = runs_length[self.run]
    if self.run in runs_ancre: 
      for i,x in enumerate(runs_ancre[self.run]): self.ancre[i] = x

  def read_header(self):
    self.f.seek(0,0)
    self.run_info()
    self.hycal_info()
    self.gem_info()
    self.epics_info()
    self.ancre = {0:self.f.tell()}

  def run_info(self):
    self.f.seek(4,0)
    if self.f.read(4)=='\x03\xe0\xe0\xe0':
      buff = unpack('iddd',self.f.read(calcsize('iddd')))
      self.run = buff[0]
      self.beam_charge = buff[1]
      self.dead_count = buff[2]
      self.ungated_count = buff[3]
    else: 
      self.run = -1
      self.f.seek(4,0)

  def hycal_info(self):
    if self.f.read(4)=='\x04\xe0\xe0\xe0':
      size = unpack('i',self.f.read(calcsize('i')))[0]
      if size>1733: size = unpack('i',self.f.read(calcsize('i')))[0]
      self.ped_mean, self.ped_sigma, self.factor, self.base_factor, self.base_gain = [],[],[],[],[]
      for i in range(size):
        buff = unpack('=ddddi',self.f.read(calcsize('=ddddi')))
        self.ped_mean.append(buff[0])
        self.ped_sigma.append(buff[1])
        self.factor.append(buff[2])
        self.base_factor.append(buff[3])
        self.base_gain.append(list(unpack('d'*buff[4],self.f.read(calcsize('d')*buff[4]))))
    else: self.f.seek(-4,1)

  def  gem_info(self):
    x = self.f.read(4)
    if x=='\x05\xe0\xe0\xe0': 
      size = unpack('i',self.f.read(calcsize('i')))[0]
    elif x=='\x00\x05\xe0\xe0':
      self.f.seek(5,1)
    else:
      self.f.seek(-4,1)
      return 
    size = unpack('i',self.f.read(calcsize('i')))[0]
    self.gem_addr0, self.gem_addr1, self.gem_ped = [],[],[]
    for i in range(size):
      buff = unpack('iii',self.f.read(calcsize('iii')))
      self.gem_addr0.append(buff[0])
      self.gem_addr1.append(buff[1])
      self.gem_ped.append(list(unpack('ff'*buff[2],self.f.read(calcsize('ff')*buff[2]))))
      
  def epics_info(self):
    x = self.f.read(4)
    if x=='\x02\xe0\xe0\xe0': 
      size = unpack('i',self.f.read(calcsize('i')))[0]
    elif x=='\x0e\x02\xe0\xe0':
      self.f.seek(5,1)
    else:
      self.f.seek(-4,1)
      return 
    size = unpack('i',self.f.read(calcsize('i')))[0]
    self.epics_name, self.epics_value = [],[]
    for i in range(size):
      name_size = unpack('i',self.f.read(calcsize('i')))[0]
      buff = unpack('='+str(name_size)+'sif',self.f.read(calcsize('='+str(name_size)+'sif')))
      self.epics_name.append(buff[0])
      self.epics_value.append(buff[2])
    self.epics_iev = 0

  def __iter__(self): return self

  def close(self): self.f.close()

  def eclear(self):
    del self.e
    self.e = event()

  def new_event(self):
    self.eclear()
    self.read(self.lload)

  def next(self):
    if self.end(): raise StopIteration
    else:
      self.new_event()
      return self.e

  def end(self):
    b = unpack('I',self.f.read(4))[0]==0xe0e0e003
    self.f.seek(-4,1)
    return b

  def goto(self,n):
    min_ancre = min(self.ancre.keys(),key=lambda x: n/self.step-x if n/self.step>x else float('inf'))
    if n<self.i or min_ancre>self.i/self.step:
      self.f.seek(self.ancre[min_ancre])
      self.i = min_ancre*self.step-1
    else: self.f.seek(self.cs['start'],0)
    for i in range(n-self.i-1): self.read()
    if n!=-1: self.new_event()
    else: self.e = event()

  def read(self,lload=[]):
    self.cs['start'] = self.f.tell()
    if unpack('I',self.f.read(4))[0]==0xe0e0e001:
      self.isepics = 1
      self.read_epics()
      return
    self.f.seek(-4,1)
    if unpack('I',self.f.read(4))[0]!=0xe0e0e000: return
    self.isepics = 0
    for x in fmt_chao:
      if (x[1]=='n') or (x[1]=='' and ('' in lload)):
        buff = unpack(x[0],self.f.read(calcsize(x[0])))
        for i,y in enumerate(x[2]): self.e[y] = buff[i]
      elif (x[1] in lload) and x[1]!='gem': 
        for y in x[2]: self.e[y] = []
        buff = unpack(self.e['n_'+x[1]]*x[0],self.f.read(self.e['n_'+x[1]]*calcsize(x[0])))
        for j in range(self.e['n_'+x[1]]):
          for i,y in enumerate(x[2]): self.e[y].append(buff[i+j*len(x[2])])
      elif (x[1] in lload) and x[1]=='gem':
        for y in x[2]:
          for z in y[2]: self.e[z] = []
        for i in range(self.e['n_'+x[1]]):
          for y in x[2]:
            if y[1] in ['','n']:
              buff = unpack(y[0],self.f.read(calcsize(y[0])))
              for j,z in enumerate(y[2]): self.e[z].append(buff[j])
            else:
              for z in y[2]: self.e[z].append([])
              print y,self.e['n_'+y[1]]
              buff = unpack(self.e['n_'+y[1]][0]*y[0],self.f.read(self.e['n_'+y[1]][0]*calcsize(y[0])))
              for j in range(self.e['n_'+y[1]][0]):
                for k,z in enumerate(y[2]): self.e[z][-1].append(buff[k+j*len(y[2])])
      elif x[1]=='gem':
        self.cs[x[1]] = self.f.tell()
        #self.f.seek(calcsize('=BBBIfff')*self.e['n_'+x[1]],1)
        for i in range(self.e['n_'+x[1]]):
          for y in x[2]:
            if y[1]=='': self.f.seek(calcsize(y[0]),1)
            elif y[1]=='n':
              buff = unpack(y[0],self.f.read(calcsize(y[0])))
              for j,z in enumerate(y[2]): self.e[z] = buff[j]
            else: self.f.seek(calcsize(y[0])*self.e['n_'+y[1]],1)
      else:
        self.cs[x[1]] = self.f.tell()
        size = 1 if x[1]=='' else self.e['n_'+x[1]]
        self.f.seek(calcsize(x[0])*size,1)
    self.cs['end'] = self.f.tell()
    self.i += 1
    if self.i%self.step==0:self.ancre[self.i/self.step] = self.cs['start']

  def read_epics(self):
    self.epics['iev'] = unpack('i',self.f.read(calcsize('i')))[0]
    n = unpack('i',self.f.read(calcsize('i')))[0]
    buff = unpack(n*'f',self.f.read(calcsize(n*'f')))
    for i in range(n): self.epics[self.epics_name[i]] = buff[i]
    self.i+=1

  def load(self,name):
    self.f.seek(self.cs[name],0)
    x = fmt_chao[[y[1] for y in fmt_chao].index(name)]
    if name=='':
      buff = unpack(x[0],self.f.read(calcsize(x[0])))
      for i,y in enumerate(x[2]): self.e[y] = buff[i]
    elif name!='gem':
      for y in x[2]: self.e[y] = []
      buff = unpack(self.e['n_'+x[1]]*x[0],self.f.read(self.e['n_'+x[1]]*calcsize(x[0])))
      for j in range(self.e['n_'+x[1]]):
        for i,y in enumerate(x[2]): self.e[y].append(buff[i+j*len(x[2])])
    else:
      for y in x[2]:
        for z in y[2]: self.e[z] = []
      for i in range(self.e['n_'+x[1]]):
        for y in x[2]:
          if y[1] in ['','n']:
            buff = unpack(y[0],self.f.read(calcsize(y[0])))
            for j,z in enumerate(y[2]): self.e[z].append(buff[j])
          else:
            for z in y[2]: self.e[z].append([])
            buff = unpack(self.e['n_'+y[1]]*y[0],self.f.read(self.e['n_'+y[1]]*calcsize(y[0])))
            for j in range(self.e['n_'+y[1]]):
              for k,z in enumerate(y[2]): self.e[z][-1].append(buff[k+j*len(y[2])])
    self.f.seek(self.cs['end'],0)

  def __len__(self): 
    if self.n!=-1: return self.n-self.i-1
    else: -1



    
def read_vector(f,ft,tr=False):
  n = unpack('I',f.read(4))[0]
  m = len(ft)
  l = list(unpack(n*ft,f.read(n*calcsize(ft))))
  if m==1: return l
  else: 
    if not tr: return [[l[m*i+j] for j in range(m)] for i in range(n)]
    else: return [[l[m*i+j] for i in range(n)] for j in range(m)]
  
load_runs_length()

class chao2:
  def __init__(self,path='',run=0,folder=replayf+'/event_sel'):
    if path=='' and run==0: 
      print 'no path or run given'
      return
    elif path=='': path = folder+'/prad_'+str(run).zfill(6)+'_sel.dst'
    self.f = open(path,'rb')
    self.epics = event()
    self.e = event()
    self.i = 0
    self.n = -1
    if 'sel' in path and run in runs_length_sel: self.n = runs_length_sel[run]
    elif 'raw' in path and run in runs_length_raw: self.n = runs_length_raw[run]
    self.eheader = []
    self.read_header()
    self.load_epics()

  def read_header(self):
    self.f.seek(0,0)
    self.header = list(unpack('HHIQ',self.f.read(16)))

  def load_epics(self):
    self.epics_names = readfile(os.environ['PRAD_PATH']+'/config/epics_channels.conf')

  def __iter__(self): return self

  def close(self): self.f.close()

  def next(self):
    if self.end(): 
      raise StopIteration
    else:
      self.read()
      return self.e

  def goto(self,n):
    while not self.end() and self.i<n-1: 
      self.read(0)
    while not self.end() and self.i<n: 
      self.read()

  def end(self):
    return self.f.tell()>=self.header[3]

  def read(self,flag=1):
    self.eheader = list(unpack('HHI',self.f.read(8)))
    if flag==0: 
      self.f.seek(self.eheader[2],1)
      if self.eheader[1]==0:
        self.i+=1
    else:
      if self.eheader[1]==0:
        self.i+=1
        self.read_event()
      elif self.eheader[1]==1:
        self.read_epics()

  def read_epics(self):
    self.epics.iev = unpack('i',self.f.read(4))[0]
    l = read_vector(self.f,'f')
    for x,y in zip(self.epics_names,l): self.epics[x] = y

  def read_event(self):
    l = list(unpack('=iBBQ',self.f.read(14)))
    for x,y in zip(['iev','typ','trigger','time'],l): self.e[x] = y
    [self.e.adc_id,self.e.adc_val] = read_vector(self.f,'HH',1)
    [self.e.tdc_id,self.e.tdc_val] = read_vector(self.f,'HH',1)
    ngem = unpack('I',self.f.read(4))[0]
    self.e.gem_id,self.e.gem_val = [],[]
    for i in range(ngem):
      self.e.gem_id.append(list(unpack('BBB',self.f.read(3))))
      self.e.gem_val.append(read_vector(self.f,'f'))
    self.e.dsc = read_vector(self.f,'II')
    
  def __len__(self): return self.n
