class Counter:
  def __init__(self,binning=[],names=[],version=0,path=''):
    self.version = version
    self.binning = binning
    self.names = [str(x) for x in names]
    self.counts = {x:[0 for i in range(len(binning)+2)] for x in names}
    if path!='': self.read(path)

  def new(name):
    self.names.append(str(name))
    self.counts[name] = [0 for i in range(len(binning)+2)]

  def count(self,name,value):
    if not self.version: return
    if self.version==2 and name not in self.names: self.new(name)
    if not isinstance(value,list): value = [value]
    for v in value:
      if v<self.binning[0]: self.counts[name][0]+=1
      elif v>=self.binning[-1]: self.counts[name][-2]+=1
      else:
        self.counts[name][-1]+=1
        for i in range(len(self.binning)-1):
          if self.binning[i]<=v<self.binning[i+1]:
            self.counts[name][i+1]+=1
            break

  def cut(self,name,cbool,value):
    if not self.version: return not cbool
    if not cbool: return True
    self.count(name,value)
    return False

  def write(self,path):
    if not self.version: return
    f = open(path,'w')
    f.write('binning '+' '.join([str(x) for x in self.binning])+'\n')
    for x in self.names:
      f.write(x+' '+' '.join([str(z) for z in self.counts[x]])+'\n')
    f.close()

  def read(self,path):
    f = open(path)
    l = [x.split(' ') for x in f.readlines()]
    self.binning = [float(x) for x in l[0][1:]]
    for x in l[1:]:
      self.names.append(x[0])
      self.counts[x[0]] = [int(y) for y in x[1:]]
    f.close()
    
