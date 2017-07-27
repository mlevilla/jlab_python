

class Mdict(dict):
  def __init__(self):
    super(Mdict,self).__init__()

  def __getattr__(self,name):
    return self.__getitem__(name)

  def __setattr__(self,name,value):
    self.__setitem__(name,value)
