from params import * 

edge = get_regions('inner')
erange_1GeV = [50,250,750,1500]
erange_2GeV = [100,500,1500,2500]
regions = ['lg_out','lg_in','pwo_out','pwo_center','pwo_in']
deads = [1,126,169,415,486,732,775,900,1213,1230,1245,1266,1322,1325,1832,1835,1891,1912,1927,1944]

load_neighbors()
load_names()
module_names2 = {x:y for x,y in module_names.items() if x<181 or (x<721 and (x-1)%30>23) or (1000<x and (x-1001)%34<17 and (x-1001)/34<17)}

index_map = readfile(conff+'/index_module.txt',[int,int],out=dict)

def module2index(cid):
  lg,row,col = rowcolid(cid)
  sector = which_sector(cid,lg,row,col)
  if sector==1: return row/3*8+col/3
  if sector==2: return 16+(col-24)/3*8+row/3
  if sector==3: return (29-row)/3*8+(29-col)/3
  if sector==4: return 16+(5-col)/3*8+(29-row)/3
  if sector==0:
    if col<17 and row<17: return 32+row/2*9+col/2
    if row<17: return 32+row/2*9+(33-col)/2
    if col<17: return 32+(33-row)/2*9+col/2
    return 32+(33-row)/2*9+(33-col)/2
  return -1

def module2index_asym(cid):
  lg,row,col = rowcolid(cid)
  sector = which_sector(cid,lg,row,col)
  if sector==1: return row/3*8+col/3
  if sector==2: return 16+(col-24)/3*8+row/3
  if sector==3: return 32+(29-row)/3*8+(29-col)/3
  if sector==4: return 48+(5-col)/3*8+(29-row)/3
  if sector==0:
    return 64+row/2*17+col/2
    if row<17: return 64+row/2*9+(33-col)/2
    if col<17: return 64+(33-row)/2*9+col/2
    return 64+(33-row)/2*9+(33-col)/2
  return -1

def module2index2(cid):
  lg,row,col = rowcolid(cid)
  if lg and (row<6 or col>23) and row<24: return cid
  elif lg: return 30*(29-row)+(29-col)+1
  elif row<17 and col<17: return cid
  elif col<17: return 34*(33-row)+col+1001
  elif row<17: return 34*row+(33-col)+1001
  else: return 34*(33-row)+(33-col)+1001


def correct_module_id(cid,x,y):
  tx = (x-module_pos[cid][0])/cell_size[int(cid<1000)][0]
  ty = (y-module_pos[cid][1])/cell_size[int(cid<1000)][1]
  cid0 = cid
  if abs(tx)<0.5 and abs(ty)<0.5: return cid,tx,ty
  cid = which_module(x,y,neighbors[cid])
  if cid==-1: return cid0,tx,ty
  tx = (x-module_pos[cid][0])/cell_size[int(cid<1000)][0]
  ty = (y-module_pos[cid][1])/cell_size[int(cid<1000)][1]
  return cid,tx,ty  

def regions_id(j):
  if 0<=j<=7 or 24<=j<=31: i2 = 0
  elif 8<=j<=23: i2 = 1
  elif 32<=j<=40 or (j-32)%9==0: i2 = 2
  elif j in [102,103,111]: i2 = 4
  else: i2 = 3
  return i2
  
