from misc import *
from jsub import *

time_start = time.time()

# folders
homef = os.environ['HOME']
workf = os.environ['WORK']
pythonf = homef+'/python'
conff = pythonf+'/config'
dataf = pythonf+'/data'
execf = pythonf+'/exec'
srcf = pythonf+'/src'
replayfm = workf+'/replay'
jobf = workf+'/jobs'
replayf = '/work/hallb/prad/replay/event_sel'
pradf = os.environ['PRAD_PATH']
pradsimf = os.environ['PRAD_SIM']


# hycal infos
cell_size = [[20.77,20.75],[38.15,38.15]]
  
def get_names(t='dict'):
  if t=='primex': l = ['G'+str(i) for i in range(1001)]+['W'+str(i) for i in range(1,1157)]
  elif t=='dict': l = {int(x[0]=='W')*1000+int(x[1:]):x for x in readfile(conff+'/module_names_prim.txt')}
  else: l = readfile(conff+'/module_names_'+t+'.txt')
  return l

module_names = {}
def load_names():
  global module_names
  for x,y in get_names().items(): module_names[x] = y

def get_positions():
  l = readfile(conff+'/module_coord.txt',[str,float,float])
  return {int(x[0][0]=='W')*1000+int(x[0][1:]):x[1:] for x in l}

module_pos,module_theta = {},{}
def load_positions(z=5722.):
  global module_pos, module_theta
  for x,y in get_positions().items(): 
    module_pos[x] = y
    module_theta[x] = degrees(atan2((y[0]**2+y[1]**2)**0.5,z)) if x>1000 else degrees(atan2((y[0]**2+y[1]**2)**0.5,z-100))

def get_neighbors():
  l = readfile(conff+'/module_coord.txt',[str]+[int]*8,cols=[0]+range(3,11))
  return {int(x[0][0]=='W')*1000+int(x[0][1:]):[y for y in x[1:] if y!=-1] for x in l}

neighbors = {}
def load_neighbors():
  global neighbors
  for x,y in get_neighbors().items(): neighbors[x] = y

hycal_trgeff = {}
def load_hycal_trgeff():
  global hycal_trgeff
  l = readfile(os.environ['PRAD_PATH']+'/database/hycal_trgeff.dat',[str,float],start=1,out=dict)
  for x,y in l.items():
    hycal_trgeff[int(x[0]=='W')*1000+int(x[1:])] = y

live_charge = {}
def load_live_charge():
  global live_charge
  l = readfile(conff+'/beam_charge.dat',[int,float],cols=[0,3],start=1,out=dict)
  for x,y in l.items(): live_charge[x] = y

thickness, dthickness = {},{}
def load_thickness():
  global thickness, dthickness
  l = readfile(conff+'/target_info.txt',[int,float,float],cols=[0,5,6],out=dict,start=1)
  for x,y in l.items():
    thickness[x] = y[0]
    dthickness[x] = y[1]

def get_regions(name):
  if name=='trans_pwo': return range(1001,1069)+range(2089,2157)+[34*i+1001 for i in range(2,32)]+[34*i+1002 for i in range(2,32)]+[34*i+1033 for i in range(2,32)]+[34*i+1034 for i in range(2,32)]
  elif name=='trans_lg': return  range(156,176)+range(726,746)+[186+30*i for i in range(18)]+[205+30*i for i in range(18)]
  elif name=='trans':
    return get_regions('trans_lg')+get_regions('trans_pwo')
  elif name=='inner': return range(1526,1530)+range(1628,1632)+[1560,1563,1594,1597]
  elif name=='inner0': return [1527,1528,1560,1563,1594,1597,1629,1630]
  elif name=='inner2': return range(1526,1530)+range(1628,1632)+[1560,1563,1594,1597]+range(1491,1497)+range(1661,1667)+[1525+34*i for i in range(4)]+[1530+34*i for i in range(4)]
  elif name=='outer': return range(1,31)+range(871,901)+[31+30*i for i in range(28)]+[60+30*i for i in range(28)]
  elif name=='outer2': return  range(1,61)+range(841,901)+[61+30*i for i in range(26)]+[62+30*i for i in range(26)]+[89+30*i for i in range(26)]+[90+30*i for i in range(26)]
  elif name=='center_pwo': return list(set(range(1001,2157))-set(get_regions('trans_pwo'))-set(get_regions('inner2')))
  elif name=='center_lg': return list(set(range(901))-set(get_regions('trans_lg'))-set(get_regions('outer2')))

def is_in_sector(cid,sector=set([0,1,2,3,4])):
  if cid>1000: return (0 in sector)
  if cid<=180 and (cid-1)%30<24: return (1 in sector) # top
  if cid<=720 and (cid-1)%30>=24: return (2 in sector) # right
  if cid>720 and (cid-1)%30>=6: return (3 in sector) # bottom
  if cid>180 and (cid-1)%30<6: return (4 in sector) # left
  return False

def is_in_region(cid,region):
  if region=='lg': return cid<1000
  if region=='pwo': return cid>1000
  if region=='transition_lg': return 156<=cid<=175 or 726<=cid<=745 or (180<cid<720 and (cid%30==6 or cid%30==25))
  if region=='transition_pwo': return 1001<=cid<=1034 or 2123<=cid<=2156 or (1034<cid<2157 and ((cid-1000)%34==0 or (cid-1000)%34==1))
  if region=='transition': return 156<=cid<=175 or 726<=cid<=745 or (180<cid<720 and (cid%30==6 or cid%30==25)) or 1001<=cid<=1034 or 1123<=cid<=1156 or (1034<cid<1157 and (cid%34==0 or cid%34==1))

def which_module(x,y):
  if len(module_pos)==0:
    sys.sterr.write('positions not loaded')
    return -1
  for i,pos in module_pos.items():
    lg = int(i<1000)
    if abs(x-pos[0])<cell_size[lg][0]/2. and abs(y-pos[1])<cell_size[lg][1]/2.:
      return i
  return -1


# periods
def get_periods(t):
  periods_gem_calib  = [[[793,794,795,796,847,849,850,851,852,853,854,855,982,983,987,988]],[[],[],[]],[[],[],[],[],[]],[[]],[[],[],[]]]
  periods_calib = [[[889,890,893,894,895,896,916,918,919,923,924,925,926,927,928,929,932,933,934,935,946,947,948,949,950,951,952,953,955,956,957,958,960,961,962,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979]],[[],[],[]],[[],[],[],[],[]],[[]],[[],[],[]]]
  #removed 1287, 1301
  periods_prod = [[[]],[[1059,1062,1071,1072,1075,1076,1077,1079,1081,1082],[1083,1085,1086,1087,1092,1093,1094,1096,1099,1100,1101,1104,1105,1106],[1136,1137,1155,1158,1159,1160,1201,1202,1203,1205,1206,1207]],[[1238,1240,1241,1242,1243,1287],[1288,1290,1291,1292,1293,1295,1296],[1302,1303,1304,1307,1308,1309,1310,1311,1313],[1314,1315,1316,1319,1320,1322,1323,1325,1328],[1331,1332,1334,1336,1337,1338,1340,1341]],[[1362,1370,1371,1372,1373,1374,1382,1383,1384,1385,1394,1395,1396,1397,1400,1401,1404,1405]],[[1412,1413,1415,1416,1420,1421,1422,1423,1424,1426,1427,1429,1430,1443,1444,1445,1446,1449,1451,1452,1453],[1455,1456,1457,1458,1460,1461,1462,1463,1466,1467,1468,1470,1471,1472,1474,1475,1476,1478,1479,1480,1482,1484,1485,1487,1488,1489,1491,1492,1495],[1506,1507,1508,1510,1511,1512,1513,1516]]]
  #removed 1153, 1197, 1318, 1386
  periods_empty = [[[]],[[1080],[1084,1088,1103],[1127,1129,1135,1153,1154,1161,1163,1198]],[[1237],[1289,1294],[1306,1312],[1317,1324,1326,1327],[1329,1330,1339]],[[1369,1381,1387,1399]],[[1406,1411,1418,1425,1431,1448,1454],[1459,1465,1469,1473,1477,1481,1486,1490,1496],[1509,1514,1515]]]
  periods_total_empty = [[[]],[[],[1091,1123],[1125,1126,1147]],[[],[],[],[],[1342]],[[]],[[1417],[],[]]]
  periods_carbon = [[[]],[[],[],[]],[[],[],[],[],[1345]],[[]],[[],[1498,1501],[]]]
  if t=='calib': return periods_calib
  elif t=='gem_calib': return periods_gem_calib
  elif t=='prod': return periods_prod
  elif t=='empty': return periods_empty
  elif t=='total_empty': return periods_total_empty
  elif t=='carbon': return periods_carbon
  elif t=='all': 
    p = [[sorted(u1+v1+a1+b1+c1+d1) for u1,v1,a1,b1,c1,d1 in zip(u2,v2,a2,b2,c2,d2)] for u2,v2,a2,b2,c2,d2 in zip(periods_gem_calib,periods_calib,periods_prod,periods_empty,periods_total_empty,periods_carbon)]
    for x in p:
      for y in x: y.sort()
    return p
    
def get_assoc_empty(lprod,lempty):
  l = []
  for x in lprod:
    dmin,imin = 100,-1
    for i,y in enumerate(lempty):
      if abs(x-y)<dmin:
        dmin = abs(x-y)
        imin = i
    l.append(imin)
  return l
  
# run infos
def get_period_number(run,periods):
  for i,p in enumerate(periods):
    for j,s in enumerate(p):
      for k,r in enumerate(s):
        if run==r: return [i,j,k]
        if run<r: break
  return [-1,-1,-1]

runs_length_sel, runs_length_raw = {},{}
def load_runs_length():
  global runs_length_sel,runs_length_raw
  for x,y in readfile(conff+'/runs_sel_length.txt',[int]*2,out=dict).items(): runs_length_sel[x]=y
  for x,y in readfile(conff+'/runs_raw_length.txt',[int]*2,out=dict).items(): runs_length_raw[x]=y
  

# bad events
bad_events = {}
def load_bad_events():
  global bad_events
  for x,y in readfile(conff+'/bad_events.txt',[int,int,int],out=dict).items(): 
    if x in bad_events: bad_events[x].append(y)
    else: bad_events[x] = [y]

def event_is_bad(iev,run):
  if run not in bad_events: return False
  for x in bad_events[run]:
    if x[0]<=iev<=x[1]: return True
  return False


def correct_w_time(l,time):
  i = 0
  while i<len(l) and time>=l[i][0]: i+=1
  if i==len(l): return [l[-1][1],l[-1][2]]
  return [l[i-1][j]+(l[i][j]-l[i-1][j])/(l[i][0]-l[i-1][0])*(time-l[i-1][0]) for j in [1,2]]
    
pedestal,pedestal_sigma,lms_factors,lms_gains,linfactor,ecalib,calib_gains = {},{},{},{},{},{},{}
beam_energy,hycal_center,gem1_center,gem2_center = [0.],[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]

def load_run_params(run,fcalib=pradf+'/database/calibration',finfo=pradf+'/database/baseinfo'):
  [i,j,_] = get_period_number(run,get_periods('all'))
  [names,ped,dped,lms] = readfile(finfo+'/db_prad_baseinfo_'+str(run)+'.dat',[str,float,float,float],start=1,tr=True)
  ids = [int(x[1:])+1000*(x[0]=='W') for x in names if x[:3]!='LMS']
  [names2,gains,e,a,_,_,lms3] = readfile(fcalib+'/calibration_'+str(i)+'_'+str(j+1)+'.dat',[str,float,float,float,float,float,float],tr=True)
  ids2 = [int(x[1:])+1000*(x[0]=='W') for x in names2]
  hycal_centers = readfile(conff+'/hycal_center.txt',[int,float,float,float],out=dict)
  gem1_centers = readfile(conff+'/gem1_center.txt',[int,float,float,float],out=dict)
  gem2_centers = readfile(conff+'/gem2_center.txt',[int,float,float,float],out=dict)
  global pedestal,pedestal_sigma,lms_factors,lms_gains,calib_gains,ecalib,linfactor, beam_energy, hycal_center, gem1_center, gem2_center
  for w,x,y,z in zip(ids,ped,dped,lms): 
    pedestal[w] = x
    pedestal_sigma[w] = y
    lms_factors[w] = z
  for v,w,x,y,z in zip(ids2,gains,e,a,lms3):
    calib_gains[v] = w
    ecalib[v] = x
    linfactor[v] = y
    lms_gains = z/lms_factors[v] if lms_factors[v]!=0 else 0
  if run<1362: beam_energy[0] = 1100.44 # 1099.65
  else: beam_energy[0] = 2142. # 2140.56
  for i in range(3): 
    hycal_center[i] = hycal_centers[run][i]
    gem1_center[i] = gem1_centers[run][i]
    gem2_center[i] = gem2_centers[run][i]
  
tagger_timings = {1:[0,0,0],2:[0,0,0],5:[0,0,0]}
def load_tagger_timings(run):
  d = readfile(conff+'/tagger_timings.txt',[int]+[float]*9,out=dict)
  global tagger_timings
  if run not in d: 
    print 'no tagger timing information for run', run
  else:
    for i in range(3): 
      tagger_timings[1][i] = d[run][i]
      tagger_timings[2][i] = d[run][i+3]
      tagger_timings[5][i] = d[run][i+6]


def get_runs_between(run_min,run_max,t='all'):
  return [r for r in unnest(get_periods(t)) if run_min<=r<=run_max]

def get_periods_between(run_min,run_max,t='all'):
  l = get_periods(t)
  periods = []
  for period in l:
    for subperiod in period:
      x = [r for r in subperiod if run_min<=r<=run_max]
      if x: periods.append(x)
  return periods

def get_bg_runs(run1=None,run2=None):
  l = [[[1059,1062,1071,1072,1075,1076,1077,1079],[1080]],[[1081,1082,1083],[1080,1084]],[[1085,1086,1087],[1084,1088]],[[1092,1093,1094,1096,1099,1100,1101],[1088,1103]],[[1104,1105,1106],[1103]],[[1136,1137],[1127,1129,1135]],[[1155,1158,1159,1160],[1154,1161,1163]],[[1201,1202,1203,1205,1206,1207],[1198]],[[1238,1240,1241,1242,1243],[1237]],[[1287,1288,1290,1291,1292,1293],[1289,1294]],[[1295,1296,1302,1303,1304],[1294,1306]],[[1307,1308,1309,1310,1311],[1306,1312]],[[1313,1314,1315,1316],[1312,1317]],[[1319,1320,1322,1323,1325,1328],[1317,1324,1326,1327]],[[1331,1332,1334,1336,1337,1338,1340,1341],[1329,1330,1339]],[[1362],[1369]],[[1370,1371,1372,1373,1374],[1381]],[[1382,1383,1384,1385],[1381,1387]],[[1394,1395,1396,1397],[1399]],[[1400,1401,1404,1405],[1399,1406]],[[1412,1413,1415,1416],[1411,1418]],[[1420,1421,1422,1423,1424],[1418,1425]],[[1426,1427,1429,1430],[1425,1431]],[[1443,1444,1445,1446],[1431,1448]],[[1449,1451,1452,1453],[1448,1454]],[[1455,1456,1457,1458],[1454,1459]],[[1460,1461,1462,1463],[1459,1465]],[[1466,1467,1468],[1465,1469]],[[1470,1471,1472],[1469,1473]],[[1474,1475,1476],[1473,1477]],[[1478,1479,1480],[1477,1481]],[[1482,1484,1485],[1481,1486]],[[1487,1488,1489],[1486,1490]],[[1491,1492,1495],[1490,1496]],[[1506,1507,1508],[1496,1509]],[[1510,1511,1512,1513],[1509,1514]],[[1516],[1515]]]
  if run1 is None: return l
  elif run2 is None:
    for x in l:
      if run1 in x[0]: return x[1]
  else:
    l2 = [[[y for y in x[0] if run1<=y<=run2],x[1]] for x in l]
    l2 = [x for x in l2 if x[0]!=[]]
    return l2
