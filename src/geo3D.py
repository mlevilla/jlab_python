import sys
from math import *

org = [0,0,0]

# any dimension
def diff(a,b): return [xa-xb for xa,xb in zip(a,b)]
def plus(a,b): return [xa+xb for xa,xb in zip(a,b)]
def mult(f,a): return [x*f for x in a]
def scalar(a,b): return sum([xa*xb for xa,xb in zip(a,b)])
def dist(a,o=org): return scalar(diff(a,o),diff(a,o))**0.5
# 3D
def product(a,b): return [a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3)]
def mixte(a,b,c): return scalar(product(a,b),c)

def proj(a,p,o=org,w=False):
  # intersection between straight line OA and plane P
  # (P) p[0]*x+p[1]*y+p[2]*z + p[3] = 0
  v = diff(a,o)
  if scalar(v,p[:3])==0:
    if w: sys.stderr.write('vector parallel to plane\n')
    return [None,None,None]
  for i in range(3):
    if p[i]!=0:
      h = [-float(p[3])/p[i] if j==i else 0 for j in range(3)]
      break
  u = diff(h,o)
  return plus(mult(scalar(u,p[:3])/scalar(v,p[:3]),v),o)

def projp(a,p): 
  # perpendicular projection of point A over plane P
  return proj(plus(a,p[:3]),p,a)

def proj2d(a,p,o=org):
  # coordinate of A in standard P reference
  c = projp(o,p)
  if p[1]!=0 or p[2]!=0: b = projp([1.,0,0],p)
  else: b = projp([0,1.,0],p)
  ux = diff(b,c)
  ux = mult(1./dist(ux),ux)
  uy = product(p[:3],ux)
  a2 = diff(projp(a,p),c)
  return [scalar(a2,ux),scalar(a2,uy)]

def plane(a,b,c):
  # plane defined by A,B and C
  m = [a,b,c]
  x = [1.,1.,1.]
  p = [determinant([[m[i][j] if j!=k else x[i] for j in range(3)] for i in range(3)]) for k in range(3)]+[-determinant(m)]
  f = dist(p[:3])
  return [-copysign(1.,p[3])*x/f for x in p]

def angle(a,b=None,o=org): 
  # angle AOB
  if b==None: b = [0,0,a[2]]
  p = plane(a,b,o)
  [u,v] = [proj2d(x,p,o) for x in [a,b]]
  return (atan2(u[1],u[0])-atan2(v[1],v[0]))

def ftheta(a,o=org):
  u = diff(a,o)
  return atan2((u[0]**2+u[1]**2)**0.5,u[2])
 
def fphi(a,o=org):
  u = diff(a,o)
  return atan2(u[1],u[0])

def clst_p2l(a,b,o):
  # closest point from the straight line AB to O
  return plus(mult(scalar(diff(b,a),diff(o,a))/dist(diff(b,a))**2,diff(b,a)),a)

def clst_l2l(a1,b1,a2,b2):
  # closest points C1,C2 from the straight lines A1B1 and A2B2
  a1b1,a2b2,a1a2 = diff(b1,a1),diff(b2,a2),diff(a2,a1)
  [k1,k2] = solve_matrix([[dist(a1b1)**2,-scalar(a1b1,a2b2)],[scalar(a1b1,a2b2),-dist(a2b2)**2]],[scalar(a1b1,a1a2),scalar(a2b2,a1a2)])
  if k1==None: return [None,None]
  return [plus(mult(k1,a1b1),a1),plus(mult(k2,a2b2),a2)]
  

def solve_matrix(a,b,w=False):
  d = determinant(a)
  n = len(a)
  if d==0:
    if w: sys.stderr.write('matrix not inversible\n')
    return [None for x in a]
  return [determinant([[a[i][j] if j!=k else b[i] for j in range(n)] for i in range(n)])/d for k in range(n)]
  
def determinant(m):
  if len(m)==1: return m[0][0]
  n = len(m)
  mk = [[[m[i][j] for j in range(n) if j!=k] for i in range(n) if i!=0] for k in range(n)]
  return sum([(-1)**k*m[0][k]*determinant(mk[k]) for k in range(n)])
