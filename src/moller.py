from misc import *
from ROOT import TLorentzVector
# from scipy.integrate import quad
# from scipy.interpolate import interp1d

m_e = 0.51099893e-3
degrad = pi/180.
gb = 389.379338e-6
alpha = 1/137.036

# lvac = readfile(os.environ['HOME']+'/python/vpol.dat',[float,float],start=1,tr=True)
# fvac = interp1d(lvac[0],lvac[1])

def moller_energy(theta,Ei=1.1,a=None,ml=m_e):
  if a==None: a = (Ei-ml)/(Ei+ml)
  return ml*(1 + a*cos(theta)**2)/(1 - a*cos(theta)**2)

def moller_center(x1,y1,x2,y2):
  # ax+by+c = 0
  a = y2-y1
  b = x1-x2
  c = -a*x1-b*y1
  return [-a*c/(a**2+b**2), -b*c/(a**2+b**2)]

def complementary_angle(theta,Ei=1.1,ml=m_e):
  Em = moller_energy(theta,Ei,ml=ml)
  return asin((Em**2-ml**2)**0.5/((Ei-Em)**2-ml**2)**0.5*sin(theta))

def moller_q2(theta,Ei=1.1,a=None,ml=m_e,Em=None):
  if Em==None: Em = moller_energy(theta,Ei,a,ml)
  return 2*Ei*Em - 2*ml**2 - 2*((Ei**2-ml**2)*(Em**2-ml**2))**0.5*cos(theta)

def moller_dq2_domega(theta,Ei=1.1,a=None,ml=m_e,Em=None):
  if a==None: a = (Ei-ml)/(Ei+ml)
  return -Ei*ml*(a**2*cos(theta)**4 - 4.*a*cos(theta)**2 + 4.*a*cos(theta) - 1.)/(pi*(a*cos(theta)**2 - 1)**2)


class MollerEvent:
  def __init__(self,theta1=1.*degrad,Ei=1.1,phi1=0.,ml=m_e,Eg_cut=1e-3):
    self.ml = ml
    self.Eg_cut = Eg_cut
    self.vli1 = TLorentzVector(0.,0.,(Ei**2-ml**2)**0.5,Ei)
    self.vli2 = TLorentzVector(0.,0.,0.,ml)
    Ef1 = moller_energy(theta1,Ei,None,ml)
    Pf1 = (Ef1**2-ml**2)**0.5
    self.vlf1 = TLorentzVector(Pf1*sin(theta1)*cos(phi1),Pf1*sin(theta1)*sin(phi1),Pf1*cos(theta1),Ef1)
    self.vlf2 = self.vli1+self.vli2-self.vlf1
    self._Ei = Ei
    self._Ef1 = Ef1
    self._Ef2 = self.vlf2.E()

  def __repr__(self):
    return '(MollerEvent: Ei={0}, theta1={1}, E1={2}, theta2={3}, E2={4})'.format(self.Ei,round(self.vlf1.Theta()/degrad,3),round(self.Ef1,3),round(self.vlf2.Theta()/degrad,3),round(self.Ef2,3))
    
  def _get_Ei(self): return self.vli1.E()
  def _get_Ef1(self): return self.vlf1.E()
  def _get_Ef2(self): return self.vlf2.E()
  def Q2(self): return -(self.vli1-self.vlf1).M2()
  
  Ei  = property(_get_Ei)
  Ef1  = property(_get_Ef1)
  Ef2  = property(_get_Ef2)

  def dsigma_domega(self,corr='',unit='GeV2'):
    t = (self.vlf1 - self.vli1).M2()
    u = (self.vlf1 - self.vli2).M2()
    s  = (self.vli1 + self.vli2).M2()
    Adir = (s - 2*self.ml**2)**2 + (u - 2*self.ml**2)**2 + 4*self.ml**2*t
    Aex = (s - 2*self.ml**2)**2 + (t - 2*self.ml**2)**2 + 4*self.ml**2*u
    Aint = -(s - 2*self.ml**2)*(s - 6*self.ml**2)
    if t==0 or u==0: return 0.
    ds = 2*alpha**2*(cos(self.vlf1.Theta())/(self.Ei + self.ml - (self.Ei - self.ml)*cos(self.vlf1.Theta())**2)**2)*(Adir/t**2 + Aex/u**2 - 2*Aint/(t*u))
    if corr=='brem': ds*=(1+self.d_brem_mol()+self.d_vac()+self.d_vertex())
    if unit=='barn': ds*= gb 
    return ds

  def dsigma_dtheta(self,corr='',unit='GeV2'): return self.dsigma_domega(corr,unit)*sin(self.vlf1.Theta())*2*pi
    
  def d_vac(self): return fvac(self.Q2())

  def d_vertex(self): return 2*(alpha/pi)*(3*log(self.Q2()/self.ml**2)/2 -2)

  def B(self,x,v1,v2):
    p_x = v1*x + v2*(1-x)
    return (v1*v2)*(log(4*self.Eg_cut**2/(p_x*p_x)) + p_x.E()*log((p_x.E() - p_x.Vect().Mag())/(p_x.E() + p_x.Vect().Mag()))/(p_x.Vect().Mag()))/(p_x*p_x)/(4*pi)
  
  def d_brem_mol(self):
    d_l1i_l1i = 0.5*(log(2.*self.Eg_cut/self.ml) + self.Ei*log(self.ml/(self.Ei + (self.Ei**2 - self.ml**2)**0.5))/(self.Ei**2 - self.ml**2)**0.5)/pi
    d_l1f_l1f = 0.5*(log(2.*self.Eg_cut/self.ml) + self.Ef1*log(self.ml/(self.Ef1 + (self.Ef1**2 - self.ml**2)**0.5))/(self.Ef1**2 - self.ml**2)**0.5)/pi
    d_l1i_l1f = quad(lambda x: self.B(x,self.vli1,self.vlf1),0.,1.)[0]
    d_l1i_l2i = quad(lambda x: self.B(x,self.vli1,self.vli2),0.,1.)[0]
    d_l1i_l2f = quad(lambda x: self.B(x,self.vli1,self.vlf2),0.,1.)[0]
    d_l1f_l2i = quad(lambda x: self.B(x,self.vlf1,self.vli2),0.,1.)[0]
    d_l1f_l2f = quad(lambda x: self.B(x,self.vlf1,self.vlf2),0.,1.)[0]
    d_l2i_l2i = (log(2*self.Eg_cut/self.ml)-1)/2/pi
    d_l2i_l2f = quad(lambda x: self.B(x,self.vli2,self.vlf2),0.,1.)[0]
    d_l2f_l2f =  0.5*(log(2*self.Eg_cut/self.ml) + self.Ef2*log(self.ml/(self.Ef2 + (self.Ef2**2 - self.ml**2)**0.5))/(self.Ef2**2 - self.ml**2)**0.5)/pi;
    return -2*alpha*(d_l1i_l1i - 2*d_l1i_l1f + d_l1f_l1f + d_l2i_l2i - 2*d_l2i_l2f + d_l2f_l2f + 2*d_l1i_l2i - 2*d_l1i_l2f - 2*d_l1f_l2i + 2*d_l1f_l2f)
