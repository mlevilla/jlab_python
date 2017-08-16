from misc_root import *
# from scipy.special import spence
# from scipy.integrate import quad
# from scipy.interpolate import interp1d

m_e = 0.51099893e-3
m_p = 938.272046e-3
mu_p = 2.79284736
alpha = 1/137.036
gb = 389.379338e-6
degrad = pi/180

# lvac = readfile(os.environ['HOME']+'/python/vpol.dat',[float,float],start=1,tr=True)
# fvac = interp1d(lvac[0],lvac[1])

def GE(q2,flag=2,a11=0.,b11=0.,b12=0.,b13=0.,mn=m_p):
  tau = abs(q2)/(4.*mn**2)
  if flag==1: return 1. # point-like proton
  elif flag==2: return 1./(1.+abs(q2)/0.71)**2
  elif flag==3: return (1. - 0.24*tau)/(1. + 10.98*tau + 12.82*tau**2 + 21.97*tau**3); # Kelly parametrization
  elif flag==4: return (1. - 0.299*tau)/(1. + 11.11*tau + 14.11*tau**2 + 15.7*tau**3) # Puckett parametrization
  elif flag==5: return (1. + 2.90966*tau - 1.11542229*tau**2 + 3.866171e-2*tau**3)/(1. + 14.5187212*tau + 40.88333*tau**2 + 99.999998*tau**3 + 4.579e-5*tau**4 + 10.3580447*tau**5)
  elif flag==6: return (1. + a11*tau)/(1. + b11*tau + b12*tau**2 + b13*tau**3) 

def GM(q2,flag=2,a21=0.,b21=0.,b22=0.,b23=0.,mn=m_p,mu=mu_p):
  tau = abs(q2)/(4.*mn**2)
  if flag==1: return mu
  elif flag==2: return mu*GE(q2,mn=mn)
  elif flag==3: return mu*(1. + 0.12*tau)/(1. + 10.98*tau + 18.86*tau**2 + 6.55*tau**3) #Kelly parametrization
  elif flag==4: return mu*(1. + 0.081*tau)/(1. + 11.15*tau + 18.45*tau**2 + 5.31*tau**3) #Puckett parametrization
  elif flag==5: return mu*(1. - 1.43573*tau + 1.19052066*tau**2 + 2.5455841e-1*tau**3)/(1. + 9.70703681*tau + 3.7357e-4*tau**2 + 6.0e-8*tau**3) + 9.9527277*tau**4 + 12.7977739*tau**5
  elif flag==6: return mu*(1. + a21*tau)/(1. + b21*tau + b22*tau**2 + b23*tau**3)
  
def F1(q2,flag=2,a11=0.,b11=0.,b12=0.,b13=0.,a21=0.,b21=0.,b22=0.,b23=0.,mn=m_p,mu=mu_p):
  if flag==1: return 1. # point-like proton
  tau = abs(q2)/(4.*mn**2)
  return (GE(q2,flag,a11,b11,b12,b13,mn)+tau*GM(q2,flag,a21,b21,b22,b23,mn))/(1+tau)

def F2(q2,flag,a11=0.,b11=0.,b12=0.,b13=0.,a21=0.,b21=0.,b22=0.,b23=0.,mn=m_p,mu=mu_p):
  if flag==1: return mu-1
  tau = abs(q2)/(4.*mn**2)
  return (GM(q2,flag,a21,b21,b22,b23,mn,mu)-GE(q2,flag,a11,b11,b12,b13,mn))/(1+tau)

def LVector(l): return TLorentzVector(l[0],l[1],l[2],l[3])

def ep_energy_el(thetaf,Ei=1.1,ml=m_e,mn=m_p): return ((Ei+mn)*(mn*Ei+ml**2) + (mn**2-(ml*sin(thetaf))**2)**0.5*(Ei**2-ml**2)*cos(thetaf))/((Ei+mn)**2 - (Ei**2-ml**2)*cos(thetaf)**2)

def ep_energy_br(thetaf,Eg,thetag,phig,en_sign,Ei=1.1,ml=m_e,mn=m_p):
  A = (Ei**2 - ml**2)**0.5*cos(thetaf) - Eg*(cos(thetaf)*cos(thetag) + sin(thetaf)*sin(thetag)*cos(phig))
  B = Ei + mn - Eg
  C = Eg*(Ei + mn - (Ei**2 - ml**2)**0.5*cos(thetag)) - mn*Ei - ml**2    
  Ef = (B*C + en_sign*A*(ml**2*(A**2 - B**2) + C**2)**0.5)/(A**2 - B**2)

  if not ((ml**2*(A**2 - B**2) + C**2 < 0.) or (abs(A**2 - B**2) < 1.e-12) or (abs(A*(Ef**2 - ml**2)**0.5 - B*Ef - C) > 1.e-9) or (Ef < ml or Ef > Ei - Eg)): return Ef

def ep_q2(theta,Ei=1.1,ml=m_e,mn=m_p,Eep=None):
  if Eep==None: Eep = ep_energy_el(theta,Ei,ml,mn)
  return 2*Ei*Eep - 2*ml**2 - 2*((Ei**2-ml**2)*(Eep**2-ml**2))**0.5*cos(theta)

def ep_dq2_domega(theta,Ei=1.1,ml=m_e,mn=m_p):
  return 2*mn*(Ei**2 - ml**2)*(mn**2 - ml**2*sin(theta)**2)**(-0.5)*((mn**2 - ml**2*sin(theta)**2)**0.5*((Ei + mn)*(Ei*mn + ml**2) + (Ei**2 - ml**2)*(mn**2 - ml**2*sin(theta)**2)**0.5*cos(theta))*cos(theta) + 0.5*(ml**2*cos(theta)**2 + (mn**2 - ml**2*sin(theta)**2))*((Ei + mn)**2 + (-Ei**2 + ml**2)*cos(theta)**2))/(pi*((Ei + mn)**2 + (-Ei**2 + ml**2)*cos(theta)**2)**2)

def mott(theta=None,Ei=1.1,Ef=None,Q2=None,unit='barn',ml=m_e,mn=m_p):
  if Ef is None: Ef = ep_energy_el(theta,Ei,ml,mn)
  pi = (Ei**2-ml**2)**0.5
  pf = (Ef**2-ml**2)**0.5
  if Q2 is None: Q2 = pf**2*sin(theta)**2+(pi-pf*cos(theta))**2-(Ei-Ef)**2
  if Q2==0: return 0.
  d = (Ef*pi)/(Ei*pf)
  ds = (alpha/(2*Ei))**2*(1+Q2/(4*Ei*Ef))/(Q2/(4*Ei*Ef))**2/d*mn*(Ef**2-ml**2)/(mn*Ei*Ef+ml**2*(Ef-Ei-mn))
  if unit=='barn': ds *= gb
  return ds

def born(theta=None,Ei=1.1,Ef=None,Q2=None,unit='barn',ml=m_e,mn=m_p,flag=2):
  if Ef is None: Ef = ep_energy_el(theta,Ei,ml,mn)
  if Q2 is None: 
    pi = (Ei**2-ml**2)**0.5
    pf = (Ef**2-ml**2)**0.5
    Q2 = pf**2*sin(theta)**2+(pi-pf*cos(theta))**2-(Ei-Ef)**2
  tau = -Q2/(4*mn**2)
  eps = 1/(1 - 2*(1 + tau)*(Q2 + 2*ml**2)/(4*Ei*Ef + Q2))
  return mott(Ei=Ei,Ef=Ef,Q2=Q2,ml=ml,mn=mn,unit=unit)*(eps*GE(Q2,flag)**2+tau*GM(Q2,flag)**2)/(eps*(1+tau))

class EpEvent:
  def __init__(self,theta=1.*degrad,Ei=1.1,phi=0.,ml=m_e,mn=m_p,Eg_cut=1e-3):
    self.ml = ml
    self.mn = mn
    self.Eg_cut = Eg_cut
    self.vli = TLorentzVector(0.,0.,(Ei**2-ml**2)**0.5,Ei)
    self.vpi = TLorentzVector(0.,0.,0.,mn)
    Ef = ep_energy_el(theta,Ei,ml,mn)
    Pf = (Ef**2-ml**2)**0.5
    self.vlf = TLorentzVector(Pf*sin(theta)*cos(phi),Pf*sin(theta)*sin(phi),Pf*cos(theta),Ef)
    self.vpf = self.vli+self.vpi-self.vlf
    self._Ei = Ei
    self._Ef = Ef
    self._Ep = self.vpf.E()
    self.theta = theta

  def _get_Ei(self): return self.vli.E()
  def _get_Ef(self): return self.vlf.E()
  def _get_Ep(self): return self.vpf.E()
  def Q2(self): return -(self.vli-self.vlf).M2()
  def tau(self): return -self.Q2()/(4*self.mn**2)
  def epsilon(self): return 1/(1 - 2*(1 + self.tau())*(self.Q2() + 2*self.ml**2)/(4*self.Ei*self.Ef + self.Q2()))

  Ei  = property(_get_Ei)
  Ef  = property(_get_Ef)
  Ep  = property(_get_Ep)

  def mott(self,corr='',unit='barn'):
    ds = mott(theta=self.theta,Ei=self.Ei,Ef=self.Ef,mn=self.mn,ml=self.ml,unit=unit)
    if corr=='brem': ds*=(1+self.d_brem_ee()+self.d_brem_pp()+self.d_brem_ep()+self.d_vac()+self.d_vertex()+self.d_prime())
    return ds

  def born(self,corr='',flag=2,unit='barn'):
    tau = self.tau()
    eps = self.epsilon()
    Q2 = self.Q2()
    return self.mott(corr,unit)*(eps*GE(Q2,flag)**2+tau*GM(Q2,flag)**2)/(eps*(1+tau))

  def dsigma_domega(self,corr='',unit='GeV2',flag=2):
    return self.born(corr,flag,unit)

  def dsigma_dtheta(self,corr='',unit='GeV2'): 
    return self.dsigma_domega(corr,unit)*sin(self.vlf.Theta())*2*pi

  def GE(self):
    return ((1+self.tau())/self.mott())**0.5

  def GE_brem(self):
    return ((1+self.tau())/self.mott_brem())**0.5

  def GE_wGM(self,flag=2):
    return (self.GE()**2-self.tau/self.eps()*GM(self.Q2(),flag)**2)**0.5

  def GE_brem_wGM(self,flag=2):
    return (self.GE_brem()**2-self.tau/self.eps()*GM(self.Q2(),flag)**2)**0.5

  def B(self,x,v1,v2):
    p_x = v1*x + v2*(1-x)
    return (v1*v2)*(log(4*self.Eg_cut**2/(p_x*p_x)) + p_x.E()*log((p_x.E() - p_x.Vect().Mag())/(p_x.E() + p_x.Vect().Mag()))/(p_x.Vect().Mag()))/(p_x*p_x)/(4*pi)


  def d_brem_ee(self):
    d_li_li = 0.5*(log(2*self.Eg_cut/self.ml) + self.Ei*log(self.ml/(self.Ei + (self.Ei**2 - self.ml**2)**0.5))/(self.Ei - self.ml**2)**0.5)/pi
    d_li_lf = quad(lambda x: self.B(x,self.vli,self.vlf),0.,1.)[0]
    d_lf_lf = 0.5*(log(2*self.Eg_cut/self.ml) + self.Ef*log(self.ml/(self.Ef + (self.Ef**2 - self.ml**2)**0.5))/(self.Ef**2 - self.ml**2)**0.5)/pi
    
    return -2*alpha*(d_li_li - 2*d_li_lf + d_lf_lf)

  def d_brem_pp(self):
    d_pi_pi = (log(2*self.Eg_cut/self.mn) - 1)/(2*pi)
    d_pi_pf = quad(lambda x: self.B(x,self.vpi,self.vpf),0., 1.)[0]
    d_pf_pf = 0.5*(log(2*self.Eg_cut/self.mn) + self.Ep*log(self.mn/(self.Ep + (self.Ep**2 - self.mn**2)**0.5))/(self.Ep**2 - self.mn**2)**0.5)/pi    
    return -2*alpha*(d_pi_pi - 2*d_pi_pf + d_pf_pf)

  def d_brem_ep(self):
    d_li_pi = quad(lambda x: self.B(x,self.vli,self.vpi),0., 1.)[0]
    d_li_pf = quad(lambda x: self.B(x,self.vli,self.vpf),0., 1.)[0]
    d_lf_pi = quad(lambda x: self.B(x,self.vlf,self.vpi),0., 1.)[0]
    d_lf_pf = quad(lambda x: self.B(x,self.vlf,self.vpf),0., 1.)[0]
    return 4*alpha*(d_li_pi - d_li_pf - d_lf_pi + d_lf_pf)

  def d_vac(self): return fvac(self.Q2())[0]

  def d_vertex(self): return (alpha/pi)*(3*log(self.Q2()/self.ml**2)/2 -2)

  def d_prime(self):
    return -(alpha/pi)*(log(self.Ei/self.Ef)*log(self.Q2()**2/(4*self.mn**2*self.Ei*self.Ef)) + 2.*spence(1-0.5*self.mn/self.Ei) - 2*spence(1-0.5*self.mn/self.Ef))
