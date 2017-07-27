from math import log,cos

def energy_loss(theta=0.,Ei=1100.):
  ZoverA = [13./26.98, 10.6/21.8 , 0.49919]
  density = [2.7, 0.1117433, 1.205e-3]
  I = [166.e-6, 106.6e-6, 85.7e-6]
  length = [0.2, 1.5, 50];
  de = 5.0989 * 1.e-25
  Na = 6.02214086e23
  me = 0.5109989181
  gamma = Ei/me
  eDep = 0.
  for i in range(3):
    dedx = 0.5 * de * Na * density[i] * ZoverA[i] * (2 * log(2*me/I[i]) + 3*log(gamma) - 1.95)
    eDep+=(dedx*length[i]/cos(theta))
  return eDep
