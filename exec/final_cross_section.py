#!/usr/bin/env python
from params import *

largs = [('files',['',''],'','files'),
         ('n','','','suffix'),
         ('o',workf,'','output folder'),
         ('ebeam',2.142,'','beam_energy'),
         ('ff',5,'','form_factor'),
         ('int',[0.8,2.0],'','integration range')]

[files,suffix,outdir,ebeam,ff,intrange] = ArgParser(largs).argl

from misc_root import *
from esepp import *


fout = TFile(outdir+'/final_cross_section{0}.root'.format(suffix),'recreate')
fdata = TFile(files[0])
fsim = TFile(files[1])

thbin = get_binning(fdata.theta.Get('hth_eff_ep'))
nbin = len(thbin)-1
max_range = [i for i in range(nbin) if thbin[i]<=intrange[1]<thbin[i+1]][0]

gthdata = [[fdata.theta.Get('grth'+['','2'][j]+'exp_'+['','eff_'][j]+'ee'+str(i)) for i in range(1,4)] for j in range(2)]
gthsim = [[fsim.theta.Get('grth'+['','2'][j]+'sim_ee'+str(i)) for i in range(1,4)] for j in range(2)]
gq2data = [[fdata.q2.Get('grq2'+['','2'][j]+'exp_'+['','eff_'][j]+'ee'+str(i)) for i in range(1,4)] for j in range(2)]
gq2sim = [[fsim.q2.Get('grq2'+['','2'][j]+'sim_ee'+str(i)) for i in range(1,4)] for j in range(2)]

thbin_data = [x for x in gthdata[1][0].GetX()]
q2bin_data = [x for x in gq2data[1][0].GetX()]
rbin_data = [x/0.197327**2 for x in gq2data[1][0].GetX()]

# combination of the two cross_section (bin by bin and intgrated)
ydatathcomb =  [gthdata[0][2].GetY()[i] for i in range(max_range)]+[gthdata[1][2].GetY()[i] for i in range(max_range,gthdata[0][2].GetN())]
ydatathcomb_err =  [gthdata[0][2].GetEY()[i] for i in range(max_range)]+[gthdata[1][2].GetEY()[i] for i in range(max_range,gthdata[0][2].GetN())]
ysimthcomb =  [gthsim[0][2].GetY()[i] for i in range(max_range)]+[gthsim[1][2].GetY()[i] for i in range(max_range,gthsim[0][2].GetN())]
ysimthcomb_err =  [gthsim[0][2].GetEY()[i] for i in range(max_range)]+[gthsim[1][2].GetEY()[i] for i in range(max_range,gthsim[0][2].GetN())]
gthdata_comb = tgraph(thbin_data,ydatathcomb,dy=ydatathcomb_err,name='gthdata_comb',title=';#theta (deg);d#sigma_{ep}/d#sigma_{ee}')
gthsim_comb = tgraph(thbin_data,ysimthcomb,dy=ysimthcomb_err,name='gthsim_comb',title=';#theta (deg);d#sigma_{ep}/d#sigma_{ee}')

# theoretical calculation
ev_ep = [[EpEvent((thbin[i]+(thbin[i+1]-thbin[i])/500*j)*degrad,ebeam) for j in range(501)] for i in range(nbin)] 
cs_ep0 = [[ev_ep[i][j].dsigma_domega(unit='barn',flag=ff) for j in range(501)] for i in range(nbin)]
cs_ep = [sum(cs_ep0[i][j]*sin((thbin[i]+(thbin[i+1]-thbin[i])/500*j)*degrad)*(thbin[i+1]-thbin[i])/500*degrad for j in range(501))/(cos(thbin[i]*degrad)-cos(thbin[i+1]*degrad)) for i in range(nbin)]
thbin_theo = [sum((thbin[i]+(thbin[i+1]-thbin[i])/500*j)*cs_ep0[i][j] for j in range(501))/cs_ep[i]/501 for i in range(nbin)]
q2bin_theo = [sum(ev_ep[i][j].Q2()*cs_ep0[i][j] for j in range(501))/cs_ep[i]/501 for i in range(nbin)]

# cross-section ratio
lsigma = [[root_op([x1,x2,cs_ep],[],['/','*']) for x1,x2 in zip(y1,y2)] for y1,y2 in zip(gthdata,gthsim)]

# strings for titles
thn = ';#theta (deg);'
q2n = ';Q^{2} (GeV^{2});'
xsn = '#frac{d#sigma}{d#Omega} (b/sr)'

gth_sigma = [[tgraph(thbin_data,lsigma[j][i][0],dy=lsigma[j][i][1],name='gth'+str(j)+'_ee'+str(i+1),title=thn+xsn) for i in range(3)] for j in range(2)]
gq2_sigma = [[tgraph(q2bin_data,lsigma[j][i][0],dy=lsigma[j][i][1],name='gth'+str(j)+'_ee'+str(i+1),title=q2n+xsn) for i in range(3)] for j in range(2)]

gth_comb = tgraph(thbin_data,lsigma[0][2][0][:max_range]+lsigma[1][2][0][max_range:],dy=lsigma[0][2][1][:max_range]+lsigma[1][2][1][max_range:],name='gth_comb_ee3',title=thn+xsn)
gq2_comb = tgraph(q2bin_data,lsigma[0][2][0][:max_range]+lsigma[1][2][0][max_range:],dy=lsigma[0][2][1][:max_range]+lsigma[1][2][1][max_range:],name='gq2_comb_ee3',title=q2n+xsn)


# GE and other quantities
ev_data = [EpEvent(th*degrad,ebeam) for th in thbin_data]
gsigma_omott = tgraph(q2bin_data,gq2_comb,op=['/',[ev.mott() for ev in ev_data]],name='gq2_over_mott',title=q2n+'#frac{d#sigma}{d#sigma_{Mott}}')

ge = [(y/ev.mott()/(1+ev.tau())-ev.tau()/ev.epsilon()*GM(ev.Q2(),flag=ff))**0.5 for y,ev in zip(gq2_comb.GetY(),ev_data)]
dge = [x*dy/y if y!=0 else 0. for x,y,dy in zip(ge,gq2_comb.GetY(),gq2_comb.GetEY())]

gq2_ge = tgraph(q2bin_data,ge,dy=dge,name='gq2_ge_ee3',title=q2n+'G_{E}')
gr_ge = tgraph(rbin_data,ge,dy=dge,name='gr_ge_ee3',title=';r^{-2} (fm^{-2});G_{E}')

ffit_q2 = TF1('ffit_q2','[0]-[1]^2*x/6/0.197327^2',0,0.05)
ffit_r = TF1('ffit_r','[0]-[1]^2*x/6',0,1)

gq2_ge.Fit(ffit_q2)
gq2_ge.Fit(ffit_q2)

gth_ep = tgraph(thbin_theo,cs_ep,name='gthep_theo',title=thn+xsn)
gq2_ep = tgraph(q2bin_theo,cs_ep,name='gq2ep_theo',title=q2n+xsn)

# comparison data/sim
gthcomp = [tgraph(thbin_data,gthdata[i][2],op=['/',gthsim[i][2]],name='gth'+['','2'][i]+'comp',title=thn) for i in range(2)]
gthcomp_comb = tgraph(thbin_data,gthdata_comb,op=['/',gthsim_comb],name='gthcomb_comp',title=thn)

writelisto([gthdata,gthsim,gthdata_comb,gthsim_comb,gth_sigma,gth_comb,gth_ep,gthcomp,gthcomp_comb],fout,['theta'],yrg=[0.,5.])
writelisto([gq2data,gq2sim,gq2_sigma,gq2_comb,gq2_ep,gsigma_omott,gq2_ge,gr_ge],fout,['q2'],yrg=[0.,5.])

fout.Close()

