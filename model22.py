#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel

# two color parameter model

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel, sivel_err, x1, x1_err = sivel.sivel(data)

# dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

# dic_meta=cPickle.load(open("META.pkl"))

# sivel=[]
# sivel_err=[]


# for sn in data['snlist']:
#    sn = 'PTF12iiq'
#    if sn in dic_meta.keys() and sn in dic_phreno.keys():
#       meta = dic_meta[sn]
#       vSiII_6355_lbd=0.
#       vSiII_6355_lbd_err=0.
#       counter  = 0
#       for sp in dic_phreno[sn]["spectra"]:
#          if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase']) < 2.5 and numpy.isfinite(dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]):
#             vSiII_6355_lbd += dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
#             vSiII_6355_lbd_err += 1/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
#             print dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]
#             counter +=1
#       if counter !=0:
#          sivel.append(vSiII_6355_lbd / vSiII_6355_lbd_err)
#          sivel_err.append(1./numpy.sqrt(vSiII_6355_lbd_err)) 
#       else:
#          sivel.append(float('nan'))
#          sivel_err.append(float('nan'))
#    else:
#       sivel.append(float('nan'))
#       sivel_err.append(float('nan'))


# sivel = numpy.array(sivel)
# sivel_err = numpy.array(sivel_err)

use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
x1=x1[use]
x1_err = x1_err[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

snname = numpy.array(data['snlist'])[use]

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn
data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_renorm, 'sivel_err': sivel_err,'x1_obs': x1, 'x1_err':x1_err}

Delta_simplex = numpy.zeros(nsne-1)
# Delta_simplex = numpy.zeros(nsne)+1./nsne
# k_simplex = numpy.zeros(nsne)
R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()

init = [{'EW' : EW_renorm, \
         'sivel': sivel_renorm,\
         'x1': x1,\
         'c_raw' : numpy.zeros(5), \
         'alpha_raw' : numpy.zeros(5), \
         'beta_raw' : numpy.zeros(5), \
         'eta_raw' : numpy.zeros(5), \
         'zeta' : numpy.zeros(5), \
         'L_sigma_raw': numpy.zeros(5)+0.03*100, \
         'gamma01': 61./5,\
         'gamma02': 50./5,\
         'gamma03': 40./5,\
         'gamma04': 30./5,\
         'gamma05': 20./5,\
         'gamma11': -12./5,\
         'gamma12': -14./5,\
         'gamma13': -16./5,\
         'gamma14': -14./5,\
         'gamma15': -14/5,\
         'mag_int_raw': mag_renorm, \
         'L_Omega': numpy.identity(5), \
         'Delta_unit':R_simplex, \
         'Delta_scale': 15./4, \
         'k_unit': R_simplex, \
         'k1_unit': R_simplex, \
         'R_unit': R_simplex, \
         # 'rho11': 2./5,\
         # 'rho12': -2./5,\
         # 'rho13': -2./5,\
         'rho11': 0./5,\
         'rho12': 0./5,\
         'rho13': -0./5,\
         # 'rho14': 1./5,\
         # 'rho15':  1./5,\
         } \
        for _ in range(8)]

sm = pystan.StanModel(file='gerard22.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=5000, chains=8,control=control,init=init, thin=1)


output = open('temp22.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit