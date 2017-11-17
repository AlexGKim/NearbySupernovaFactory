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

sivel, sivel_err, _, _, _, _, _ = sivel.sivel(data)

use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn
data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_renorm, 'sivel_err': sivel_err}

# pystan.misc.stan_rdump(data, 'data.R')
# wefew

Delta_simplex = numpy.zeros(nsne-1)
# Delta_simplex = numpy.zeros(nsne)+1./nsne
# k_simplex = numpy.zeros(nsne)
R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()

init = [{'EW' : EW_renorm, \
         'sivel': sivel_renorm,\
         'c_raw' : numpy.zeros(5), \
         'alpha_raw' : numpy.zeros(5), \
         'beta_raw' : numpy.zeros(5), \
         'eta_raw' : numpy.zeros(5), \
         'ev_sig': 0.1, \
         'ev': numpy.array([1.,1,1,1,1])/numpy.sqrt(5),\
         # 'L_sigma_raw': numpy.zeros(5)+0.03*100, \
         'gamma01': 6./5,\
         'gamma02': 4./5,\
         'gamma03': 3./5,\
         'gamma04': 2./5,\
         'gamma05': 1./5,\
         # 'mag_int_raw': numpy.mean(mag_renorm,axis=1), \
         # 'L_Omega': numpy.identity(5), \
         'Delta_unit':R_simplex, \
         'Delta_scale': 15./4, \
         'k_unit': R_simplex, \
         'R_unit': R_simplex, \
         'rho11': -17./5,\
         'rho12': 0.*3./5,\
         'rho13': 0./5,\
         'rho14': 0.*3./5,\
         'rho15': 0.*3./5,\
         } \
        for _ in range(8)]


sm = pystan.StanModel(file='fix1.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=5000, chains=8,control=control,init=init, thin=1)


output = open('fix1.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit