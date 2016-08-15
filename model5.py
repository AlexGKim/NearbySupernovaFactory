#!/usr/bin/env python

import pickle
import numpy
import pystan

# two color parameter model

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov}


Delta_simplex = numpy.zeros(nsne)+1./nsne
k_simplex = numpy.zeros(nsne)
# R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne


init = [{'EW' : EW_renorm, \
         'c1': -3.2,\
         'c2': -2.8,\
         'c3': -1.8,\
         'c4': -1.5,\
         'c5': -1.1,\
         # roughly the peak of one-color
         'alpha1': 0.003,\
         'alpha2': 0.001,\
         'alpha3': 0.0007,\
         'alpha4': 0.0008,\
         'alpha5': 0.0022,\
         # roughly the peak of one-color
         'beta1': 0.034,\
         'beta2': 0.028,\
         'beta3': 0.027,\
         'beta4': 0.023,\
         'beta5': 0.022,\
         # roughly the peak of one-color
         'gamma01': 5.,\
         'gamma02': 3.1,\
         'gamma03': 2.4,\
         'gamma04': 1.8,\

         # 'alpha': numpy.zeros(5), \
         # 'beta':numpy.zeros(5), \
         # 'gamma0': 0.1, 'gamma_': numpy.zeros(4), \
         'mag_int': mag_renorm, \
         'L_sigma': numpy.zeros(5)+0.03, \
         # 'L_Omega': numpy.identity(5), \
         'Delta_unit':Delta_simplex, 'Delta_scale': nsne/8.,\
         'k_unit': k_simplex, \
         'R':k_simplex+0.1, \
         'rho11': 3.6,\
         'rho12': 1.8,\
         'rho13': 1.4,\
         'rho14': 0.9
         # 'rho1':numpy.zeros(5)\
         } \
        for _ in range(4)]

sm = pystan.StanModel(file='gerard5.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=8000, chains=4,control=control,init=init, thin=4)
print fit

output = open('temp5.pkl','wb')
pickle.dump(fit.extract(), output)
output.close()
