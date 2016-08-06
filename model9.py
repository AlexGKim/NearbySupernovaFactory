#!/usr/bin/env python

import pickle
import numpy
import pystan

# one color parameter model two dust

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

init = [{'EW' : EW_renorm, \
         'c': numpy.zeros(5)-.1,\
         # 'alpha': numpy.array([0.0031, 0.0005,0.0006,0.0007,0.0021]), \
         'alpha1': 0.001,\
         'alpha2': -0.0016,\
         'alpha3': -0.0015,\
         'alpha4': -0.0015,\
         'alpha5': 0.0021,\
         # 'beta':numpy.array([0.0345, 0.0274, 0.0274, 0.0223, 0.0213]), \
         'beta1': 0.0125,\
         'beta2': 0.0055,\
         'beta3': 0.006,\
         'beta4': 0.001,\
         'beta5': 0.0213,\
         # 'gamma_': numpy.array([5.0,3.1,2.4,1.8]), \
         'gamma01': 1.95,\
         'gamma02': 3.1,\
         'gamma03': -0.7,\
         'gamma04': -1.3,\
         'gamma1_': numpy.zeros(4)+2., \
         'prob0': 0.3,\
         'mag_int': mag_renorm, \
         'L_sigma': numpy.zeros(5)+0.05, \
         # 'L_Omega': numpy.identity(5), \
         'Delta_unit':Delta_simplex, 'Delta_scale': nsne/8.,\
         'k': k_simplex} \
        for _ in range(4)]

sm = pystan.StanModel(file='gerard9.stan')
control = {'stepsize':1.}
fit = sm.sampling(data=data, iter=2000, chains=4,control=control,init=init,thin=2)
print fit

output = open('temp9.pkl','wb')
pickle.dump(fit.extract(), output)
output.close()
