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
         # 'c': numpy.array([0.06, 0.05,0.01,0.01,-0.01]),\
         # roughly the peak of one-color
         'c1': 0.00,\
         'c2': 0.00,\
         'c3': 0.00,\
         'c4': 0.00,\
         'c5': 0.00,\
         # roughly the peak of one-color
         # 'alpha1': 0.0031+0.0001,\
         # 'alpha2': 0.0005+0.0001,\
         # 'alpha3': 0.0006+0.0001,\
         # 'alpha4': 0.0007+0.0001,\
         # 'alpha5': 0.0021+0.0001,\
         'alpha1': 0.000,\
         'alpha2': 0.000,\
         'alpha3': 0.000,\
         'alpha4': 0.000,\
         'alpha5': 0.000,\
         # roughly the peak of one-color
         'beta1': 0.0,\
         'beta2': 0.0,\
         'beta3': 0.0,\
         'beta4': 0.0,\
         'beta5': 0.0,\
         # roughly the peak of one-color
         # 'gamma01': 4.9882,\
         # 'gamma02': 3.0604,\
         # 'gamma03': 2.387,\
         # 'gamma04': 1.7696,\
         'gamma01': 4.9882,\
         'gamma03': 3.0604,\
         'gamma04': 2.387,\
         'gamma05': 1.7696,\
         # from init the best values for gamma 1
         'gamma11': 2.4,\
         'gamma12': 0.7,\
         'gamma13': 0.3,\
         'gamma14': -0.2,\
         'gamma15': -0.2,\
         # 'gamma1_': numpy.zeros(4)+2., \
         'prob0': 0.3,\
         'mag_int': mag_renorm, \
         'L_sigma': numpy.zeros(5)+0.05, \
         # 'L_Omega': numpy.identity(5), \
         'Delta_unit':Delta_simplex, 'Delta_scale': nsne/8.,\
         'k_unit': Delta_simplex, 'k_scale': nsne/8., 'k_zero': 1./nsne} \
        for _ in range(4)]

sm = pystan.StanModel(file='gerard9.stan')
control = {'stepsize':1.}
fit = sm.sampling(data=data, iter=1000, chains=4,control=control,init=init,thin=1)
print fit

output = open('temp9.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
