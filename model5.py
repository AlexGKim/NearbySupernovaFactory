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

Delta_simplex = numpy.zeros(nsne-1)
# Delta_simplex = numpy.zeros(nsne)+1./nsne
# k_simplex = numpy.zeros(nsne)
R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne

init = [{'EW' : EW_renorm, \
         'c_raw' : numpy.zeros(5), \
         'alpha_raw' : numpy.zeros(5), \
         'beta_raw' : numpy.zeros(5), \
         'L_sigma_raw': numpy.zeros(5)+0.03*100, \
         # roughly the peak of one-color
         # 'gamma01': 5.2*1.,\
         # 'gamma03': 3.4*1.,\
         # 'gamma03': 3.4*1.,\
         # 'gamma04': 3.0*1.,\
         # 'gamma05': 2.4*1.,\
         'gamma01': 6./5,\
         'gamma02': 4./5,\
         'gamma03': 3./5,\
         'gamma04': 2./5,\
         'gamma05': 1./5,\
         'mag_int_raw': mag_renorm, \
         'L_Omega': numpy.identity(5), \
         'Delta_unit':R_simplex, \
         'Delta_scale': 2.5, \
         'k_unit': R_simplex, \
         # 'k_scale': 20, \
         'R_unit': R_simplex, \
         # 'R_scale': 20, \
         # 'rho11': 4.4/100.,\
         # 'rho13': 2.6/100.,\
         # 'rho14': 2.2/100.,\
         # 'rho15': 1.8/100.,\
         'rho11': 17./5,\
         'rho12': 0.*3./5,\
         'rho13': 0./5,\
         'rho14': 0.*3./5,\
         'rho15': 0.*3./5,\
         # 'rho1':numpy.zeros(5)\
         } \
        for _ in range(4)]

sm = pystan.StanModel(file='gerard5.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=2000, chains=4,control=control,init=init, thin=1)

output = open('temp5.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit
