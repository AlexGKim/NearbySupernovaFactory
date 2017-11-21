#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel

f = open('fix3.pkl','rb')
(fit,_) = pickle.load(f)

gamma0 = numpy.median(fit['gamma'],axis=0)
gamma1 = numpy.median(fit['rho1'],axis=0)
gamma0_cov = numpy.cov(fit['gamma'],rowvar=False)
gamma1_cov = numpy.cov(fit['rho1'],rowvar=False)


_, gamma0_ev = numpy.linalg.eigh(gamma0_cov)

# print gamma0_ev[:,0]

ev_median = numpy.median(fit['ev']*numpy.sign(fit['ev'][:,4][:,None]),axis=0)
ev_median = ev_median/numpy.linalg.norm(ev_median)

evsig_median = numpy.median(fit['ev_sig'])


gamma0_eval = numpy.array(fit['gamma'])
for i in xrange(gamma0_eval.shape[0]):
   for j in xrange(5):
      gamma0_eval[i,j]=numpy.dot(fit['gamma'][i,:], gamma0_ev[:,j])

# gamma0_min = numpy.min(gamma0_eval,axis=0)-0.05
# gamma0_max = numpy.max(gamma0_eval,axis=0)+0.05

# for i in xrange(5):
#    plt.hist(gamma0_eval[:,i])
#    plt.show()

_, gamma1_ev = numpy.linalg.eigh(gamma1_cov)

gamma1_eval = numpy.array(fit['rho1'])
for i in xrange(gamma1_eval.shape[0]):
   for j in xrange(5):
      gamma1_eval[i,j]=numpy.dot(fit['rho1'][i,:], gamma1_ev[:,j])

gamma0median = numpy.median(gamma0_eval,axis=0)
gamma1median = numpy.median(gamma1_eval,axis=0)

gamma1_min = (numpy.min(gamma1_eval,axis=0)-gamma1median)*3+gamma1median
gamma1_max = (numpy.max(gamma1_eval,axis=0)-gamma1median)*3+gamma1median


# two color parameter model

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel, sivel_err, x1, x1_err, _, _, _ = sivel.sivel(data)

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

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn
data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_renorm, 'sivel_err': sivel_err,'x1_obs': x1, 'x1_err':x1_err, \
   'rho1_ev':gamma1_ev, 'rho1_min':gamma1_min, 'rho1_max': gamma1_max }

# pystan.misc.stan_rdump(data, 'data.R')
# wefew

Delta_simplex = numpy.zeros(nsne-1)
# Delta_simplex = numpy.zeros(nsne)+1./nsne
# k_simplex = numpy.zeros(nsne)
R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()

numpy.random.seed(100)
ruv = []
for _ in range(8):
   temp = numpy.random.uniform(-1,1,5)
   ruv.append(temp/numpy.linalg.norm(temp))

init = [{'EW' : EW_renorm, \
         'sivel': sivel_renorm,\
         'x1': x1,\
         'c_raw' : numpy.zeros(5), \
         'alpha_raw' : numpy.zeros(5), \
         'beta_raw' : numpy.zeros(5), \
         'eta_raw' : numpy.zeros(5), \
         'zeta' : numpy.zeros(5), \
         'ev_sig': evsig_median, \
         'ev': ev_median,\
         # 'L_sigma_raw': numpy.zeros(5)+0.03*100, \
         'gamma01': gamma0[0],\
         'gamma02': gamma0[1],\
         'gamma03': gamma0[2],\
         'gamma04': gamma0[3],\
         'gamma05': gamma0[4],\
         'mag_int_raw': numpy.mean(mag_renorm,axis=1), \
         # 'L_Omega': numpy.identity(5), \
         'Delta_unit':R_simplex, \
         'Delta_scale': 15./4, \
         'k_unit': R_simplex, \
         'R_unit': R_simplex, \
         'rho11': gamma1median[0],\
         'rho12': gamma1median[1],\
         'rho13': gamma1median[2],\
         'rho14': gamma1median[3],\
         'rho15': gamma1median[4]\
         } \
        for _ in range(8)]


sm = pystan.StanModel(file='fix3_x1.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=5000, chains=8,control=control,init=init, thin=1)


output = open('fix3_x1.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit