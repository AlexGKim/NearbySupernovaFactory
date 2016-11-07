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

sivel, sivel_err,_,_ = sivel.sivel(data)


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

# Truth

Delta = numpy.random.uniform(low=0.2, high=0.2, size=nsne)
c = numpy.zeros(5)
alpha = numpy.array([0.005, 0.002,0.002,0.002,0.003])
beta = numpy.array([0.03, 0.03,0.03,0.02,0.02])
eta = numpy.array([-0.0002,0,0.0005,0.0005,-0.0003])
delta = numpy.array([-0.737, 0.129,1.961,0.598,-2.545])
gamma = numpy.array([62.5,50.7,37.7,28.9,20.6])
gamma1 = numpy.array([-9.93,-11.8,-14.2,-13.4,-12.1])
sigma = numpy.array([0.06,0.03,0.02,0.01,0.04])

EW = numpy.array(EW_renorm)
sivel = numpy.array(sivel_renorm)
D = numpy.random.uniform(low=0.02, high=0.02, size=nsne)
k0 = numpy.random.uniform(low=0.2, high=0.2, size=nsne)
k1 = numpy.random.uniform(low=0.2, high=0.2, size=nsne)

mn = Delta[:,None] + alpha[None,:] * EW[:,0][:,None] + \
   beta[None,:] * EW[:,1][:,None] + eta[None,:]*sivel[:,None] + delta[None,:]*D[:,None]

intrinsic=[]
for i in xrange(nsne):
   e = numpy.random.normal(0,1,size=5)
   intrinsic.append(mn[i,:] + e*sigma)

intrinsic = intrinsic + gamma[None,:]*k0[:,None] + gamma1[None,:]*k1[:,None]
mag_obs=[]
EW_obs=[]
sivel_obs=[]
for i in xrange(nsne):
   mag_obs.append(numpy.random.multivariate_normal(intrinsic[i,:],mag_cov[i]))
   EW_obs.append(numpy.random.multivariate_normal(EW[i,:],EW_cov[i]))
   sivel_obs.append(numpy.random.normal(sivel[i],sivel_err[i]))

mag_obs=numpy.array(mag_obs)
EW_obs=numpy.array(EW_obs)
sivel_obs = numpy.array(sivel_obs)

data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_renorm, 'sivel_err': sivel_err}


output = open('simdata.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()