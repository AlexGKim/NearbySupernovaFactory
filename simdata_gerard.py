#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel
import sncosmo

# two color parameter model
numpy.random.seed(seed=1)

synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))


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
av = numpy.random.exponential(scale=0.1, size=nsne)
rv = numpy.random.lognormal(mean=0.916, sigma=0.1, size=nsne)

cov =numpy.array([[0.0038, 0.001,-0.0002,0,0.0003],[0.001,0.0011,0.0002,0.,-0.0006],[-0.0002,0.0002,0.0004,0,-0.0002],\
   [0,0,0,0.0002,0.0002],[0.0003,-0.0006,-0.0002,0.0002,0.002]])


EW = numpy.array(EW_renorm)
sivel = numpy.array(sivel_renorm)
D = numpy.random.uniform(low=0.02, high=0.02, size=nsne)
k0 = numpy.random.uniform(low=0.2, high=0.2, size=nsne)
k1 = numpy.random.uniform(low=0.2, high=0.2, size=nsne)

mn = Delta[:,None] + alpha[None,:] * EW[:,0][:,None] + \
   beta[None,:] * EW[:,1][:,None] + eta[None,:]*sivel[:,None] + delta[None,:]*D[:,None]

intrinsic=[]
for i in xrange(nsne):
    dust = sncosmo.F99Dust(r_v=rv[i])
    dust.set(ebv=av[i]/rv[i])
    x0 = numpy.random.lognormal(mean=0, sigma=0.05)
    x1 = numpy.random.lognormal(mean=0, sigma=0.2)
    c = numpy.random.lognormal(mean=0, sigma=0.02)
    model0 = sncosmo.Model(source='salt2')
    model0.parameters[2]=x0
    model0.parameters[3]=x1
    model0.parameters[4]=c
    model = sncosmo.Model(source='salt2', effects=[dust], effect_names=['host'], effect_frames=['rest'])
    model.parameters[2]=x0
    model.parameters[3]=x1
    model.parameters[4]=c
    intrinsic.append(numpy.random.multivariate_normal(mn[i,:]+2.5*numpy.log10(model0.bandflux(synbands,0.)/model.bandflux(synbands,0.)),cov))

intrinsic = numpy.array(intrinsic)
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

data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_obs, 'EW_obs': EW_obs, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_obs, 'sivel_err': sivel_err}


output = open('simdata_gerard.pkl','wb')
pickle.dump(data, output, protocol=2)
output.close()