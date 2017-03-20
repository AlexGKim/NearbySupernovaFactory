#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel as sivel2

for seed in xrange(10):
    numpy.random.seed(seed=seed)

    # two color parameter model

    pkl_file = open('gege_data.pkl', 'r')
    data = pickle.load(pkl_file)
    pkl_file.close()


    f = open('temp23.pkl','rb')
    (fit,_) = pickle.load(f)
    f.close()

    sivel,sivel_err,x1,x1_err,zcmb,zerr = sivel2.sivel(data)



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

    c = numpy.zeros(5)
    alpha = numpy.array([0.005, 0.002,0.002,0.002,0.003])
    beta = numpy.array([0.03, 0.03,0.03,0.02,0.02])
    eta = numpy.array([0.0001,1.5e-4,0.0005,0.0005,-0.00003])
    gamma = numpy.array([64.1,51.9,38.5,29.3,20.8])
    gamma1 = numpy.array([-10.5,-11.9,-14.6,-13.9,-12.9])

    delta = numpy.array([2.13, -2.30,-1.87,0.38,2.09])
    delta = delta.transpose()
    y=numpy.array([numpy.dot(delta, gamma), numpy.dot(delta, gamma1)])
    A=numpy.zeros((2,2))
    A[0,0] = numpy.dot(gamma,gamma)
    A[1,1] = numpy.dot(gamma1,gamma1)
    A[0,1] = numpy.dot(gamma,gamma1)
    A[1,0] = A[0,1];
    x = numpy.dot(numpy.linalg.inv(A),y)
    print delta
    delta= delta - x[0]*gamma - x[1]*gamma1;
    print delta


    cov =numpy.array([[0.002, 0.0,-0.0,0,0.00],[0.000,0.0008,0.0000,0.,-0.0000],[-0.0000,0.0000,0.0003,0,-0.0000],\
       [0,0,0,0.0002,0.0000],[0.0000,-0.0000,-0.0000,0.0000,0.001]])


    EW = numpy.array(EW_renorm)
    sivel = numpy.array(sivel_renorm)

    Delta = numpy.median(fit['Delta'],axis=0)
    D = numpy.median(fit['R'],axis=0)
    k0 = numpy.median(fit['k'],axis=0)
    k1 = numpy.median(fit['k1'],axis=0)
    # D = numpy.random.normal(0,0.002, size=nsne)+numpy.random.exponential(0.0045, size=nsne)
    # k0 = numpy.random.normal(0,0.001, size=nsne)+numpy.random.exponential(0.004, size=nsne)
    # k1 = numpy.random.normal(0,0.002, size=nsne)+numpy.random.exponential(0.003, size=nsne)

    mn = Delta[:,None] + alpha[None,:] * EW[:,0][:,None] + \
       beta[None,:] * EW[:,1][:,None] + eta[None,:]*sivel[:,None] + delta[None,:]*D[:,None]

    intrinsic=[]
    for i in xrange(nsne):
       intrinsic.append(numpy.random.multivariate_normal(mn[i,:],cov))

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

    data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_obs, 'EW_obs': EW_obs, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
       'sivel_obs': sivel_obs, 'sivel_err': sivel_err}


    output = open('simdata{}.pkl'.format(seed),'wb')
    pickle.dump(data, output, protocol=2)
    output.close()