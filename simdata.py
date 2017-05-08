#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel as sivel2
import flip
import matplotlib.pyplot as plt
import corner

# two color parameter model

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()


f = open('temp25.pkl','rb')
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

# c = numpy.zeros(5)
# alpha = numpy.array([0.005, 0.002,0.002,0.002,0.003])
# beta = numpy.array([0.03, 0.03,0.03,0.02,0.02])
# eta = numpy.array([0.0001,1.5e-4,0.0005,0.0005,-0.00003])
# gamma = numpy.array([64.1,51.9,38.5,29.3,20.8])
# gamma1 = numpy.array([-10.5,-11.9,-14.6,-13.9,-12.9])
# delta = numpy.array([2.13, -2.30,-1.87,0.38,2.09])
# delta = delta.transpose()
# y=numpy.array([numpy.dot(delta, gamma), numpy.dot(delta, gamma1)])
# A=numpy.zeros((2,2))
# A[0,0] = numpy.dot(gamma,gamma)
# A[1,1] = numpy.dot(gamma1,gamma1)
# A[0,1] = numpy.dot(gamma,gamma1)
# A[1,0] = A[0,1];
# x = numpy.dot(numpy.linalg.inv(A),y)
# delta= delta - x[0]*gamma - x[1]*gamma1;
# cov =numpy.array([[0.002, 0.0,-0.0,0,0.00],[0.000,0.0008,0.0000,0.,-0.0000],[-0.0000,0.0000,0.0003,0,-0.0000],\
#    [0,0,0,0.0002,0.0000],[0.0000,-0.0000,-0.0000,0.0000,0.001]])
c = numpy.median(fit['c'],axis=0)
alpha = numpy.median(fit['alpha'],axis=0)
beta = numpy.median(fit['beta'],axis=0)
eta = numpy.median(fit['eta'],axis=0)
gamma = numpy.median(fit['gamma'],axis=0)
gamma1 = numpy.median(fit['gamma1'],axis=0)

# EW = numpy.array(EW_renorm)
# sivel = numpy.array(sivel_renorm)
EW = numpy.median(fit['EW'],axis=0)
sivel = numpy.median(fit['sivel'],axis=0)

# for i in xrange(5):
#   plt.hist(fit['rho1'][:,i])
#   plt.show()


temprho1,tempR = flip.flip(fit)

H,edges  = numpy.histogramdd(temprho1,bins=30)
w = numpy.argmax(H)
w = numpy.unravel_index(w,H.shape)
rho1=numpy.array([edges[0][w[0]],edges[1][w[1]],edges[2][w[2]],edges[3][w[3]],edges[4][w[4]]])

print rho1

# R=[]
# R2=[]
# for i in xrange(tempR.shape[1]):
  # H,edges  = numpy.histogram(tempR[:,i],bins=20)
  # R.append(edges[numpy.argmax(H)])
  # R2.append(numpy.mean(temprho1*tempR[:,i][:,None]/rho1))
  # plt.hist(tempR[:,i],bins=20)
  # plt.show()

# R=numpy.array(R)
# R2=numpy.array(R2)

# plt.hist([R,R2,numpy.median(tempR,axis=0)],label=['mode','avg','median'])
# plt.legend()
# plt.show()

# wfwe

R = numpy.median(tempR,axis=0)

Delta = numpy.median(fit['Delta'],axis=0)
k0 = numpy.median(fit['k'],axis=0)
k1 = numpy.median(fit['k1'],axis=0)

cov = numpy.zeros((5,5))
for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
    cov= cov+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

cov/= len(fit['L_Omega'])

mn = Delta[:,None] + c[None,:] + alpha[None,:] * EW[:,0][:,None] + \
   beta[None,:] * EW[:,1][:,None] + eta[None,:]*sivel[:,None] + rho1[None,:]*R[:,None]

for seed in xrange(10):
    numpy.random.seed(seed=seed)


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

    # # renormalize the data
    EW_mn = EW_obs.mean(axis=0)
    EW_renorm = (EW_obs - EW_mn)

    mag_mn = mag_obs.mean(axis=0)
    mag_renorm  = mag_obs-mag_mn

    sivel_mn = sivel.mean()
    sivel_renorm = sivel-sivel_mn

    data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
       'sivel_obs': sivel_renorm, 'sivel_err': sivel_err}


    output = open('simdata{}.pkl'.format(seed),'wb')
    pickle.dump(data, output, protocol=2)
    output.close()