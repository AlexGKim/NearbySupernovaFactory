#!/usr/bin/env python

import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo
import scipy
import cPickle
import matplotlib as mpl
import sivel

f = open('temp23_sim0.pkl','rb')
(fit,_) = pickle.load(f)

mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])

print mega.shape
for index in xrange(1,10):
    f = open('temp23_sim{}.pkl'.format(index),'rb')
    (fit,_) = pickle.load(f)
    mega_ = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])

    mega = numpy.concatenate((mega,mega_),axis=1)

mega = numpy.transpose(mega)



f = open('temp23.pkl','rb')
(fit_input,_) = pickle.load(f)


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel,sivel_err,x1,x1_err,zcmb,zerr = sivel.sivel(data)

use = numpy.isfinite(sivel)


filts = ['U','B','V','R','I']

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

EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn

nsne, nmags = mag_obs.shape
color_obs = numpy.zeros((nsne,nmags-1))
color_obs[:,0] = mag_renorm[:,0]- mag_renorm[:,2]
color_obs[:,1] = mag_renorm[:,1]- mag_renorm[:,2]
color_obs[:,2] = mag_renorm[:,3]- mag_renorm[:,2]
color_obs[:,3] = mag_renorm[:,4]- mag_renorm[:,2]

EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

pkl_file.close()


trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))

filts = ['U','B','V','R','I']
c = numpy.median(fit_input['c'],axis=0)
alpha = numpy.median(fit_input['alpha'],axis=0)
beta = numpy.median(fit_input['beta'],axis=0)
eta = numpy.median(fit_input['eta'],axis=0)
gamma = numpy.median(fit_input['gamma'],axis=0)
gamma1 = numpy.median(fit_input['gamma1'],axis=0)
delta = numpy.median(fit_input['rho1'],axis=0)


cov = numpy.zeros((5,5))
for x1, x2 in zip(fit_input['L_Omega'], fit_input['L_sigma']):
    cov= cov+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

cov/= len(fit_input['L_Omega'])

for index in xrange(5):
    figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index), r"$\alpha_{}$".format(index),\
                    r"$\beta_{}$".format(index),r"$\eta_{}$".format(index),r"$\gamma^0_{}$".format(index),\
                    r"$\gamma^1_{}$".format(index),r"$\delta_{{{}}}$".format(index), r"$\sigma_{}$".format(index)],\
                    truths=[c[index], alpha[index],beta[index],eta[index],gamma[index], gamma1[index],delta[index],numpy.sqrt(cov[index,index])])
    figure.suptitle(filts[index],fontsize=28)
    pp = PdfPages('output23_sim/coeff{}.pdf'.format(index))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()



