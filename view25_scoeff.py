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
import flip


simperc = []

# mpl.rcParams['font.size'] = 18
tag="_null"
f = open('temp25_sim'+tag+'0.pkl','rb')
(fit,_) = pickle.load(f)

simperc.append(numpy.percentile(numpy.sum(fit['rho1']**2,axis=1),(68,90,95)))

mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])

for index in xrange(1,20):
    f = open('temp25_sim'+tag+'{}.pkl'.format(index),'rb')
    (fit,_) = pickle.load(f)
    simperc.append(numpy.percentile(numpy.sum(fit['rho1']**2,axis=1),(68,90,95)))

    mega_ = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])

    mega = numpy.concatenate((mega,mega_),axis=1)

mega = numpy.transpose(mega)

simperc=numpy.array(simperc)


f = open('temp25.pkl','rb')
(fit_input,_) = pickle.load(f)

dataperc = numpy.percentile(numpy.sum(fit_input['rho1']**2,axis=1),(68,90,95))

print (simperc[:,1]< dataperc[1]).sum(), simperc.shape[0], (simperc[:,1]< dataperc[1]).sum()*1./simperc.shape[0]
wefwef

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


temprho1,_ = flip.flip(fit_input)


H,edges  = numpy.histogramdd(temprho1,bins=30)
w = numpy.argmax(H)
w = numpy.unravel_index(w,H.shape)
delta=numpy.array([edges[0][w[0]],edges[1][w[1]],edges[2][w[2]],edges[3][w[3]],edges[4][w[4]]])
numpy.set_printoptions(precision=2)
print delta

if tag == '_null':
    delta[:]=0

cov = numpy.zeros((5,5))
for x1, x2 in zip(fit_input['L_Omega'], fit_input['L_sigma']):
    cov= cov+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

cov/= len(fit_input['L_Omega'])

for index in xrange(5):
    figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(filts[index]), r"$\alpha_{}$".format(filts[index]),\
                    r"$\beta_{}$".format(filts[index]),r"$\eta_{}$".format(filts[index]),r"$\gamma^0_{}$".format(filts[index]),\
                    r"$\gamma^1_{}$".format(filts[index]),r"$\delta_{{{}}}$".format(filts[index]), r"$\sigma_{}$".format(filts[index])],\
                    truths=[c[index], alpha[index],beta[index],eta[index],gamma[index], gamma1[index],delta[index],numpy.sqrt(cov[index,index])],label_kwargs={'fontsize':22})
    figure.suptitle(filts[index],fontsize=28)

    for ax in figure.get_axes():
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)

    pp = PdfPages('output25_sim'+tag+'/coeff{}.pdf'.format(index))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()


# figure = corner.corner(fit['gamma'],labels=[r"${\gamma^0}_0$",r"${\gamma^0}_1$",r"${\gamma^0}_2$",r"${\gamma^0}_3$",r"${\gamma^0}_4$"],label_kwargs={'fontsize':22})
# pp = PdfPages('output25_sim'+tag+'/gamma_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# figure = corner.corner(fit['gamma1'],labels=[r"${\gamma^1}_0$",r"${\gamma^1}_1$",r"${\gamma^1}_2$",r"${\gamma^1}_3$",r"${\gamma^1}_4$"])
# pp = PdfPages('output25_sim'+tag+'/gamma1_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()

# doublerho = numpy.concatenate((fit['rho1'],-fit['rho1']),axis=0)


H,edges  = numpy.histogramdd(numpy.transpose(mega[:,:,6]),bins=30)
w = numpy.argmax(H)
w = numpy.unravel_index(w,H.shape)
print numpy.array([edges[0][w[0]],edges[1][w[1]],edges[2][w[2]],edges[3][w[3]],edges[4][w[4]]])

perc = []
for i in xrange(10,25):
    H,edges  = numpy.histogramdd(numpy.transpose(mega[:,:,6]),bins=i)
    indeces=numpy.zeros(5,dtype='int')
    for ind in xrange(5):
      indeces[ind] = numpy.digitize(0,edges[ind])
    indeces=indeces-1
    # print indeces
    # print H.sum()
    # print H[indeces[0],indeces[1],indeces[2],indeces[3],indeces[4]]
    wbig = H > H[indeces[0],indeces[1],indeces[2],indeces[3],indeces[4]]
    print 'confidence non-zero ',H[wbig].sum() / H.sum(), i,  H[indeces[0],indeces[1],indeces[2],indeces[3],indeces[4]]
    # print 'confidence non-zero ',1./(H >0).sum(), (H >0).sum()
    perc.append(H[wbig].sum() / H.sum())

figure = corner.corner(numpy.transpose(mega[:,:,6]),labels=[r"${\delta}_U$",r"${\delta}_B$",r"${\delta}_V$",r"${\delta}_R$",r"${\delta}_I$"], \
    truths=-delta)
pp = PdfPages('output25_sim'+tag+'/delta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
# plt.close()


# import flip

# temprho1,_ = flip.flip(fit)

# figure = corner.corner(temprho1,labels=[r"${\delta}_{U}$",r"${\delta}_{B}$",r"${\delta}_{V}$",r"${\delta}_{R}$",r"${\delta}_{I}$"],\
#   truths=[0,0,0,0,0])
# pp = PdfPages('output25_sim'+tag+'/delta_corner_flipped.pdf')
# plt.tight_layout()
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# delta_ratio = fit_input['rho1']/fit_input['rho1'][:,4][:,None]


# figure = corner.corner(delta_ratio[:,0:4],labels=[r'\delta_U/\delta_I',r'\delta_B/\delta_I',r'\delta_V/\delta_I',r'\delta_R/\delta_I'])
# pp = PdfPages('output25_sim'+tag+'/delta_ratio.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()
