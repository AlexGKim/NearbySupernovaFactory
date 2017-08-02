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



# statistics of first two vectors relative to third
def meanfrac(fit):
    perpfrac=[]
    perpabs=[]
    for i in xrange(fit['gamma'].shape[0]):
      norm_gamma = numpy.dot(fit['gamma'][i,:], fit['gamma'][i,:])
      norm_gamma1 = numpy.dot(fit['gamma1'][i,:], fit['gamma1'][i,:])
      cross = numpy.dot(fit['gamma'][i,:], fit['gamma1'][i,:])

      a = numpy.array([[norm_gamma,cross],[cross,norm_gamma1]])
      y=numpy.array([numpy.dot(fit['rho1'][i,:],fit['gamma'][i,:]), numpy.dot(fit['rho1'][i,:],fit['gamma1'][i,:])])
      ans_ = numpy.linalg.solve(a,y)
      ans = fit['rho1'][i,:]-ans_[1]*fit['gamma1'][i,:] - ans_[0]*fit['gamma'][i,:]
      perpfrac.append(numpy.sum(ans **2)/numpy.sum(fit['rho1'][i,:] **2))
      perpabs.append(numpy.sum(ans **2))

    perpfrac=numpy.array(perpfrac)
    perpfracmode = numpy.percentile(perpfrac,(50,50-34,50+34))

    perpabs = numpy.array(perpabs)
    perpabsmode = numpy.percentile(perpabs,50+34)
    # bins=numpy.arange(0,1,.002)
    # plt.hist(perpfrac,bins)
    # plt.show()
    # wefwe
    # H,edges  = numpy.histogram(perpfrac,bins=bins)
    # w = numpy.argmax(H)
    # w = numpy.unravel_index(w,H.shape)
    # perpfracmode=edges[w]

    # bins=numpy.arange(0,2,.02)
    # H,edges  = numpy.histogram(perpabs,bins=bins)
    # w = numpy.argmax(H)
    # w = numpy.unravel_index(w,H.shape)
    # perpabsmode=edges[w]

    # plt.hist(perpabs,bins=20)
    # plt.show()
    return perpfracmode, perpabsmode


simperc = []
simmedfrac = []
simmedabs=[]
tag="_null"
tag=""
f = open('temp25_sim'+tag+'0.pkl','rb')
(fit,_) = pickle.load(f)



# bins=numpy.arange(0,4,0.1)
# H,edges  = numpy.histogram(numpy.sum(fit['rho1']**2,axis=1),bins=bins)
# w = numpy.argmax(H)
# w = numpy.unravel_index(w,H.shape)
# simperc.append(edges[w])

# simperc.append(numpy.mean(numpy.sum(fit['rho1']**2,axis=1)*numpy.sum(fit['R']**2,axis=1)))

simperc.append(numpy.percentile(numpy.sum(fit['rho1']**2,axis=1),(50,50-34,50+34)))

dum=meanfrac(fit)
simmedfrac.append(dum[0])
simmedabs.append(dum[1])

mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])

for index in xrange(1,50):
    f = open('temp25_sim'+tag+'{}.pkl'.format(index),'rb')
    (fit,_) = pickle.load(f)

    # H,edges  = numpy.histogram(numpy.sum(fit['rho1']**2,axis=1),bins=bins)
    # w = numpy.argmax(H)
    # w = numpy.unravel_index(w,H.shape) 
    # simperc.append(edges[w])
    # simperc.append(numpy.mean(numpy.sum(fit['rho1']**2,axis=1)*numpy.sum(fit['R']**2,axis=1)))
    simperc.append(numpy.percentile(numpy.sum(fit['rho1']**2,axis=1),(50,50-34,50+34)))


    dum=meanfrac(fit)
    simmedfrac.append(dum[0])
    simmedabs.append(dum[1])

    mega_ = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])
    mega = numpy.concatenate((mega,mega_),axis=1)

mega = numpy.transpose(mega)

simperc=numpy.array(simperc)
simmedfrac= numpy.array(simmedfrac)
simmedabs= numpy.array(simmedabs)

f = open('temp25.pkl','rb')
(fit_input,_) = pickle.load(f)

dataperc = numpy.percentile(numpy.sum(fit_input['rho1']**2,axis=1),(50,50-34,50+34))
datamedfrac, datamedabs = meanfrac(fit_input)

asort = numpy.argsort(numpy.append(simperc[:,0],dataperc[0]))
resort = numpy.argsort(asort)

pp = PdfPages('output25_sim'+tag+'/pvalues_delta2.pdf')
plt.errorbar(resort[:-1],simperc[:,0],yerr=[simperc[:,0]-simperc[:,1],simperc[:,2]-simperc[:,0]],fmt='.',color='blue')
plt.errorbar(resort[-1],dataperc[0],yerr=[[dataperc[0]-dataperc[1]],[dataperc[2]-dataperc[0]]],fmt='.',color='red')
plt.xlim([-1,len(resort)])
plt.xlabel('Sorted index')
plt.ylabel(r'$\delta^2$')
plt.savefig(pp,format='pdf')
plt.clf()
pp.close()

pp = PdfPages('output25_sim'+tag+'/pvalues_relperp.pdf')
asort = numpy.argsort(numpy.append(simmedfrac[:,0],datamedfrac[0]))
resort = numpy.argsort(asort)

plt.errorbar(resort[:-1],simmedfrac[:,0],yerr=[simmedfrac[:,0]-simmedfrac[:,1],simmedfrac[:,2]-simmedfrac[:,0]],fmt='.',color='blue')
plt.errorbar(resort[-1],datamedfrac[0],yerr=[[datamedfrac[0]-datamedfrac[1]],[datamedfrac[2]-datamedfrac[0]]],fmt='.',color='red')
plt.xlim([-1,len(resort)])
plt.xlabel('Sorted index')
plt.ylabel(r'$\delta_\perp^2 / \delta^2$')
plt.savefig(pp,format='pdf')
plt.clf()
pp.close()


# bins=numpy.arange(0,4,0.1)
# H,edges  = numpy.histogram(numpy.sum(fit_input['rho1']**2,axis=1),bins=bins)
# w = numpy.argmax(H)
# w = numpy.unravel_index(w,H.shape)
# dataperc=edges[w]


# print delta

# plt.hist(numpy.sum(fit_input['rho1']**2,axis=1),log=True,bins=bins)
# plt.show()


# wefwe
# dataperc = numpy.mean(numpy.sum(fit_input['rho1']**2,axis=1)*numpy.sum(fit_input['R']**2,axis=1))


# print datamedfrac
# print simmedfrac


print (simperc< dataperc).sum(), simperc.shape[0], (simperc< dataperc).sum()*1./simperc.shape[0]
print (simmedfrac > datamedfrac).sum(), simmedfrac.shape[0]
print (simmedabs > datamedabs).sum(), simmedabs.shape[0]

# pp = PdfPages('output25_sim'+tag+'/pvalues_delta2.pdf')
# plt.hist(simperc)
# plt.axvline(dataperc)
# plt.xlabel(r'$\delta^2$')
# plt.savefig(pp,format='pdf')
# plt.clf()
# pp.close()


# pp = PdfPages('output25_sim'+tag+'/pvalues_relperp.pdf')
# plt.hist(simmedfrac,log=True)
# plt.axvline(datamedfrac)
# plt.xlabel(r'$\delta^2_\perp/\delta^2$')
# plt.savefig(pp,format='pdf')
# plt.clf()

# pp.close()
# pp = PdfPages('output25_sim'+tag+'/pvalues_perp.pdf')

# plt.hist(simmedabs)
# plt.axvline(datamedabs)
# plt.xlabel(r'$\delta^2_\perp$')
# plt.savefig(pp,format='pdf')
# plt.clf()
# pp.close()

# wefw


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

# var = fit['rho1'][:,1:]/fit['rho1'][:,0][:,None]
# figure = corner.corner(var,labels=[r"${\delta}_B/{\delta}_U$",r"${\delta}_V/{\delta}_U$",r"${\delta}_R/{\delta}_U$",r"${\delta}_I/{\delta}_U$"],truths=delta[1:]/delta[0], \
#     range=(0.9,0.9,0.9,0.9))
# pp = PdfPages('output25_sim'+tag+'/deltaratio_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()

# wefwe


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


H,edges  = numpy.histogramdd(numpy.transpose(mega[:,:,6]),bins=50)
w = numpy.argmax(H)
w = numpy.unravel_index(w,H.shape)
print numpy.array([edges[0][w[0]],edges[1][w[1]],edges[2][w[2]],edges[3][w[3]],edges[4][w[4]]])

wefwe

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
