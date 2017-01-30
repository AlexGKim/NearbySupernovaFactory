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


mpl.rcParams['font.size'] = 18


ext=''

f = open('temp18'+ext+'.pkl','rb')
(fit18,_) = pickle.load(f)
f = open('temp18'+ext+'_free.pkl','rb')
(fit18_free,_) = pickle.load(f)
f = open('temp18'+ext+'_salt2.pkl','rb')
(fit18_salt2,_) = pickle.load(f)
f = open('temp19'+ext+'.pkl','rb')
(fit19,_) = pickle.load(f)
f = open('temp20'+ext+'.pkl','rb')
(fit20,_) = pickle.load(f)
f = open('temp20'+ext+'_salt2.pkl','rb')
(fit20_salt2,_) = pickle.load(f)
f = open('temp20_ccm.pkl','rb')
(fit20_ccm,_) = pickle.load(f)
f = open('temp11'+ext+'.pkl','rb')
(fit11,_) = pickle.load(f)




plt.hist(fit18_salt2['AV'].flatten(),histtype='step', stacked='false', \
    label=r'$R=const$ SALT2',bins=numpy.arange(-.25,1.8,0.1),normed=True)
plt.hist(fit18['AV'].flatten(),histtype='step', stacked='false', \
    label=r'$R=const$ Hsiao',bins=numpy.arange(-.25,1.8,0.1),normed=True)
# plt.hist(fit18_free['AV'].flatten(),histtype='step', stacked='false', \
#     label=r'$R=const$ free',bins=numpy.arange(-.25,1.8,0.1),normed=True)
# plt.hist(fit19['AV'].flatten(),histtype='step', stacked='false', \
#     label=r'$\ln{R}\sim \mathcal{N}$',bins=numpy.arange(-.25,1.8,0.1),normed=True)
# plt.hist(fit20['AV'].flatten(),histtype='step', stacked='false', \
#     label=r'$R$-free',bins=numpy.arange(-.25,1.8,0.1),normed=True)
# plt.hist(fit20_salt2['AV'].flatten(),histtype='step', stacked='false', \
#     label=r'$R$-free S2',bins=numpy.arange(-.25,1.8,0.1),normed=True)

plt.xlabel(r'$A_V$')
plt.legend()
plt.tight_layout()
pp = PdfPages('output18'+ext+'/AVs_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(numpy.median(fit18_salt2['AV'],axis=0),label=r'$R=const$ SALT2',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
plt.hist(numpy.median(fit18['AV'],axis=0),label=r'$R=const$ Hsiao',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
# plt.hist(numpy.median(fit19['AV'],axis=0),label=r'$\ln{R}\sim \mathcal{N}$',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
# plt.hist(numpy.median(fit20['AV'],axis=0),label=r'$R$-free',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
# plt.hist(numpy.median(fit20_salt2['AV'],axis=0),label=r'$R$-free S2',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
# plt.hist(numpy.median(fit20_ccm['AV'],axis=0),label=r'$R$-free CCM',histtype='step',bins=numpy.arange(-0.2,1.8,0.1))
plt.hist(numpy.median(fit11['gamma'][:,2][:,None]*fit11['k'],axis=0),label=r'$\gamma_0 k_0$',histtype='step',bins=numpy.arange(-0.2,1.8,0.1),normed=True)
plt.xlabel(r'$A_V$')
plt.legend()
plt.tight_layout()
pp = PdfPages('output18'+ext+'/AVs_mode_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# plt.hist(numpy.median(fit20_ccm['AV']-fit20['AV'],axis=0),histtype='step')
# # plt.hist(numpy.median(fit20_ccm['AV'],axis=0),label=r'$R$-free CCM',histtype='step',bins=numpy.arange(-0.2,1.8,0.1))
# plt.xlabel(r'$A_V(CCM) - A_V(F99)$')
# pp = PdfPages('output18'+ext+'/AVs_mode_hist_comp.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


(ymin, ymax) = numpy.percentile(fit18['AV'],(50-34,50+34),axis=0)
f18perc=(ymax-ymin)/2
(ymin, ymax) = numpy.percentile(fit19['AV'],(50-34,50+34),axis=0)
f19perc=(ymax-ymin)/2
(ymin, ymax) = numpy.percentile(fit20['AV'],(50-34,50+34),axis=0)
f20perc=(ymax-ymin)/2
(ymin, ymax) = numpy.percentile(fit20_ccm['AV'],(50-34,50+34),axis=0)
f20_ccmperc=(ymax-ymin)/2
(ymin, ymax) = numpy.percentile(fit11['gamma'][:,2][:,None]*fit11['k'],(50-34,50+34),axis=0)
f11perc=(ymax-ymin)/2

print "{:6.2f} {:6.2f} {:6.2f} {:6.2f}".format(f18perc.mean(),f19perc.mean(),f20perc.mean(),f20_ccmperc.mean(),f11perc.mean())


n=10000
choice = numpy.random.randint(0,len(fit19['lnRV_mn']),size=n)
raw = numpy.random.normal(0,1,size=n)
fit19RV = numpy.exp(fit19['lnRV_mn'][choice] + raw*fit19['lnRV_sig'][choice])
shit=numpy.percentile(1./fit18_salt2['RVinv'],(50,50-34,50+34))
print shit[0],shit[0]-shit[1],shit[2]-shit[0]
shit= numpy.percentile(1./fit18['RVinv'],(50,50-34,50+34))
print shit[0],shit[0]-shit[1],shit[2]-shit[0]
plt.hist(1./fit18_salt2['RVinv'].flatten(),histtype='step', stacked='false', \
    label=r'$R=const$ SALT2',bins=numpy.arange(2,4.2,0.05),normed=True)
plt.hist(1./fit18['RVinv'].flatten(),histtype='step', stacked='false', \
    label=r'$R=const$ Hsiao',bins=numpy.arange(2,4.2,0.05),normed=True)
# plt.hist(1./fit18_free['RVinv'].flatten(),histtype='step', stacked='false', \
#     label=r'$R=const$ free',bins=numpy.arange(0,8.2,0.1),normed=True)
# plt.hist(fit19RV,histtype='step', stacked='false', \
#     label=r'$\ln{R}\sim \mathcal{N}$',bins=numpy.arange(0,8.2,0.1),normed=True)
# plt.hist(numpy.median(fit20['RV'],axis=0),histtype='step', stacked='false', \
#     label=r'$R$-free',bins=numpy.arange(0,8.2,0.25),normed=True)
# plt.hist(numpy.median(fit20_salt2['RV'],axis=0),histtype='step', stacked='false', \
#     label=r'$R$-free S2',bins=numpy.arange(0,8.2,0.25),normed=True)
# plt.hist(numpy.median(fit20_ccm['RV'],axis=0),histtype='step', stacked='false', \
#     label=r'$R$-free CCM',bins=numpy.arange(0,8.2,0.1),normed=True)
plt.xlabel(r'$R_V$')
plt.legend()
plt.tight_layout()
pp = PdfPages('output18'+ext+'/RVs_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



# plt.hist(numpy.median(fit20_ccm['RV']-fit20['RV'],axis=0),histtype='step', stacked='false',normed=True,color='red')
# plt.xlabel(r'$R_V(CCM)-R_V(F99)$')
# pp = PdfPages('output18'+ext+'/RVs_hist_comp.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.hist(numpy.median(fit20_ccm['AV']-fit20['AV'],axis=0),histtype='step')
# # plt.hist(numpy.median(fit20_ccm['AV'],axis=0),label=r'$R$-free CCM',histtype='step',bins=numpy.arange(-0.2,1.8,0.1))
# plt.xlabel(r'$A_V(CCM) - A_V(F99)$')
# pp = PdfPages('output18'+ext+'/AVs_mode_hist_comp.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


(y, ymin, ymax) = numpy.percentile(fit20['RV'],(50,50-34,50+34),axis=0)
(x, xmin, xmax) = numpy.percentile(fit20['AV'],(50,50-34,50+34),axis=0)


plt.errorbar(x, y, xerr=[x-xmin,xmax-xmin],yerr=[y-ymin,ymax-y],fmt='.',alpha=0.15,color='b')
plt.scatter(x, y,marker='o',alpha=0.3,s=4,c='b')
# plt.xlim((-0.5,0.75))
# plt.ylim((-4,6))
plt.xlabel(r'$A_V^F$')
plt.ylabel(r'$R_V^F$')
plt.tight_layout()
pp = PdfPages('output18'+ext+'/AVRV.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# (y, ymin, ymax) = numpy.percentile(fit20_ccm['RV']-fit20['RV'],(50,50-34,50+34),axis=0)
# (x, xmin, xmax) = numpy.percentile(fit20_ccm['AV']-fit20['AV'],(50,50-34,50+34),axis=0)

# dx = (xmax-xmin)/2
# print numpy.sum(x/dx**2)/numpy.sum(1/dx**2)
# print numpy.sqrt(1/numpy.sum(1/dx**2))

# dy = (ymax-ymin)/2
# print numpy.sum(y/dy**2)/numpy.sum(1/dy**2)
# print numpy.sqrt(1/numpy.sum(1/dy**2))

# plt.errorbar(x, y, xerr=[x-xmin,xmax-xmin],yerr=[y-ymin,ymax-ymin],fmt='.',alpha=0.15,color='b')
# plt.scatter(x, y,marker='o',alpha=0.3,s=4,c='b')
# plt.xlim((-0.5,0.75))
# plt.ylim((-4,6))
# plt.xlabel(r'$A_V^C - A_V^F$')
# plt.ylabel(r'$R_V^C-R_V^F$')
# pp = PdfPages('output18'+ext+'/CCMF99.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()