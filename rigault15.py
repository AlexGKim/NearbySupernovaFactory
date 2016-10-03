#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle
import sivel

f = open('forAlex_tmp16sept2016_localhost_parameters.dat', 'r')
for i in xrange(1):
    f.readline()

data = [x.split() for x in f.readlines()]
f.close()
data=numpy.array(data)
rnames= data[:,0]
data=numpy.delete(data,0,axis=1)
rdata=data.astype('float')

f = open('temp15.pkl','rb')
(fit, _) = pickle.load(f)
f.close()


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel, sivel_err, _, _ = sivel.sivel(data)

head = ['lmass','lsfr','lssfr','gmass']


use = numpy.isfinite(sivel)

names= numpy.array(data['snlist'])[use]

i = numpy.intersect1d(names, rnames, assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')
indr = numpy.zeros(len(i),dtype='int')

for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    indr[j] = numpy.where(rnames == i[j])[0]


(x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'],(50,50-34,50+34),axis=0)
fig, axes = plt.subplots(nrows=4)
for i in xrange(4):
  axes[i].errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
  axes[i].set_ylabel(r'$A_{\delta I}$')
  axes[i].set_xlabel(head[i])
  axes[i].set_ylim((-0.1,0.03))
fig.subplots_adjust(hspace=.3)
fig.set_size_inches(8,11)
pp = PdfPages("output15/rigault.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# wefwe


# i=2
# plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
# plt.ylabel(r'$E_\delta(B-V)$')
# plt.xlabel(r'local sSFR')
# pp = PdfPages("output15/rigaultssfr.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

wm = numpy.where(rdata[indr,9] < 10)[0]
# print r"${:9.3f} \pm {:9.3f}$".format((fit['rho1'][:,4][:,None]*fit['R'][:,wm]).mean(),(fit['rho1'][:,1][:,None]*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'][:,inda[wm]],(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
low1 = numpy.sum(x/dx**2)/dx2
dlow1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 


wm = numpy.where(rdata[indr,9] > 10)[0]
# print r"${:9.3f} \pm {:9.3f}$".format((fit['rho1'][:,4][:,None]*fit['R'][:,wm]).mean(),(fit['rho1'][:,1][:,None]*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'][:,inda[wm]],(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
high1 = numpy.sum(x/dx**2)/dx2
dhigh1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 


i=3
(x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'],(50,50-34,50+34),axis=0)
plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.axhline(low1, xmax=10)
plt.axhline(low1+dlow1, xmax=10,linestyle='--')
plt.axhline(low1-dlow1, xmax=10,linestyle='--')
plt.axhline(high1, xmin=10)
plt.axhline(high1+dhigh1, xmin=10,linestyle='--')
plt.axhline(high1-dhigh1, xmin=10,linestyle='--')
plt.ylabel(r'$A_{\delta I}$')
plt.xlabel(r'Host Mass ($M_\cdot$')
plt.ylim((-0.1,0.03))
pp = PdfPages("output15/rigault3.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# i=3

# plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
# plt.ylabel(r'$E_\delta(B-V)$')
# plt.xlabel(r'$\log{(M/M_\odot)}$')
# plt.axhline(low1, xmax=10)
# plt.axhline(low1+dlow1, xmax=10,linestyle='--')
# plt.axhline(low1-dlow1, xmax=10,linestyle='--')
# plt.axhline(high1, xmin=10)
# plt.axhline(high1+dhigh1, xmin=10,linestyle='--')
# plt.axhline(high1-dhigh1, xmin=10,linestyle='--')
# pp = PdfPages("output15/rigault3.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()
wm = numpy.where(rdata[indr,9] < 10)[0]
low = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()
wm = numpy.where(rdata[indr,9] > 10)[0]
hig = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()

plt.hist([low,hig],20,label=['low mass','high mass'],normed=True,range=(-0.04,0.02))
plt.xlabel(r'$A_{\delta V}$')
plt.legend()
pp = PdfPages("output15/rigault2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# import scipy.stats
# wm = numpy.where(rdata[:,9] < 10)[0]
# low = numpy.median(fit['rho1'][:,4]*fit['R'][:,wm],axis=0)

# wm = numpy.where(rdata[:,9] > 10)[0]
# hig = numpy.median(fit['rho1'][:,4]*fit['R'][:,wm],axis=0)

# ans= scipy.stats.ks_2samp(hig,low)
# print ans

