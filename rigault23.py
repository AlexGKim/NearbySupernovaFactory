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


f = open('temp23.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel, sivel_err, _, _, _, _ = sivel.sivel(data)

head = ['lmass','lsfr','lssfr','gmass']


use = numpy.isfinite(sivel)

names= numpy.array(data['snlist'])[use]

# Figure out intersection between the two lists
i = numpy.intersect1d(names, rnames, assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')
indr = numpy.zeros(len(i),dtype='int')

for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    indr[j] = numpy.where(rnames == i[j])[0]

# the intrinsic parameter is A_I

(x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'],(50,50-34,50+34),axis=0)

# print the table of intrinsic parameters
# for x1,x2,x3,x4 in zip(names,x,xmin,xmax):
#   print x1,x2,x3,x4

# wefwe

fig, axes = plt.subplots(nrows=4)
for i in xrange(4):
  axes[i].errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
  axes[i].set_ylabel(r'$A_{\delta I}$')
  axes[i].set_xlabel(head[i])
  axes[i].set_ylim((-0.015,0.06))
fig.subplots_adjust(hspace=.3)
fig.set_size_inches(8,11)
pp = PdfPages("output23/rigault.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# wefwe


# i=2
# plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
# plt.ylabel(r'$E_\delta(B-V)$')
# plt.xlabel(r'local sSFR')
# pp = PdfPages("output23/rigaultssfr.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# Calculate the weighted mean of low and high mass hosts
wm = numpy.where(rdata[indr,9] < 10)[0]
x_=x[inda[wm]]
dx = (xmax[inda[wm]]-xmin[inda[wm]])/2.
dx2 = numpy.sum(1/dx**2)
low1 = numpy.sum(x_/dx**2)/dx2
dlow1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(low1, dlow1) 

wm = numpy.where(rdata[indr,9] > 10)[0]
x_=x[inda[wm]]
dx = (xmax[inda[wm]]-xmin[inda[wm]])/2.
dx2 = numpy.sum(1/dx**2)
high1 = numpy.sum(x_/dx**2)/dx2
dhigh1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(high1, dhigh1) 

i=3
# (x, xmin, xmax) = numpy.percentile(fit['rho1'][:,4][:,None]*fit['R'],(50,50-34,50+34),axis=0)
plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.plot([5,10],[low1,low1],color='black')
plt.plot([5,10],[low1+dlow1,low1+dlow1],linestyle='--',color='black')
plt.plot([5,10],[low1-dlow1,low1-dlow1],linestyle='--',color='black')
plt.plot([10,12],[high1,high1],color='black')
plt.plot([10,12],[high1+dhigh1,high1+dhigh1],linestyle='--',color='black')
plt.plot([10,12],[high1-dhigh1,high1-dhigh1],linestyle='--',color='black')
plt.ylabel(r'$A_{\delta I}$',fontsize=20)
plt.xlabel(r'Host Mass ($\log(M_\odot)$)',fontsize=20)
plt.ylim((-0.02,0.05))
plt.xlim((6,12))
pp = PdfPages("output23/rigault3.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


wm = numpy.where(rdata[indr,9] < 10)[0]
low = (fit['rho1'][:,4][:,None]*fit['R'][:,wm]).flatten()
wm = numpy.where(rdata[indr,9] > 10)[0]
hig = (fit['rho1'][:,4][:,None]*fit['R'][:,wm]).flatten()

plt.hist([low,hig],20,label=['low mass','high mass'],color=['blue','yellow'],normed=True,range=(-0.05,0.05))
plt.xlabel(r'$A_{\delta I}$',fontsize=20)
plt.xlim((-0.05,0.05))
plt.legend()
pp = PdfPages("output23/rigault2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wm = numpy.where(rdata[indr,6] < -11)[0]
low = (fit['rho1'][:,4][:,None]*fit['R'][:,wm]).flatten()
wm = numpy.where(rdata[indr,6] > -11)[0]
hig = (fit['rho1'][:,4][:,None]*fit['R'][:,wm]).flatten()

plt.hist([low,hig],20,label=['low lssfr','high lssfr'],normed=True,range=(-0.03,0.06))
plt.xlabel(r'$A_{\delta I}$')
plt.legend()
pp = PdfPages("output23/rigault4.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# Calculate the weighted mean of llsfr
wm = numpy.where(rdata[indr,6] < -11)[0]
x_=x[inda[wm]]
dx = (xmax[inda[wm]]-xmin[inda[wm]])/2.
dx2 = numpy.sum(1/dx**2)
low1 = numpy.sum(x_/dx**2)/dx2
dlow1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(low1, dlow1) 

wm = numpy.where(rdata[indr,6] > -11)[0]
x_=x[inda[wm]]
dx = (xmax[inda[wm]]-xmin[inda[wm]])/2.
dx2 = numpy.sum(1/dx**2)
high1 = numpy.sum(x_/dx**2)/dx2
dhigh1=  1/numpy.sqrt(dx2)
print r"${:9.4f} \pm {:9.4f}$".format(high1, dhigh1) 


import scipy.stats
wm = numpy.where(rdata[:,9] < 10)[0]
low = numpy.median(fit['rho1'][:,4][:,None]*fit['R'][:,wm],axis=0)

wm = numpy.where(rdata[:,9] > 10)[0]
hig = numpy.median(fit['rho1'][:,4][:,None]*fit['R'][:,wm],axis=0)

ans= scipy.stats.ks_2samp(hig,low)
print ans

