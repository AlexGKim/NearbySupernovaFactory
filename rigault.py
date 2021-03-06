#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle

f = open('forAlex_tmp16sept2016_localhost_parameters.dat', 'r')
for i in xrange(1):
    f.readline()

data = [x.split() for x in f.readlines()]
f.close()
data=numpy.array(data)
rnames= data[:,0]
data=numpy.delete(data,0,axis=1)
rdata=data.astype('float')

f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

head = ['lmass','lsfr','lssfr','gmass']

dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

dic_meta=cPickle.load(open("META.pkl"))

sivel=[]
sivel_err=[]
for sn in data['snlist']:
   if sn in dic_meta.keys() and sn in dic_phreno.keys():
      meta = dic_meta[sn]
      vSiII_6355_lbd=0.
      vSiII_6355_lbd_err=0.
      counter  = 0
      for sp in dic_phreno[sn]["spectra"]:
         if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase']) < 2.5  and numpy.isfinite(dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]):
            vSiII_6355_lbd += dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
            vSiII_6355_lbd_err += 1/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
            counter +=1
      if counter !=0:
         sivel.append(vSiII_6355_lbd / vSiII_6355_lbd_err)
         sivel_err.append(1./numpy.sqrt(vSiII_6355_lbd_err))
      else:
         sivel.append(float('nan'))
         sivel_err.append(float('nan'))
   else:
      sivel.append(float('nan'))
      sivel_err.append(float('nan'))

sivel = numpy.array(sivel)
sivel_err = numpy.array(sivel_err)

use = numpy.isfinite(sivel)

names= numpy.array(data['snlist'])[use]

i = numpy.intersect1d(names, rnames, assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')
indr = numpy.zeros(len(i),dtype='int')

for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    indr[j] = numpy.where(rnames == i[j])[0]


(x, xmin, xmax) = numpy.percentile(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'],(50,50-34,50+34),axis=0)
fig, axes = plt.subplots(nrows=4)
for i in xrange(4):
  axes[i].errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
  axes[i].set_ylabel(r'$E_\delta(B-V)$')
  axes[i].set_xlabel(head[i])
fig.subplots_adjust(hspace=.3)
fig.set_size_inches(8,11)
pp = PdfPages("output11/rigault.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

i=3
plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.ylabel(r'$E_\delta(B-V)$')
plt.xlabel(r'$\log{(M/M_\odot)}$')
pp = PdfPages("output11/rigault3.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

i=2
plt.errorbar(rdata[indr,3*i],x[inda],xerr=[rdata[indr,3*i+2], rdata[indr,3*i+1]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.ylabel(r'$E_\delta(B-V)$')
plt.xlabel(r'local sSFR')
pp = PdfPages("output11/rigaultssfr.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wm = numpy.where(rdata[:,9] < 10)[0]
print r"${:9.3f} \pm {:9.3f}$".format((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).mean(),(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]),(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 
low = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()


wm = numpy.where(rdata[:,9] > 10)[0]
print r"${:9.3f} \pm {:9.3f}$".format((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).mean(),(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]),(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 
hig = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()


plt.hist([low,hig],label=['low mass','high mass'],normed=True,range=(-0.03,0.1))
plt.xlabel(r'$E_\delta(B-V)$')
plt.legend()
pp = PdfPages("output11/rigault2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


import scipy.stats
wm = numpy.where(rdata[:,9] < 10)[0]
low = numpy.median(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm],axis=0)

wm = numpy.where(rdata[:,9] > 10)[0]
hig = numpy.median(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm],axis=0)

ans= scipy.stats.ks_2samp(hig,low)
print ans

