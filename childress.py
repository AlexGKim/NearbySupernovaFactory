#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle
import scipy.stats

f = open('MJC_compile_SNdata.pkl','r')
gal= pickle.load(f)

f = open('temp25.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

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

i = numpy.intersect1d(names, gal.keys(), assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')

mass=[]
emass = []
for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    emass.append(gal[i[j]]['eMass'])
    mass.append(gal[i[j]]['Mass'])

mass = numpy.array(mass)
emass= numpy.array(emass)

# scaled = fit['rho1'][:,0][:,None]*(fit['R']-fit['R'][:,0][:,None])
scaled = fit['rho1'][:,0][:,None]*fit['R']
scaled = scaled[:,inda]

ksarr=[]
pvarr=[]
wlow = numpy.where(mass < 10)[0]
wbig = numpy.where(mass >= 10)[0]
for i in xrange(scaled.shape[0]):
  temp1,temp2= scipy.stats.ks_2samp(scaled[i,wlow],scaled[i,wbig])
  print temp1,temp2
  ksarr.append(temp1)
  pvarr.append(temp2)
ksarr=numpy.array(ksarr)
pvarr=numpy.array(pvarr)

plt.hist(pvarr)
plt.xlabel(r'$p$-value')

pp = PdfPages("output25/childress_pvalue.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(x, xmin, xmax) = numpy.percentile(scaled,(50,50-34,50+34),axis=0)

plt.errorbar(mass,x,xerr=[emass, emass], yerr=[x-xmin,xmax-x],fmt='o')
plt.ylabel(r'$A_{\delta U}-A_{\delta U}|_0$')
plt.xlabel(r'$\log{(M_{host}/M_{\odot})}$')

pp = PdfPages("output25/childress.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

done

ux = numpy.array([6,10])
wm = numpy.where(mass < 10)[0]
temp = scaled[:,wm].flatten()
temp=temp[temp !=0]
(x, xmin, xmax) = numpy.percentile(temp,(50,50-34,50+34))

print r"${:9.4f}^{{ {:9.4f} }}_{{ {:9.4f} }}$".format(x,xmax-x,x-xmin) 
plt.plot(ux, [x,x],color='black')
plt.plot(ux, [xmax,xmax],color='red')
plt.plot(ux, [xmin,xmin],color='red')


ux = numpy.array([10,13])
wm = numpy.where(mass[1:] > 10)[0]
temp = scaled[:,wm].flatten()
temp=temp[temp !=0]
(x, xmin, xmax) = numpy.percentile(temp,(50,50-34,50+34))
print r"${:9.4f}^{{ {:9.4f} }}_{{ {:9.4f} }}$".format(x,xmax-x,x-xmin) 
plt.plot(ux, [x,x],color='black')
plt.plot(ux, [xmax,xmax],color='red')
plt.plot(ux, [xmin,xmin],color='red')

plt.xlim((6,13))
# plt.ylim((-0.1,0.1))
pp = PdfPages("output25/childress.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wm = numpy.where(mass[1:] < 10)[0]
low = (((fit['rho1'][:,0])[:,None])*fit['R'][:,wm]).flatten()
lowlim = numpy.percentile(scaled,(50,50-34,50+34),axis=0)
lowmn = lowlim[2,:]+lowlim[1,:]/2
lowsig = lowlim[2,:]-lowlim[1,:]/2
print "${:6.1e} \pm {:6.1e}$".format(numpy.sum(lowmn/lowsig**2)/numpy.sum(1/lowsig**2), 1/numpy.sqrt(numpy.sum(1/lowsig**2)))
wm = numpy.where(mass[1:] > 10)[0]
hig = (((fit['rho1'][:,0])[:,None])*fit['R'][:,wm]).flatten()
higlim = numpy.percentile(scaled,(50,50-34,50+34),axis=0)
higmn = higlim[2,:]+higlim[1,:]/2
higsig = higlim[2,:]-higlim[1,:]/2
print "${:6.1e} \pm {:6.1e}$".format(numpy.sum(higmn/higsig**2)/numpy.sum(1/higsig**2), 1/numpy.sqrt(numpy.sum(1/higsig**2)))



# import scipy.stats
# ans= scipy.stats.ks_2samp(hig,low)


# print '{:6.2e} {:6.2e}'.format(low.mean(), low.std())
# print '{:6.2e} {:6.2e}'.format(hig.mean(), hig.std())


plt.hist([lowlim[1],higlim[1]],label=['low mass','high mass'],normed=True,bins=20)
plt.xlabel(r'$A_{\delta U}-A_{\delta U}|_0$')
plt.legend()
pp = PdfPages("output25/childress2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



