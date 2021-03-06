#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

rc('text', usetex=True)


f = open('apjl495662t1_mrt.txt', 'r')

for i in xrange(32):
    f.readline()

data = [x.split() for x in f.readlines()]
f.close()

filts = ['U','B','R','i']
elam = numpy.array([3656.,4353,6349,7625])
exv=[]
e_exv=[]
for f in filts:
    a1 = 0.
    a2 = 0.
    n = 0
    for d in data:
        if (f == d[2].strip() and numpy.abs(float(d[1]) <= 5)):
          err2 = float(d[4])**2 + float(d[8])**2
          a1 = a1 + (float(d[3]) - float(d[5])) - (float(d[7])-float(d[9]))+ float(d[10])
          a2 = a2 + err2
          n = n+1
    exv.append(a1/n)
    e_exv.append(numpy.sqrt(a2/n))

exv= numpy.array(exv)
e_exv = numpy.array(e_exv)

filts=['U','B','R','i']
for i in xrange(4):
    print "$E({}-V) = {:6.2f} \pm {:6.2f}$".format(filts[i],exv[i],e_exv[i])

f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

gammam1 = numpy.median(fit['gamma']/(fit['gamma'][:,2][:,None]) -1,axis=0)
gammam1= numpy.delete(gammam1,2)
deltam1 = numpy.median(fit['rho1']/(fit['rho1'][:,2] [:,None])-1,axis=0)
deltam1 = numpy.delete(deltam1,2)

fisher = numpy.array([[sum(gammam1**2/e_exv), sum(gammam1*deltam1/e_exv)],[sum(gammam1*deltam1/e_exv), sum(deltam1**2/e_exv)]])
u=numpy.linalg.inv(fisher)
g = numpy.array([sum(exv*gammam1/e_exv), sum(exv*deltam1/e_exv)])
ans = numpy.dot(u,g)
print 'best fit parameters'
print '${0[0]:6.2f}$, ${0[1]:6.2f}$'.format(ans)
print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in u])


fiterr = numpy.sqrt(gammam1**2 * u[0,0] + 2*deltam1*gammam1*u[0,1] + deltam1**2 * u[1,1])
ebvext= gammam1[1]*ans[0]
ebvint = deltam1[1]*ans[1]
print "$E(B-V)={:6.2f} \\pm {:6.2f}$, $E_\delta(B-V)={:6.2f} \\pm {:6.2f}$".format(ebvext,numpy.sqrt(u[0,0])*gammam1[1], ebvint,numpy.sqrt(u[1,1])*numpy.abs(deltam1[1]))


f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()
dum = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],(50,50-34,50+34),axis=0)
imax = numpy.argmin(dum,axis=1)[0]
print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
imax = numpy.argmax(dum,axis=1)[0]
print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
dum = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],(50,50-34,50+34),axis=0)
imax = numpy.argmin(dum,axis=1)[0]
print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
imax = numpy.argmax(dum,axis=1)[0]
print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])

import sncosmo
def lnprob(p, x, y ,yerr):
   # p is ebv and r_v
   f99 = sncosmo.F99Dust(r_v =p[1])
   f99.set(ebv=p[0])
   A_ = f99.propagate(x,1.) / f99.propagate(numpy.array([5287.48667023]),1.)[0]
   A_ = -2.5*numpy.log10(A_)
   dum = y-A_
   ans = -0.5*numpy.sum((dum/yerr)**2)
   return ans

efflam = numpy.array([ 3693.16777627,  4369.37505509,  6319.19906153,7610.89305298])
ndim, nwalkers = 2, 100
pos = [numpy.array([1.3, 1.4]) + 1e-4*numpy.random.randn(ndim) for i in range(nwalkers)]

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(efflam, exv, e_exv))
sampler.run_mcmc(pos,10)
samples = sampler.chain[:, 5:, :].reshape((-1, ndim))
(y,ymin,ymax)  = numpy.percentile(samples[:,1],(50,50-34,50+34))
ebv=y
print 'Fitzpatrick fit'
print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)
(y,ymin,ymax)  = numpy.percentile(samples[:,0],(50,50-34,50+34))
print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)
rvs=[y]

plt.errorbar(elam,exv,yerr=[e_exv,e_exv],fmt='.',label='SN2014J',color='black')
plt.errorbar(elam, gammam1*ans[0]+deltam1*ans[1],yerr=[fiterr,fiterr],label='Best-fit Model',color='red',fmt='.')


lambdas = numpy.arange(3500.,8000,100)
rvs=[1.4]
ebv=1.37

for rv in rvs:
    f99 = sncosmo.F99Dust(r_v =rv)
    f99.set(ebv=ebv)
    A_ = f99.propagate(lambdas,1.) / f99.propagate(numpy.array([5477]),1.)[0]
    A_=-2.5*numpy.log10(A_)
    # A_ = sncosmo._extinction.ccm89(lambdas, 1.37, rv) - sncosmo._extinction.ccm89(numpy.array([5477.]), 1.37, rv)[0]
    # norm  = sncosmo._extinction.ccm89(numpy.array([5477.]), 1., rv)
    # A_ = A_/norm[0]-1
    plt.plot(lambdas,A_,label=r"Amanullah et al. (2014)")

plt.legend()
plt.xlim((3200,8000))
plt.ylim((-2,3))
plt.ylabel(r'$E_o(X-V)$')
plt.xlabel(r'Wavelength (\AA)')
pp = PdfPages('output11/sn2014j.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
