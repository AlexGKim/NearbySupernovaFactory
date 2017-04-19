#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import f99_band
import matplotlib as mpl
mpl.rcParams['font.size'] = 16


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


# in fundamental units
tvec = [numpy.median(fit['gamma'][:,1]/fit['gamma'][:,2]-1,axis=0),numpy.median(fit['rho1'][:,1]/fit['rho1'][:,2]-1,axis=0)]
Egamma = ans* tvec
ugamma=u*numpy.outer(tvec,tvec)
# print numpy.outer(tvec,tvec)
# print Egamma
# print ugamma


fiterr = numpy.sqrt(gammam1**2 * u[0,0] + 2*deltam1*gammam1*u[0,1] + deltam1**2 * u[1,1])
ebvext= gammam1[1]*ans[0]
ebvint = deltam1[1]*ans[1]
# print "$E(B-V)={:6.2f} \\pm {:6.2f}$, $E_\delta(B-V)={:6.2f} \\pm {:6.2f}$".format(ebvext,numpy.sqrt(u[0,0])*gammam1[1], ebvint,numpy.sqrt(u[1,1])*numpy.abs(deltam1[1]))


dum = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],(50,50-34,50+34),axis=0)
imax = numpy.argmin(dum,axis=1)[0]
# print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
imax = numpy.argmax(dum,axis=1)[0]
# print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
dum = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],(50,50-34,50+34),axis=0)
imax = numpy.argmin(dum,axis=1)[0]
# print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])
imax = numpy.argmax(dum,axis=1)[0]
# print '${:6.2f}^{{{:6.2f}}}_{{{:6.2f}}}$'.format(dum[0][imax],dum[2][imax]-dum[0][imax],dum[1][imax]-dum[0][imax])




import sncosmo
# def lnprob(p, x, y ,yerr):
#    # p is ebv and r_v
#    f99 = sncosmo.F99Dust(r_v =p[1])
#    f99.set(ebv=p[0])
#    A_ = f99.propagate(x,1.) / f99.propagate(numpy.array([5287.48667023]),1.)[0]
#    A_ = -2.5*numpy.log10(A_)
#    dum = y-A_
#    ans = -0.5*numpy.sum((dum/yerr)**2)
#    return ans

# efflam = numpy.array([ 3693.16777627,  4369.37505509,  6319.19906153,7610.89305298])
# ndim, nwalkers = 2, 100
# pos = [numpy.array([1.3, 1.4]) + 1e-4*numpy.random.randn(ndim) for i in range(nwalkers)]

# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(efflam, exv, e_exv))
# sampler.run_mcmc(pos,10)
# samples = sampler.chain[:, 5:, :].reshape((-1, ndim))
# (y,ymin,ymax)  = numpy.percentile(samples[:,1],(50,50-34,50+34))
# ebv=y
# print 'Fitzpatrick fit'
# print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)
# (y,ymin,ymax)  = numpy.percentile(samples[:,0],(50,50-34,50+34))
# print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)
# rvs=[y]

fig, axes = plt.subplots(nrows=2,sharex=True,gridspec_kw = {'height_ratios':[3, 1]}, figsize=(8,8))

axes[0].errorbar(elam,exv,yerr=[e_exv,e_exv],fmt='.',label='SN2014J',color='black')
axes[0].errorbar(elam, gammam1*ans[0]+deltam1*ans[1],yerr=[fiterr,fiterr],label='Best-fit Model',color='red',fmt='.')


lambdas = numpy.arange(3500.,8000,100)
rvs=[1.4]
ebv=1.37

# for rv in rvs:
#     f99 = sncosmo.F99Dust(r_v =rv)
#     f99.set(ebv=ebv)
#     A_ = f99.propagate(lambdas,1.) / f99.propagate(numpy.array([5477]),1.)[0]
#     A_=-2.5*numpy.log10(A_)
#     # A_ = sncosmo._extinction.ccm89(lambdas, 1.37, rv) - sncosmo._extinction.ccm89(numpy.array([5477.]), 1.37, rv)[0]
#     # norm  = sncosmo._extinction.ccm89(numpy.array([5477.]), 1., rv)
#     # A_ = A_/norm[0]-1
#     axes[0].plot(lambdas,A_,label=r"Amanullah et al. (2014)")

axes[0].legend()
axes[0].set_xlim((3200,8000))
axes[0].set_ylim((-2,3))
axes[0].set_ylabel(r'$E_o(X-V)$')
axes[1].set_xlabel(r'Wavelength (\AA)')

plt.setp(axes[0].get_xticklabels(),visible=False)

axes[1].errorbar(elam,exv-(gammam1*ans[0]+deltam1*ans[1]), \
  yerr=[numpy.sqrt(e_exv**2+fiterr**2),numpy.sqrt(e_exv**2+fiterr**2)],fmt='.')
axes[1].axhline(0,linestyle=':')
axes[1].set_ylabel(r'$\Delta E_o(X-V)$')


plt.subplots_adjust(hspace=.07)
pp = PdfPages('output11/sn2014j.pdf')

plt.savefig(pp,format='pdf')
plt.tight_layout()

pp.close()
plt.close()

# Partial derivatives with respect to av and ebv
av=0.1
ebv=0.1/2.5
A1= f99_band.A_X(r_v=av/ebv, ebv=ebv)
A2= f99_band.A_X(r_v=(av+0.01)/ebv, ebv=ebv)
dAdAv = (A2 - A1)/0.01

A3= f99_band.A_X(r_v=av/(ebv+0.001), ebv=ebv+0.001)
dAdebv = (A3 - A1)/0.001

# print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdAv)
# print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdebv)


# The equation of interest is
# gammma0 = ans00 F0 + ans01 F1 + res
# gammma0 = ans10 F0 + ans11 F1 + res
# where F are the Fitzpatrick vectors (partial derivatives above) and
# the residues are perpendicular to a and b
# Note that the gammas are really gamma_X/(gamma_B-gamma_V)

norm_dAdebv = numpy.dot(dAdebv, dAdebv)
norm_dAdAv = numpy.dot(dAdAv, dAdAv)
cross = numpy.dot(dAdebv, dAdAv)

a = numpy.array([[norm_dAdebv,cross],[cross,norm_dAdAv]])

tmat = []
res = []
c_n = []
cs = []
for s in ['gamma','rho1']:
  c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
  # print "{:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}".format(c[0],c[1],c[2],c[3],c[4])
  cs.append(c)
  c_norm = numpy.linalg.norm(c)
  c_n.append(c_norm)

  y = numpy.array([numpy.dot(c,dAdebv),numpy.dot(c,dAdAv)])
  ans = numpy.linalg.solve(a,y)

  tmat.append(ans)
  ans = c-ans[0]*dAdebv - ans[1]*dAdAv
  res.append(ans)

tmat = numpy.array(tmat)
res= numpy.array(res)

#print the matrix and the residues

print numpy.dot(tmat.T,Egamma)

newmat=numpy.array(tmat)
newmat[:]=0
for i in xrange(2):
  for j in xrange(i,2):
    for k in xrange(2):
      for l in xrange(2):
        # print i, j ,k,l,tmat.T[i,k],tmat.T[j,l]
        newmat[i,j]= newmat[i,j]+tmat.T[i,k]*tmat.T[j,l]*ugamma[k,l]
    newmat[j,i]=newmat[i,j]
print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in newmat])

