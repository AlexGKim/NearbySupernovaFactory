#!/usr/bin/env python

import pickle
import numpy
import sncosmo
import scipy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def getFitzExt(efflam,av, ebv):
  f99 = sncosmo.F99Dust(r_v =av/ebv)
  f99.set(ebv=ebv)
  A_ = f99.propagate(efflam,1.)
  A1_ = -2.5*numpy.log10(A_)
  return A1_


efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])

f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

av=0.1
ebv=0.1/2.27
A1=getFitzExt(efflam, av , ebv)

A2=getFitzExt(efflam, av+0.01 , ebv)
dAdAv = (A2 - A1)/0.01

A3=getFitzExt(efflam, av , ebv+0.001)
dAdebv = (A3 - A1)/0.001

# dAdebv_norm = numpy.linalg.norm(dAdebv)
# dAdAv_norm = numpy.linalg.norm(dAdAv)

norm_dAdebv = numpy.dot(dAdebv, dAdebv)
norm_dAdAv = numpy.dot(dAdAv, dAdAv)
cross = numpy.dot(dAdebv, dAdAv)

a = numpy.array([[norm_dAdebv,cross],[cross,norm_dAdAv]])


# dAdebv = dAdebv/dAdebv_norm
# dAdAv = dAdAv/dAdAv_norm

tmat = []
res = []
c_n = []
cs = []
for s in ['gamma','rho1']:
  c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
  cs.append(c)
  c_norm = numpy.linalg.norm(c)
  c_n.append(c_norm)
  # c=c/c_norm

  y = numpy.array([numpy.dot(c,dAdebv),numpy.dot(c,dAdAv)])
  ans = numpy.linalg.solve(a,y)

  tmat.append(ans)
  ans = c-ans[0]*dAdebv - ans[1]*dAdAv
  res.append(ans)

tmat = numpy.array(tmat)
res= numpy.array(res)

print tmat
tmat = numpy.transpose(tmat)
print numpy.linalg.norm(res,axis=1)/numpy.array(c_n)


ebv  = ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']+((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']))
ebv = numpy.array([ebv,((fit['gamma'][:,2])[:,None] * fit['k'])+((fit['rho1'][:,2])[:,None] * fit['R'])])

coeffs = []
for i in xrange(ebv.shape[1]):
  coeffs.append( numpy.polyfit(ebv[0, i,:],ebv[1,i,:], 1))
coeffs = numpy.array(coeffs)

rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)
ebvav_s = numpy.percentile(ebv,(50,50-34,50+34),axis=1)


plt.errorbar(ebvav_s[0,0,:], ebvav_s[0,1,:], \
 xerr=(ebvav_s[0,0,]-ebvav_s[1,0,], ebvav_s[2,0,]-ebvav_s[0,0,]),\
 yerr=(ebvav_s[0,1,:]-ebvav_s[1,1,:], ebvav_s[2,1,:]-ebvav_s[0,1,:]),fmt='o',alpha=0.4,color='blue')
# plt.scatter(ebvav[0,::1000],ebvav[1,::1000],marker='.')
plt.ylabel(r'$A_{V}+ const $')
plt.xlabel(r'$E(B-V) + const$')
x = numpy.array([-0.15,0.45])
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend()
pp = PdfPages("output11/avebv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

rv=(fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R'])/ \
  ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']+(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R'])

ebvav_r = numpy.percentile(rv,(50,50-34,50+34),axis=0)
plt.errorbar(ebvav_s[0,0,:], ebvav_r[0,:], \
 xerr=(ebvav_s[0,0,:]-ebvav_s[1,0,:], ebvav_s[2,0,:]-ebvav_s[0,0,:]),\
 yerr=(ebvav_r[0,:]-ebvav_r[1,:],ebvav_r[2,:]-ebvav_r[0,:]),fmt='o',alpha=0.4)
# plt.scatter(ebvav[0,::1000],ebvav[1,::1000]/ebvav[0,::1000],marker='.')
plt.ylabel(r'$R_{V} $')
plt.xlabel(r'$E(B-V) + const$')
plt.ylim((-1,5))
pp = PdfPages("output11/rv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


ebv  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']
ebv = numpy.array([ebv,(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']])

ebvav=[]
for i in xrange(ebv.shape[1]):
  tmat_=[]
  for s in ['gamma','rho1']:
    ga=fit[s][i,:]
    y = numpy.array([numpy.dot(ga,dAdebv),numpy.dot(ga,dAdAv)])
    ans = numpy.linalg.solve(a,y)
    tmat_.append(ans)
  tmat_=numpy.array(tmat_)
  tmat_ = numpy.transpose(tmat_)
  inner = []
  for j in xrange(ebv.shape[2]):
    inner.append(numpy.dot(tmat_,ebv[:,i,j]))
  ebvav.append(inner)
ebvav = numpy.array(ebvav)

coeffs = []
for i in xrange(ebv.shape[1]):
  coeffs.append(   numpy.polyfit(ebvav[i,:,0],ebvav[i,:,1], 1))

coeffs = numpy.array(coeffs)

rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

print '$R^F_V={:6.2f}_{{-{:6.2f}}}^{{+ {:6.2f}}}  $'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
print '${:6.2f} -{:6.2f} + {:6.2f}$'.format(rbv[1],rbv[1]-mrbv[1],prbv[1]-rbv[1])

# ebv  = ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['R'])
# ebv = numpy.array([ebv,((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['k'])])

ebvav=[]
for ind in xrange(ebv.shape[2]):
  ebvav.append(numpy.dot(tmat,ebv[:,:,ind]))

ebvav=numpy.array(ebvav)

ebvav_s = numpy.percentile(ebvav,(50,50-34,50+34),axis=2)

# ebv  = numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R'],axis=0)
# ebv = numpy.array([ebv,numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k'],axis=0)])
# ebvav = numpy.dot(tmat,ebv)
# plt.scatter(ebvav[:,0,::200],ebvav[:,1,::200],marker='.',color='red')
plt.errorbar(ebvav_s[0,:,0], ebvav_s[0,:,1], \
 xerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),\
 yerr=(ebvav_s[0,:,1]-ebvav_s[1,:,1], ebvav_s[2,:,1]-ebvav_s[0,:,1]),fmt='o',alpha=0.4,color='blue')
# plt.scatter(ebvav[0,::1000],ebvav[1,::1000],marker='.')
plt.ylabel(r'$A^F_{V,eff}+ const $')
plt.xlabel(r'$E^F(B-V)_{eff} + const$')
x = numpy.array([-0.15,0.45])
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R^F_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend()
pp = PdfPages("output11/avebv_synth.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# import scipy.odr.odrpack as odrpack

# def f(B, x):
#     return B[0]*x + B[1]
# linear = odrpack.Model(f)
# # mydata = odrpack.Data(x, y, wd=1./np.power(sx,2), we=1./np.power(sy,2))
# mydata = odrpack.RealData(ebvav_s[0,:,0], ebvav_s[0,:,1], sx=(ebvav_s[2,:,0]-ebvav_s[1,:,0])/2, sy=(ebvav_s[2,:,1]-ebvav_s[1,:,1])/2)

# myodr = odrpack.ODR(mydata, linear, beta0=[3, 0.])
# myoutput = myodr.run()
# myoutput.pprint()


ebvav_r = numpy.percentile(ebvav[:,1,:]/ebvav[:,0,:],(50,50-34,50+34),axis=1)
plt.errorbar(ebvav_s[0,:,0], ebvav_r[0], \
 xerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),\
 yerr=(ebvav_r[0]-ebvav_r[1],ebvav_r[2]-ebvav_r[0]),fmt='o',alpha=0.4)
# plt.scatter(ebvav[0,::1000],ebvav[1,::1000]/ebvav[0,::1000],marker='.')
plt.ylabel(r'$R^F_{V,eff} $')
plt.xlabel(r'$E^F(B-V)_{eff} + const$')
plt.ylim((-1,5))
pp = PdfPages("output11/avrv_synth.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
# def lnprob(p, ebvav):
#   dum= (ebvav[1,:] - p[1]) - (p[0] * (ebvav[0,:] - p[2]))
#   ans = -0.5 * numpy.sum(dum**2)
#   print p,ans
#   return ans

# ndim, nwalkers = 3, 3*4
# zeros = numpy.zeros(ndim)
# zeros = numpy.array([2.23,-0.25,-0.1])

# pos = [zeros + 1e-4*numpy.random.randn(ndim) for i in range(nwalkers)]
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ebvav])
# sampler.run_mcmc(pos,1000)

# samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
# (y,ymin,ymax)  = numpy.percentile(samples[:,0],(50,50-34,50+34)) 
# print y, y-ymin, ymax-y
# (y,ymin,ymax)  = numpy.percentile(samples[:,1],(50,50-34,50+34)) 
# print y, y-ymin, ymax-y
# (y,ymin,ymax)  = numpy.percentile(samples[:,2],(50,50-34,50+34)) 
# print y, y-ymin, ymax-y