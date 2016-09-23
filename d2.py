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

# f99 = sncosmo.F99Dust(r_v =2.97)
# f99.set(ebv=1.)
# A_ = f99.propagate(efflam,1.)
# A1 = -2.5*numpy.log10(A_)
# A1 = A1/A1[2]



# dum = fit['gamma'][:,None,:] + fit['R'][:,:,None]/fit['k'][:,:,None] * fit['rho1'][:,None,:]
# dum = dum/dum[:,:,2][:,:,None]
# B1=[]
# for i in xrange(5):
#   B1.append(numpy.percentile(dum[:,:,i],(50,50-34,50+34)))

# B1=numpy.array(B1)
# plt.errorbar(efflam,B1[:,0],yerr=(B1[:,0]-B1[:,1],B1[:,2]-B1[:,0]),fmt='o',color='blue',label='Model')
# plt.scatter(efflam,A1,color='red',label='Fitzpatrick')
# plt.legend()
# pp = PdfPages("output11/rvvector.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# x, xmin, xmax = numpy.percentile((dum[:,:,1]-1),(50,50-34,50+34),axis=0)
# y, ymin, ymax = numpy.percentile(1./(dum[:,:,1]-1),(50,50-34,50+34),axis=0)
# # plt.hist(numpy.median(1./(dum[:,:,1]-1),axis=0))
# plt.errorbar(x,y,xerr=(x-xmin,xmax-x),yerr=(y-ymin,ymax-y),fmt='.',alpha=0.4)
# plt.xlim((-1,2))
# plt.ylim((-7,10))
# plt.xlabel(r'$E(B-V)$')
# plt.ylabel(r'$R_V$')
# pp = PdfPages("output11/embrv.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# x, xmin, xmax = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],(50,50-34,50+34),axis=0)
# y, ymin, ymax = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=(x-xmin,xmax-x),yerr=(y-ymin,ymax-y),fmt='.',alpha=0.4)
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$E_\delta(B-V)$')
# pp = PdfPages("output11/egammaedelta.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# y, ymin, ymax = numpy.percentile(((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'])/((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']),(50,50-34,50+34),axis=0)

# plt.errorbar(x,y,xerr=(x-xmin,xmax-x),yerr=(y-ymin,ymax-y),fmt='.',alpha=0.4)
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$E_\delta(B-V)$')
# pp = PdfPages("output11/egammaedeltaratio.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


av=0.1
ebv=0.1/3.1
A1=getFitzExt(efflam, av , ebv)

A2=getFitzExt(efflam, av+0.01 , ebv)
dAdAv = (A2 - A1)/0.01

A3=getFitzExt(efflam, av , ebv+0.001)
dAdebv = (A3 - A1)/0.001


dAdebv_norm = numpy.linalg.norm(dAdebv)
dAdAv_norm = numpy.linalg.norm(dAdAv)

# dAdebv = dAdebv/dAdebv_norm
# dAdAv = dAdAv/dAdAv_norm

tmat = []
res = []
c_n = []
cs = []
for s in ['rho1','gamma']:
  c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
  cs.append(c)
  c_norm = numpy.linalg.norm(c)
  c_n.append(c_norm)
  # c=c/c_norm

  y = numpy.array([numpy.dot(c,dAdebv),numpy.dot(c,dAdAv)])
  norm_dAdebv = numpy.dot(dAdebv, dAdebv)
  norm_dAdAv = numpy.dot(dAdAv, dAdAv)
  cross = numpy.dot(dAdebv, dAdAv)

  a = numpy.array([[norm_dAdebv,cross],[cross,norm_dAdAv]])
  ans = numpy.linalg.solve(a,y)

  tmat.append(ans)
  ans = c-ans[0]*dAdebv - ans[1]*dAdAv
  res.append(ans)

tmat = numpy.array(tmat)
res= numpy.array(res)

print tmat
tmat = numpy.transpose(tmat)
print numpy.linalg.norm(res,axis=1)/numpy.array(c_n)

ebv  = ((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']).flatten()
ebv = numpy.array([ebv,((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']).flatten()])
ebvav = numpy.dot(tmat,ebv)

coeffs, cov = numpy.polyfit(ebvav[0,:],ebvav[1,:], 1,cov=True)

# ebv  = numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R'],axis=0)
# ebv = numpy.array([ebv,numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k'],axis=0)])
# ebvav = numpy.dot(tmat,ebv)
plt.scatter(ebvav[0,::1000],ebvav[1,::1000],marker='.')
plt.xlabel(r'$A_{V,eff}+ const $')
plt.ylabel(r'$E(B-V)_{eff} + const$')
x = numpy.array([-0.1,0.4])
plt.plot(x, coeffs[1]+coeffs[0]*x,label=r'$R_V={:6.2f}$'.format(coeffs[0]))
plt.legend()
pp = PdfPages("output11/avrv_synth.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
#plt.scatter(ebvav[1,:],(ebvav[1,:]+0.3)/(ebvav[0,:]+2))



# def lnprob(p, ebvav):
#   dum= ebvav[1] - (p[0]*ebvav[0] + p[1])
#   ans = -0.5 * numpy.sum(dum**2)
#   return ans

# ndim, nwalkers = 2, 3*4
# zeros = numpy.zeros(ndim)
# zeros = numpy.array([2.23,0.])

# pos = [zeros + 1e-4*numpy.random.randn(ndim) for i in range(nwalkers)]
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ebvav])
# sampler.run_mcmc(pos,1000)

# samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
# (y,ymin,ymax)  = numpy.percentile(samples[:,0],(50,50-34,50+34)) 
# print y, y-ymin, ymax-y
# (y,ymin,ymax)  = numpy.percentile(samples[:,1],(50,50-34,50+34)) 
# print y, y-ymin, ymax-y

wef