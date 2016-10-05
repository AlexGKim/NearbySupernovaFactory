#!/usr/bin/env python

import pickle
import numpy
import sncosmo
import scipy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Get the data
f = open('temp15.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

# Determine the plane approximaion for Fitzpatrick

def getFitzExt(efflam,av, ebv):
  f99 = sncosmo.F99Dust(r_v =av/ebv)
  f99.set(ebv=ebv)
  A_ = f99.propagate(efflam,1.)
  A1_ = -2.5*numpy.log10(A_)
  return A1_

efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])

# Partial derivatives with respect to av and ebv
av=0.1
ebv=0.1/2.27
A1=getFitzExt(efflam, av , ebv)

A2=getFitzExt(efflam, av+0.01 , ebv)
dAdAv = (A2 - A1)/0.01

A3=getFitzExt(efflam, av , ebv+0.001)
dAdebv = (A3 - A1)/0.001

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
for s in ['gamma','gamma1']:
  c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
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
print tmat
print numpy.linalg.norm(res,axis=1)/numpy.array(c_n)

# The matrix to transform the per-SN parameters from gamma to fitzpatrick
# A= gamma0 k0 + gamma1 k1 = ans00 F0 k0 + ans01 F1 k0 + ans10 F0 k1 + ans11 F1 k1 
#  = (ans00 k0 + ans10 k1)F0 + (ans01 k0 + ans11 k1)F1
tmat = numpy.transpose(tmat)


# Plot AV versus E(B-V) from the data

# container that contains E(B-V) and AV
ebv  = ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']+((fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None] * fit['k1']))
ebv = numpy.array([ebv,((fit['gamma'][:,2])[:,None] * fit['k'])+((fit['gamma1'][:,2])[:,None] * fit['k1'])])

# RV is calculated as a Monte Carlo, i.e RV is calculated for each link
coeffs = []
for i in xrange(ebv.shape[1]):
  coeffs.append(numpy.polyfit(ebv[0, i,:],ebv[1,i,:], 1))
coeffs = numpy.array(coeffs)

# the fit RV
rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

# the fit EBV and AV
ebvav_s = numpy.percentile(ebv,(50,50-34,50+34),axis=1)

plt.errorbar(ebvav_s[0,0,:], ebvav_s[0,1,:], \
 xerr=(ebvav_s[0,0,:]-ebvav_s[1,0,:], ebvav_s[2,0,:]-ebvav_s[0,0,:]),\
 yerr=(ebvav_s[0,1,:]-ebvav_s[1,1,:], ebvav_s[2,1,:]-ebvav_s[0,1,:]),fmt='o',alpha=0.4,color='blue')
plt.ylabel(r'$A_{V}+ const $')
plt.xlabel(r'$E(B-V) + const$')
x = numpy.array([-0.15,0.45])
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend()
pp = PdfPages("output15/avebv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# Calculation of native RV
rv=(fit['gamma'][:,2][:,None]*fit['k'] + fit['gamma1'][:,2][:,None]*fit['k1'])/ \
  ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']+(fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None] * fit['k1'])

ebvav_r = numpy.percentile(rv,(50,50-34,50+34),axis=0)
plt.errorbar(ebvav_s[0,0,:], ebvav_r[0,:], \
 xerr=(ebvav_s[0,0,:]-ebvav_s[1,0,:], ebvav_s[2,0,:]-ebvav_s[0,0,:]),\
 yerr=(ebvav_r[0,:]-ebvav_r[1,:],ebvav_r[2,:]-ebvav_r[0,:]),fmt='o',alpha=0.4)
plt.ylabel(r'$R_{V} $')
plt.xlabel(r'$E(B-V) + const$')
plt.ylim((-1,5))
pp = PdfPages("output15/rv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# Calculation of the weighted mean of RV for extreme blue and red samples
w = ebvav_s[0,0,:] < -0.05
err = (ebvav_r[2,:]-ebvav_r[1,:])/2
dum = ebvav_r[0,w]/err[w]**2
dum2= 1/err[w]**2
print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))
w = ebvav_s[0,0,:] > 0.1
err = (ebvav_r[2,:]-ebvav_r[1,:])/2
dum = ebvav_r[0,w]/err[w]**2
dum2= 1/err[w]**2
print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))

# Transform native parameters onto the Fitzpatrick plane

# The native E(B-V) parameters
ebv  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']
ebv = numpy.array([ebv,(fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None] * fit['k1']])

# For each link recalculate the transformation matrix and get the Fitzpatrick values
ebvav=[]
for i in xrange(ebv.shape[1]):
  tmat_=[]
  for s in ['gamma','gamma1']:
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

# For each link calculate the slope
coeffs = []
for i in xrange(ebv.shape[1]):
  coeffs.append(numpy.polyfit(ebvav[i,:,0],ebvav[i,:,1], 1))

coeffs = numpy.array(coeffs)

# the monte carlo regions of rv
rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

print '$R^F_V={:6.2f}_{{-{:6.2f}}}^{{+ {:6.2f}}}  $'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
print '${:6.2f} -{:6.2f} + {:6.2f}$'.format(rbv[1],rbv[1]-mrbv[1],prbv[1]-rbv[1])

# Plot syntehetic Fitzpatrix E(B-V) and AV with slope derived above

ebvav=[]
for ind in xrange(ebv.shape[2]):
  ebvav.append(numpy.dot(tmat,ebv[:,:,ind]))

ebvav=numpy.array(ebvav)

ebvav_s = numpy.percentile(ebvav,(50,50-34,50+34),axis=2)

plt.errorbar(ebvav_s[0,:,0], ebvav_s[0,:,1], \
 xerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),\
 yerr=(ebvav_s[0,:,1]-ebvav_s[1,:,1], ebvav_s[2,:,1]-ebvav_s[0,:,1]),fmt='o',alpha=0.4,color='blue')

plt.ylabel(r'$A^F_{V,eff}+ const $')
plt.xlabel(r'$E^F(B-V)_{eff} + const$')
x = numpy.array([-0.15,0.45])
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R^F_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend()
pp = PdfPages("output15/avebv_synth.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# plot Rv versus Av for the best fit

ebvav_r = numpy.percentile(ebvav[:,1,:]/ebvav[:,0,:],(50,50-34,50+34),axis=1)
plt.errorbar(ebvav_s[0,:,0], ebvav_r[0], \
 xerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),\
 yerr=(ebvav_r[0]-ebvav_r[1],ebvav_r[2]-ebvav_r[0]),fmt='o',alpha=0.4)
# plt.scatter(ebvav[0,::1000],ebvav[1,::1000]/ebvav[0,::1000],marker='.')
plt.ylabel(r'$R^F_{V,eff} $')
plt.xlabel(r'$E^F(B-V)_{eff} + const$')
plt.ylim((-1,5))
pp = PdfPages("output15/avrv_synth.pdf")
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
