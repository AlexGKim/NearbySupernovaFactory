#!/usr/bin/env python
import pickle
import numpy
import sncosmo
import scipy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import f99_band
import emcee
import matplotlib as mpl

mpl.rcParams['font.size'] = 14

# Get the data
f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

# Determine the plane approximaion for Fitzpatrick

# Partial derivatives with respect to av and ebv
av=0.1
ebv=0.1/2.5
A1= f99_band.A_X(r_v=av/ebv, ebv=ebv)

# pkl_file = open('fitz.pkl', 'r')
# a=pickle.load(pkl_file)
# pkl_file.close()

# AX = a[0]* av + a[1] * av**2 \
#   + a[2]* ebv+ a[3] * ebv**2 \
#   + a[4] * av* ebv \
#   + a[5]* av**3 \
#   + a[6] * ebv**3 \
#   + a[7] * (av**2) * ebv \
#   + a[8] * av * (ebv**2)

# plt.plot(A1-AX)

# plt.show()

# wefe
A2= f99_band.A_X(r_v=(av+0.01)/ebv, ebv=ebv)
dAdAv = (A2 - A1)/0.01

A3= f99_band.A_X(r_v=av/(ebv+0.001), ebv=ebv+0.001)
dAdebv = (A3 - A1)/0.001

print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdAv)
print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdebv)

# av=0.1
# ebv=0.1/3.1
# A1= f99_band.A_X(r_v=av/ebv, ebv=ebv)

# A2= f99_band.A_X(r_v=(av+0.01)/ebv, ebv=ebv)
# dAdAv = (A2 - A1)/0.01

# A3= f99_band.A_X(r_v=av/(ebv+0.001), ebv=ebv+0.001)
# dAdebv = (A3 - A1)/0.001

# print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdAv)
# print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdebv)

# av=0.1
# ebv=0.1/1.9
# A1= f99_band.A_X(r_v=av/ebv, ebv=ebv)

# A2= f99_band.A_X(r_v=(av+0.01)/ebv, ebv=ebv)
# dAdAv = (A2 - A1)/0.01

# A3= f99_band.A_X(r_v=av/(ebv+0.001), ebv=ebv+0.001)
# dAdebv = (A3 - A1)/0.001

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
  print "{:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}".format(c[0],c[1],c[2],c[3],c[4])
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
print res


# The matrix to transform the per-SN parameters from gamma to fitzpatrick
# A= gamma0 k0 + gamma1 k1 = ans00 F0 k0 + ans01 F1 k0 + ans10 F0 k1 + ans11 F1 k1 
#  = (ans00 k0 + ans10 k1)F0 + (ans01 k0 + ans11 k1)F1
tmat = numpy.transpose(tmat)


# Plot vectors in UVI
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
fig = plt.figure()
ax = fig.gca(projection='3d')
dum  = numpy.sqrt(cs[0][0]**2+cs[0][2]**2+cs[0][4]**2)
ax.plot([0,cs[0][0]/dum],[0,cs[0][2]/dum],[0,cs[0][4]/dum],label=r'$\gamma^0_X/(\gamma^0_B-\gamma^0_V)$')
dum  = numpy.sqrt(cs[1][0]**2+cs[1][2]**2+cs[1][4]**2)
ax.plot([0,cs[1][0]/dum],[0,cs[1][2]/dum],[0,cs[1][4]/dum],label=r'$\gamma^1_X/(\gamma^1_B-\gamma^1_V)$')
dum  = numpy.sqrt(dAdAv[0]**2+dAdAv[2]**2+dAdAv[4]**2)
ax.plot([0,dAdAv[0]/dum],[0,dAdAv[2]/dum],[0,dAdAv[4]/dum],label=r'$a(X)$',ls='--')
dum  = numpy.sqrt(dAdebv[0]**2+dAdebv[2]**2+dAdebv[4]**2)
ax.plot([0,dAdebv[0]/dum],[0,dAdebv[2]/dum],[0,dAdebv[4]/dum],label=r'$b(X)$',ls='--')
crap = dAdAv + dAdebv/2.4
dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/2.4$',ls=':')
crap = -6.8*dAdAv + dAdebv
dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$-6.8a(X)+b(X)$',ls=':',color='black')
crap = dAdAv + dAdebv/2.6
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/2.6$',ls=':',color='black')
#crap = dAdebv
#ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$b(X)$',ls=':')
ax.legend(prop={'size':14})
ax.set_xlabel(r'$U$',labelpad=18)
ax.set_ylabel(r'$V$',labelpad=18)
ax.set_zlabel(r'$I$',labelpad=18)
ax.xaxis.set_ticks(numpy.arange(-.5,1.1,.25))
ax.yaxis.set_ticks(numpy.arange(-.8,.81,.4))
ax.view_init(elev=2, azim=-114)
pp = PdfPages("output11/plane0.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output11/plane1.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# Plot vectors in BVR
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
fig = plt.figure()
ax = fig.gca(projection='3d')
dum  = numpy.sqrt(cs[0][1]**2+cs[0][2]**2+cs[0][3]**2)
ax.plot([0,cs[0][1]/dum],[0,cs[0][2]/dum],[0,cs[0][3]/dum],label=r'$\gamma^0_X/(\gamma^0_B-\gamma^0_V)$')
dum  = numpy.sqrt(cs[1][1]**2+cs[1][2]**2+cs[1][3]**2)
ax.plot([0,cs[1][1]/dum],[0,cs[1][2]/dum],[0,cs[1][3]/dum],label=r'$\gamma^1_X/(\gamma^1_B-\gamma^1_V)$')
dum  = numpy.sqrt(dAdAv[1]**2+dAdAv[2]**2+dAdAv[3]**2)
ax.plot([0,dAdAv[1]/dum],[0,dAdAv[2]/dum],[0,dAdAv[3]/dum],label=r'$a(X)$',ls='--')
dum  = numpy.sqrt(dAdebv[1]**2+dAdebv[2]**2+dAdebv[3]**2)
ax.plot([0,dAdebv[1]/dum],[0,dAdebv[2]/dum],[0,dAdebv[3]/dum],label=r'$b(X)$',ls='--')
# crap = dAdAv + dAdebv/2.4
# dum  = numpy.sqrt(crap[1]**2+crap[2]**2+crap[3]**2)
# ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$a(X)+b(X)/2.4$',ls=':')
crap = dAdAv + dAdebv/2.6
dum  = numpy.sqrt(crap[1]**2+crap[2]**2+crap[3]**2)
ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$a(X)+b(X)/2.6$',ls=':')
crap = -6.1*dAdAv + dAdebv
dum  = numpy.sqrt(crap[1]**2+crap[2]**2+crap[3]**2)
ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$-6.1a(X)+b(X)$',ls=':',color='black')
#crap = dAdebv
#ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$b(X)$',ls=':')
ax.legend(prop={'size':14})
ax.set_xlabel(r'$B$',labelpad=18)
ax.set_ylabel(r'$V$',labelpad=18)
ax.set_zlabel(r'$R$',labelpad=18)
ax.xaxis.set_ticks(numpy.arange(-.5,1.1,.25))
ax.yaxis.set_ticks(numpy.arange(-.8,.81,.4))
ax.view_init(elev=2, azim=-114)
pp = PdfPages("output11/plane0BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output11/plane1BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wef

# Plot AV versus E(B-V) from the data

# container that contains E(B-V) and AV
ebv  = ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']+((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']))
ebv = numpy.array([ebv,((fit['gamma'][:,2])[:,None] * fit['k'])+((fit['rho1'][:,2])[:,None] * fit['R'])])

ebv_mn = numpy.mean(ebv,axis=1)

ebv_cov = numpy.zeros((ebv.shape[2],2,2))
for i in xrange(ebv.shape[2]):
  ebv_cov[i,:,:] = numpy.cov(ebv[:,:,i])

# ebv_icov = numpy.zeros((ebv.shape[2],2,2))
# for i in xrange(ebv.shape[2]):
#   ebv_icov[i,:,:] = numpy.linalg.inv(ebv_cov[i,:,:])


# def lnprob(pars, ebv_mn, ebv_icov):
#   ans = 0.
#   for i in xrange(ebv_mn.shape[1]):
#     A1 = f99_band_bv.A_X(r_v=pars[0], ebv=pars[1+i])
#     vec = ebv_mn[:,i] - numpy.array([A1[0]-A1[1], A1[1]])
#     term  = numpy.dot(vec, numpy.dot(ebv_icov[i], vec))
#     ans -= 0.5 * term
#   return ans

# ndim, nwalkers = len(ebv_icov)+1, 2*(len(ebv_icov)+1)

# p0 = [numpy.random.normal(0,0.1,size=ndim) * numpy.concatenate(([1],numpy.zeros(ndim-1)+0.05)) \
#   +numpy.concatenate(([2.5],ebv_mn[0,:])) for i in range(nwalkers)]

# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ebv_mn, ebv_icov])
# sampler.run_mcmc(p0, 200)

# plt.plot(sampler.chain[1,:,0])
# samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

# wefwe

import pystan

data = {'D': len(ebv_cov),'obs': numpy.transpose(ebv_mn), 'obs_cov': ebv_cov}
# print len(ebv_cov), ebv_mn.shape, data['obs_cov'].shape

init1 = {'x0' : -0.01, 'x1': 2.5, 'ebv':ebv_mn[0,:]}#

sm = pystan.StanModel(file='rv.stan')
fit1 = sm.sampling(data=data, iter=1000, chains=4,init=[init1,init1,init1,init1])
links = fit1.extract(('x1','x0'))
rbv, mrbv, prbv= numpy.percentile([links['x1'],links['x0']],(50,50-34,50+34),axis=1)
print "Fit RV"
print "${:6.2f}_{:6.2f}^{:6.2f}$".format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
# # RV is calculated as a Monte Carlo, i.e RV is calculated for each link
# coeffs = []
# for i in xrange(ebv.shape[1]):
#   coeffs.append(numpy.polyfit(ebv[0, i,:],ebv[1,i,:], 1))
# coeffs = numpy.array(coeffs)

# # the fit RV
# rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

# the fit EBV and AV
ebvav_s = numpy.percentile(ebv,(50,50-34,50+34),axis=1)

plt.errorbar(ebvav_s[0,0,:], ebvav_s[0,1,:], \
 xerr=(ebvav_s[0,0,:]-ebvav_s[1,0,:], ebvav_s[2,0,:]-ebvav_s[0,0,:]),\
 yerr=(ebvav_s[0,1,:]-ebvav_s[1,1,:], ebvav_s[2,1,:]-ebvav_s[0,1,:]),fmt='o',alpha=0.4,color='orange')
plt.ylabel(r'$A_{V}$')
plt.xlabel(r'$E(B-V)$')
x = numpy.array([-0.15,0.45])
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
rbvs = [1.1,1.9,2.5,3.1]
for rbv_ in rbvs:
  x=[]
  y=[]
  for ebv in numpy.arange(-0.15,0.5,0.02):
    A1= f99_band.A_X(r_v=rbv_, ebv=ebv)
    x.append(A1[1]-A1[2])
    y.append(A1[2])
  plt.plot(x,y,label=r'$R^F={:6.1f}$'.format(rbv_))
plt.legend(loc=4)
plt.xlim((-0.1,0.4))
plt.ylim((-0.7,1.2))
pp = PdfPages("output11/avebv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# Calculation of native RV
ebv0=-0.080
av0=-0.155


ebv2=((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * (fit['k'])+(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * (fit['R']))
ebv2 = ebv2-ebv0

av=(fit['gamma'][:,2][:,None]*(fit['k']) + fit['rho1'][:,2][:,None]*(fit['R']))

av = av-av0

rv=av/ebv2

ebvav_r = numpy.percentile(rv,(50,50-34,50+34),axis=0)
ebvav2_s = numpy.percentile(ebv2,(50,50-34,50+34),axis=0)

plt.errorbar(ebvav2_s[0,:], ebvav_r[0,:], \
 xerr=(ebvav2_s[0,:]-ebvav2_s[1,:], ebvav2_s[2,:]-ebvav2_s[0,:]),\
 yerr=(ebvav_r[0,:]-ebvav_r[1,:],ebvav_r[2,:]-ebvav_r[0,:]),fmt='o',alpha=0.3)

xerr=(ebvav2_s[2,:]+ebvav2_s[1,:])/2
yerr=(ebvav_r[2,:]+ebvav_r[1,:])/2
so = numpy.argsort(ebvav2_s[0,:])
so_split = numpy.array_split(so,15)
bx=[]
dbx=[]
by=[]
dby=[]

# for uso in so_split:
#   bx.append(numpy.sum(ebvav2_s[0,uso]/xerr[uso]**2)/numpy.sum(1./xerr[uso]**2))
#   # dbx.append(1./numpy.sqrt(numpy.sum(1./xerr[uso]**2)))
#   # by.append(numpy.sum(ebvav_r[0,uso]/yerr[uso]**2)/numpy.sum(1./yerr[uso]**2))
#   # dby.append(1./numpy.sqrt(numpy.sum(1./yerr[uso]**2)))  
#   by.append(numpy.mean(ebvav_r[0,uso]))
#   dby.append(numpy.std(ebvav_r[0,uso]))

# plt.errorbar(bx,by,yerr=[dby,dby], fmt='o', markersize='10',color='black',elinewidth=2)
for rbv in rbvs:
  x=[]
  y=[]
  for ebv in numpy.arange(-0.05,0.6,0.02):
    A1= f99_band.A_X(r_v=rbv, ebv=ebv)
    x.append(A1[1]-A1[2])
    y.append(A1[2]/(A1[1]-A1[2]))
  plt.plot(x,y,label=r'$R^F={:6.1f}$'.format(rbv))
plt.xlim((-0.02,0.45))
plt.ylabel(r'$R^T_{V} $')
plt.xlabel(r'$E^T(B-V)$')
plt.ylim((-2,5))
plt.legend(loc=4)
pp = PdfPages("output11/rv.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# Calculation of the weighted mean of RV for extreme blue and red samples
w = ebvav2_s[0,:] < 0.05
err = (ebvav_r[2,:]-ebvav_r[1,:])/2
dum = ebvav_r[0,w]/err[w]**2
dum2= 1/err[w]**2
print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))
w = ebvav2_s[0,:] > 0.1
err = (ebvav_r[2,:]-ebvav_r[1,:])/2
dum = ebvav_r[0,w]/err[w]**2
dum2= 1/err[w]**2
print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))

# Transform native parameters onto the Fitzpatrick plane

# The native E(B-V) parameters
ebv  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']
ebv = numpy.array([ebv,(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']])

# For each link recalculate the transformation matrix and get the Fitzpatrick values
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

# For each link calculate the slope
coeffs = []
for i in xrange(ebv.shape[1]):
  coeffs.append(numpy.polyfit(ebvav[i,:,0],ebvav[i,:,1], 1))

coeffs = numpy.array(coeffs)

# the monte carlo regions of rv
rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

print '$R^F={:6.2f}_{{-{:6.2f}}}^{{+ {:6.2f}}}  $'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
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
plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R^F={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend()
pp = PdfPages("output11/avebv_synth.pdf")
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
pp = PdfPages("output11/avrv_synth.pdf")
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
