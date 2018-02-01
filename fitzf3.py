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

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

mpl.rcParams['font.size'] = 14

# Get the data
f = open('fix3_decorr.pkl','rb')
fit = pickle.load(f)
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


#residual
AX = f99_band.A_X(r_v=2.0, ebv=1./2)
print AX -A1 - dAdAv*(1.-av) - dAdebv*(1./2 - ebv)
AX = f99_band.A_X(r_v=3.5, ebv=1./3.5)
print AX -A1 - dAdAv*(1.-av) - dAdebv*(1./3.5 - ebv)

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

a = numpy.array([[norm_dAdAv,cross],[cross,norm_dAdebv]])
a_save = numpy.array(a)


a = numpy.array([[norm_dAdAv,cross],[cross,norm_dAdebv]])
a_save = numpy.array(a)

# calculate the projection onto the plane
proj=[]
kappa=[]
for i in xrange(fit['gamma'].shape[0]):
  proj2=[]
  kappa2=[]
  for s in ['gamma','rho1']:
    y = numpy.array([numpy.dot(fit[s][i,:],dAdAv), numpy.dot(fit[s][i,:],dAdebv)])
    ans = numpy.linalg.solve(a,y)
    kappa2.append(ans)
    proj2.append((numpy.linalg.norm(ans[1]*dAdebv + ans[0]*dAdAv)/numpy.linalg.norm(fit[s][i,:]))**2)
  kappa.append(kappa2)
  proj.append(proj2)

kappa = numpy.array(kappa)

print "R1 (1/kappa1) term"
dum1, dumm, dump =  numpy.percentile(kappa[:,0,0]/kappa[:,0,1],(50,50-34,50+34))
print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1,dump-dum1,dumm-dum1)
print "1/kappa2 term"
dum1, dumm, dump =  numpy.percentile(kappa[:,1,1]/kappa[:,1,0],(50,50-34,50+34))
print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1,dump-dum1,dumm-dum1)


print "projection of gamma0, gamma1 onto Fitzpatrick plane"
proj=numpy.array(proj)
dum1, dumm, dump =  numpy.percentile(proj,(50,50-34,50+34),axis=0)
print "{:.4f}^{{+{:.4f}}}_{{{:.4f}}}".format(dum1[0],dump[0]-dum1[0],dumm[0]-dum1[0])
print "{:.4f}^{{+{:.4f}}}_{{{:.4f}}}".format(dum1[1],dump[1]-dum1[1],dumm[1]-dum1[1])

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

  # y = numpy.array([numpy.dot(c,dAdebv),numpy.dot(c,dAdAv)])
  y = numpy.array([numpy.dot(c,dAdAv), numpy.dot(c,dAdebv)])
  ans = numpy.linalg.solve(a,y)

  tmat.append(ans)
  ans = c-ans[1]*dAdebv - ans[0]*dAdAv
  res.append(ans)



tmat = numpy.array(tmat)
res= numpy.array(res)

#print the matrix and the residues
print tmat
print (numpy.linalg.norm(res,axis=1)/numpy.array(c_n))**2
print res


kappa1 = tmat[0,1]/tmat[0,0]
kappa2 = tmat[1,0]/tmat[1,1]


# The matrix to transform the per-SN parameters from gamma to fitzpatrick
# A= gamma0 k0 + gamma1 k1 = ans00 F0 k0 + ans01 F1 k0 + ans10 F0 k1 + ans11 F1 k1 
#  = (ans00 k0 + ans10 k1)F0 + (ans01 k0 + ans11 k1)F1
tmat = numpy.transpose(tmat)



#range of R^F 
r1 = []
r2 = []
for ind in xrange(fit['gamma'].shape[0]):
  tmat = []
  cs = []
  for s in ['gamma','rho1']:
    c = fit[s][ind,:]/((fit[s][ind,1]-fit[s][ind,2]))
    cs.append(c)

    y = numpy.array([numpy.dot(c,dAdAv),numpy.dot(c,dAdebv)])

    ans = numpy.linalg.solve(a,y)

    tmat.append(ans)

  tmat = numpy.array(tmat)
  r1.append(tmat[0,1]/tmat[0,0])
  
r1 = numpy.array(r1)
print numpy.percentile(1/r1,(50,50-34,50+34))
# plt.hist(r1)
# plt.show()



# Plot vectors in UVI
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
fig = plt.figure()
ax = fig.gca(projection='3d')
dum  = numpy.sqrt(cs[0][0]**2+cs[0][2]**2+cs[0][4]**2)
ax.plot([0,cs[0][0]/dum],[0,cs[0][2]/dum],[0,cs[0][4]/dum],label=r'$\hat{\gamma}^0$',color='blue')
a = Arrow3D([0,cs[0][0]/dum],[0,cs[0][2]/dum],[0,cs[0][4]/dum], mutation_scale=10, 
            arrowstyle="-|>",color='blue')
ax.add_artist(a)
dum  = numpy.sqrt(cs[1][0]**2+cs[1][2]**2+cs[1][4]**2)
ax.plot([0,cs[1][0]/dum],[0,cs[1][2]/dum],[0,cs[1][4]/dum],label=r'$\hat{\gamma}^1$',color='green')
a = Arrow3D([0,cs[1][0]/dum],[0,cs[1][2]/dum],[0,cs[1][4]/dum], mutation_scale=10, 
            arrowstyle="-|>",color='green')
ax.add_artist(a)

dum  = numpy.sqrt(dAdAv[0]**2+dAdAv[2]**2+dAdAv[4]**2)
ax.plot([0,dAdAv[0]/dum],[0,dAdAv[2]/dum],[0,dAdAv[4]/dum],label=r'$\hat{a(X)}$',ls='--',color='blue')
a = Arrow3D([0,dAdAv[0]/dum],[0,dAdAv[2]/dum],[0,dAdAv[4]/dum],ls='--', mutation_scale=10, color='blue',
            arrowstyle="-|>")
ax.add_artist(a)
dum  = numpy.sqrt(dAdebv[0]**2+dAdebv[2]**2+dAdebv[4]**2)
ax.plot([0,dAdebv[0]/dum],[0,dAdebv[2]/dum],[0,dAdebv[4]/dum],label=r'$\hat{b(X)}$',ls='--', color='green')
a = Arrow3D([0,dAdebv[0]/dum],[0,dAdebv[2]/dum],[0,dAdebv[4]/dum],ls='--', mutation_scale=10, color='green',
            arrowstyle="-|>")
ax.add_artist(a)
crap = dAdAv + dAdebv*kappa1
dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/{:4.2f}$'.format(1/kappa1),ls=':',color='black')
a = Arrow3D([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum], mutation_scale=10, ls=':',color='black', linewidth=3.,
            arrowstyle="-|>",alpha=0.8)
ax.add_artist(a)


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# dum  = numpy.sqrt(cs[0][0]**2+cs[0][2]**2+cs[0][4]**2)
# ax.plot([0,cs[0][0]/dum],[0,cs[0][2]/dum],[0,cs[0][4]/dum],label=r'$\gamma^0_X/(\gamma^0_B-\gamma^0_V)$')
# dum  = numpy.sqrt(cs[1][0]**2+cs[1][2]**2+cs[1][4]**2)
# ax.plot([0,cs[1][0]/dum],[0,cs[1][2]/dum],[0,cs[1][4]/dum],label=r'$\gamma^1_X/(\gamma^1_B-\gamma^1_V)$')
# dum  = numpy.sqrt(dAdAv[0]**2+dAdAv[2]**2+dAdAv[4]**2)
# ax.plot([0,dAdAv[0]/dum],[0,dAdAv[2]/dum],[0,dAdAv[4]/dum],label=r'$a(X)$',ls='--')
# dum  = numpy.sqrt(dAdebv[0]**2+dAdebv[2]**2+dAdebv[4]**2)
# ax.plot([0,dAdebv[0]/dum],[0,dAdebv[2]/dum],[0,dAdebv[4]/dum],label=r'$b(X)$',ls='--')
# crap = dAdAv + dAdebv*kappa1
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/{:4.2f}$'.format(1/kappa1),ls=':')
# crap = kappa2*dAdAv + dAdebv
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'${:4.2f}a(X)+b(X)$'.format(kappa2),ls=':',color='black')

# crap = -16*dAdAv - dAdebv
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$-6.8a(X)-b(X)$',ls=':',color='black')

# crap = dAdAv + dAdebv/2.6
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/2.6$',ls=':',color='black')
#crap = dAdebv
#ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$b(X)$',ls=':')
ax.legend(prop={'size':14})
ax.set_xlabel(r'$\hat{U}$',labelpad=18)
ax.set_ylabel(r'$\hat{V}$',labelpad=18)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$\hat{I}$',labelpad=9,rotation=0)
ax.xaxis.set_ticks(numpy.arange(-.5,1.1,.25))
ax.yaxis.set_ticks(numpy.arange(-.8,.81,.4))
ax.view_init(elev=2, azim=-114)
pp = PdfPages("output_fix3/plane0.pdf")
# plt.tight_layout()
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output_fix3/plane1.pdf")
# plt.tight_layout()
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()

# Plot vectors in BVR
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
crap = dAdAv + dAdebv*kappa1
dum  = numpy.sqrt(crap[1]**2+crap[2]**2+crap[3]**2)
ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$a(X)+b(X)/{:4.2f}$'.format(1/kappa1),ls=':')
crap = kappa2*dAdAv + dAdebv
dum  = numpy.sqrt(crap[1]**2+crap[2]**2+crap[3]**2)
ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'${:4.2f}a(X)+b(X)$'.format(kappa2),ls=':',color='black')
#crap = dAdebv
#ax.plot([0,crap[1]/dum],[0,crap[2]/dum],[0,crap[3]/dum],label=r'$b(X)$',ls=':')
ax.legend(prop={'size':14})
ax.set_xlabel(r'$B$',labelpad=18)
ax.set_ylabel(r'$V$',labelpad=18)
ax.set_zlabel(r'$R$',labelpad=18)
ax.xaxis.set_ticks(numpy.arange(-.5,1.1,.25))
ax.yaxis.set_ticks(numpy.arange(-.8,.81,.4))
ax.view_init(elev=2, azim=-114)
pp = PdfPages("output_fix3/plane0BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output_fix3/plane1BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()





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

# import pystan

# data = {'D': len(ebv_cov),'obs': numpy.transpose(ebv_mn), 'obs_cov': ebv_cov}
# # print len(ebv_cov), ebv_mn.shape, data['obs_cov'].shape

# init1 = {'x0' : -0.01, 'x1': 2.5, 'ebv':ebv_mn[0,:]}#

# sm = pystan.StanModel(file='rv.stan')
# fit1 = sm.sampling(data=data, iter=1000, chains=4,init=[init1,init1,init1,init1])
# links = fit1.extract(('x1','x0'))
# rbv, mrbv, prbv= numpy.percentile([links['x1'],links['x0']],(50,50-34,50+34),axis=1)
# print "Fit RV"
# print "${:6.2f}_{:6.2f}^{:6.2f}$".format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
# # # RV is calculated as a Monte Carlo, i.e RV is calculated for each link
# # coeffs = []
# # for i in xrange(ebv.shape[1]):
# #   coeffs.append(numpy.polyfit(ebv[0, i,:],ebv[1,i,:], 1))
# # coeffs = numpy.array(coeffs)

# # # the fit RV
# # rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

# # the fit EBV and AV
# ebvav_s = numpy.percentile(ebv,(50,50-34,50+34),axis=1)

# plt.errorbar(ebvav_s[0,0,:], ebvav_s[0,1,:], \
#  xerr=(ebvav_s[0,0,:]-ebvav_s[1,0,:], ebvav_s[2,0,:]-ebvav_s[0,0,:]),\
#  yerr=(ebvav_s[0,1,:]-ebvav_s[1,1,:], ebvav_s[2,1,:]-ebvav_s[0,1,:]),fmt='o',alpha=0.4,color='orange')
# plt.ylabel(r'$A_{V}$')
# plt.xlabel(r'$E(B-V)$')
# x = numpy.array([-0.15,0.45])
# plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R_V={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
# rbvs = [1.1,1.9,2.5,3.1]
# for rbv_ in rbvs:
#   x=[]
#   y=[]
#   for ebv in numpy.arange(-0.15,0.5,0.02):
#     A1= f99_band.A_X(r_v=rbv_, ebv=ebv)
#     x.append(A1[1]-A1[2])
#     y.append(A1[2])
#   plt.plot(x,y,label=r'$R^F={:6.1f}$'.format(rbv_))
# plt.legend(loc=4)
# plt.xlim((-0.1,0.4))
# plt.ylim((-0.7,1.2))
# pp = PdfPages("output_fix3/avebv.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# # Calculation of native RV
# ebv0=-0.080
# av0=-0.155


# ebv2=((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * (fit['k'])+(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * (fit['R']))
# ebv2 = ebv2-ebv0

# av=(fit['gamma'][:,2][:,None]*(fit['k']) + fit['rho1'][:,2][:,None]*(fit['R']))

# av = av-av0

# rv=av/ebv2

# ebvav_r = numpy.percentile(rv,(50,50-34,50+34),axis=0)
# ebvav2_s = numpy.percentile(ebv2,(50,50-34,50+34),axis=0)

# plt.errorbar(ebvav2_s[0,:], ebvav_r[0,:], \
#  xerr=(ebvav2_s[0,:]-ebvav2_s[1,:], ebvav2_s[2,:]-ebvav2_s[0,:]),\
#  yerr=(ebvav_r[0,:]-ebvav_r[1,:],ebvav_r[2,:]-ebvav_r[0,:]),fmt='o',alpha=0.3)

# xerr=(ebvav2_s[2,:]+ebvav2_s[1,:])/2
# yerr=(ebvav_r[2,:]+ebvav_r[1,:])/2
# so = numpy.argsort(ebvav2_s[0,:])
# so_split = numpy.array_split(so,15)
# bx=[]
# dbx=[]
# by=[]
# dby=[]

# # for uso in so_split:
# #   bx.append(numpy.sum(ebvav2_s[0,uso]/xerr[uso]**2)/numpy.sum(1./xerr[uso]**2))
# #   # dbx.append(1./numpy.sqrt(numpy.sum(1./xerr[uso]**2)))
# #   # by.append(numpy.sum(ebvav_r[0,uso]/yerr[uso]**2)/numpy.sum(1./yerr[uso]**2))
# #   # dby.append(1./numpy.sqrt(numpy.sum(1./yerr[uso]**2)))  
# #   by.append(numpy.mean(ebvav_r[0,uso]))
# #   dby.append(numpy.std(ebvav_r[0,uso]))

# # plt.errorbar(bx,by,yerr=[dby,dby], fmt='o', markersize='10',color='black',elinewidth=2)
# for rbv in rbvs:
#   x=[]
#   y=[]
#   for ebv in numpy.arange(-0.05,0.6,0.02):
#     A1= f99_band.A_X(r_v=rbv, ebv=ebv)
#     x.append(A1[1]-A1[2])
#     y.append(A1[2]/(A1[1]-A1[2]))
#   plt.plot(x,y,label=r'$R^F={:6.1f}$'.format(rbv))
# plt.xlim((-0.02,0.45))
# plt.ylabel(r'$R^T_{V} $')
# plt.xlabel(r'$E^T(B-V)$')
# plt.ylim((-2,5))
# plt.legend(loc=4)
# pp = PdfPages("output_fix3/rv.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# # Calculation of the weighted mean of RV for extreme blue and red samples
# w = ebvav2_s[0,:] < 0.05
# err = (ebvav_r[2,:]-ebvav_r[1,:])/2
# dum = ebvav_r[0,w]/err[w]**2
# dum2= 1/err[w]**2
# print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))
# w = ebvav2_s[0,:] > 0.1
# err = (ebvav_r[2,:]-ebvav_r[1,:])/2
# dum = ebvav_r[0,w]/err[w]**2
# dum2= 1/err[w]**2
# print '${:6.2f} \pm {:6.2f}$'.format(dum.sum()/dum2.sum(),1./numpy.sqrt(dum2.sum()))

# # Transform native parameters onto the Fitzpatrick plane

# The native E(B-V) parameters
# ebv  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']
# ebv = numpy.array([ebv,(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']])
ebv  = fit['k']
ebv = numpy.array([ebv,fit['R']])
a=a_save
# For each link recalculate the transformation matrix and get the Fitzpatrick values
ebvav=[]
for i in xrange(ebv.shape[1]):
  tmat_=[]
  for s in ['gamma','rho1']:
    ga=fit[s][i,:]
    y = numpy.array([numpy.dot(ga,dAdAv),numpy.dot(ga,dAdebv)])
    ans = numpy.linalg.solve(a,y)
    tmat_.append(ans)
  tmat_=numpy.array(tmat_)
  tmat_ = numpy.transpose(tmat_)
  inner = []
  for j in xrange(ebv.shape[2]):
    inner.append(numpy.dot(tmat_,ebv[:,i,j]))
  ebvav.append(inner)
ebvav = numpy.array(ebvav)

# # For each link calculate the slope
# coeffs = []
# for i in xrange(ebv.shape[1]):
#   coeffs.append(numpy.polyfit(ebvav[i,:,0],ebvav[i,:,1], 1))

# coeffs = numpy.array(coeffs)

# # # the monte carlo regions of rv
# rbv, mrbv, prbv= numpy.percentile(coeffs,(50,50-34,50+34),axis=0)

# print '$R^F={:6.2f}_{{-{:6.2f}}}^{{+ {:6.2f}}}  $'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0])
# print '${:6.2f} -{:6.2f} + {:6.2f}$'.format(rbv[1],rbv[1]-mrbv[1],prbv[1]-rbv[1])

# # Plot syntehetic Fitzpatrix E(B-V) and AV with slope derived above

# ebvav=[]
# for ind in xrange(ebv.shape[2]):
#   ebvav.append(numpy.dot(tmat,ebv[:,:,ind]))

# ebvav=numpy.array(ebvav)

# ebvav_s = numpy.percentile(ebvav,(50,50-34,50+34),axis=2)

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

ebvav_s = numpy.percentile(ebvav,(50,50-34,50+34),axis=0)

w= numpy.logical_or(ebvav_s[0,:,0] < -.5, numpy.logical_and(ebvav_s[0,:,1]>0.05 , ebvav_s[0,:,0] < -0.2))
print type(data['snlist'])
print numpy.array(data['snlist'])[w]
wefwe
plt.errorbar(ebvav_s[0,:,1], ebvav_s[0,:,0], \
 xerr=(ebvav_s[0,:,1]-ebvav_s[1,:,1], ebvav_s[2,:,1]-ebvav_s[0,:,1]),\
 yerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),fmt='o',alpha=0.4,color='blue')

plt.ylabel(r'$A^F_{V,eff}+ const $')
plt.xlabel(r'$E^F(B-V)_{eff} + const$')
x = numpy.array([-0.09,0.45])
plt.plot(x,2.44*x,color='black',label="slope = 2.44")
# plt.plot(x, rbv[1]+rbv[0]*x,label=r'$R^F={:6.2f}_{{-{:6.2f}}}^{{+{:6.2f}}}$'.format(rbv[0],rbv[0]-mrbv[0],prbv[0]-rbv[0]),color='black')
plt.legend(loc=4)
pp = PdfPages("output_fix3/avebv_synth.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# # plot Rv versus Av for the best fit

# ebvav_r = numpy.percentile(ebvav[:,1,:]/ebvav[:,0,:],(50,50-34,50+34),axis=1)
# plt.errorbar(ebvav_s[0,:,0], ebvav_r[0], \
#  xerr=(ebvav_s[0,:,0]-ebvav_s[1,:,0], ebvav_s[2,:,0]-ebvav_s[0,:,0]),\
#  yerr=(ebvav_r[0]-ebvav_r[1],ebvav_r[2]-ebvav_r[0]),fmt='o',alpha=0.4)
# # plt.scatter(ebvav[0,::1000],ebvav[1,::1000]/ebvav[0,::1000],marker='.')
# plt.ylabel(r'$R^F_{V,eff} $')
# plt.xlabel(r'$E^F(B-V)_{eff} + const$')
# plt.ylim((-1,5))
# pp = PdfPages("output_fix3/avrv_synth.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# import scipy.odr.odrpack as odrpack

# def f(B, x):
#     return B[0]*x + B[1]
# linear = odrpack.Model(f)
# # mydata = odrpack.Data(x, y, wd=1./np.power(sx,2), we=1./np.power(sy,2))
# mydata = odrpack.RealData(ebvav_s[0,:,0], ebvav_s[0,:,1], sx=(ebvav_s[2,:,0]-ebvav_s[1,:,0])/2, sy=(ebvav_s[2,:,1]-ebvav_s[1,:,1])/2)

# myodr = odrpack.ODR(mydata, linear, beta0=[3, 0.])
# myoutput = myodr.run()
# myoutput.pprint()
