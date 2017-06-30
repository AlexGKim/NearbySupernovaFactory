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
f = open('temp25.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

# Determine the plane approximaion for Fitzpatrick

# Partial derivatives with respect to av and ebv
av=0.1
ebv=av/2.5
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

# # vector of difference
# av_=1.
# ebv_=av/2.
# A1_= f99_band.A_X(r_v=av_/ebv_, ebv=ebv_)
# A2_= f99_band.A_X(r_v=(av_+0.01)/ebv, ebv=ebv_)
# dAdAv_low = (A2_ - A1_)/0.01
# A3_= f99_band.A_X(r_v=av_/(ebv+0.001), ebv=ebv_+0.001)
# dAdebv_low = (A3_ - A1_)/0.001




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

tmat = []
res = []
c_n = []
cs = []
for s in ['gamma','gamma1','rho1']:
  c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
  if s == 'rho1':
    c, cmin, cmax = numpy.percentile(fit[s]/fit[s][:,0][:,None],(50,50-34,50+34),axis=0)

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
kappa3 = tmat[2,0]/tmat[2,1]

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
dum  = numpy.sqrt(cs[2][0]**2+cs[2][2]**2+cs[2][4]**2)
ax.plot([0,cs[2][0]/dum],[0,cs[2][2]/dum],[0,cs[2][4]/dum],label=r'$\delta_X/(\delta_B-\delta_V)$')
dum  = numpy.sqrt(dAdAv[0]**2+dAdAv[2]**2+dAdAv[4]**2)
ax.plot([0,dAdAv[0]/dum],[0,dAdAv[2]/dum],[0,dAdAv[4]/dum],label=r'$a(X)$',ls='--')
dum  = numpy.sqrt(dAdebv[0]**2+dAdebv[2]**2+dAdebv[4]**2)
ax.plot([0,dAdebv[0]/dum],[0,dAdebv[2]/dum],[0,dAdebv[4]/dum],label=r'$b(X)$',ls='--')
# crap = dAdAv + dAdebv*kappa1
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$a(X)+b(X)/{:4.2f}$'.format(1/kappa1),ls=':')
crap = kappa3*dAdAv + dAdebv
dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'${:4.2f}a(X)+b(X)$'.format(kappa3),ls=':',color='black')

# crap = -16*dAdAv - dAdebv
# dum  = numpy.sqrt(crap[0]**2+crap[2]**2+crap[4]**2)
# ax.plot([0,crap[0]/dum],[0,crap[2]/dum],[0,crap[4]/dum],label=r'$-6.8a(X)-b(X)$',ls=':',color='black')

# crap = dAdAv + dAdebv/2.6
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
pp = PdfPages("output25/plane0.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output25/plane1.pdf")
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
pp = PdfPages("output25/plane0BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
ax.view_init(elev=7, azim=-165)
ax.yaxis.set_ticks(numpy.arange(-.75,.76,.25))
ax.xaxis.set_ticks(numpy.arange(-.5,.76,.5))
pp = PdfPages("output25/plane1BVR.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

AVconv = numpy.array(fit['k'][0,:])
EBVconv = numpy.array(fit['k1'][0,:])

# print tmat[0:2,0:2].T
for ind2 in xrange(AVconv.shape[0]):
  # print fit['k'][ind,ind2],fit['k1'][ind,ind2]
  dum= numpy.dot(numpy.array([numpy.median(fit['k'][:,ind2]),numpy.median(fit['k1'][:,ind2])]),tmat[0:2,0:2].T)
  AVconv[ind2]=dum[0]
  EBVconv[ind2]=dum[1]


# plt.hist(AVconv)
# plt.show()

# plt.hist(EBVconv)
# plt.show()

plt.scatter(AVconv,numpy.median(fit['rho1'][:,0][:,None]*fit['R'],axis=0))
plt.show()

plt.scatter(EBVconv,numpy.median(fit['rho1'][:,0][:,None]*fit['R'],axis=0))
plt.show()

wefwe
# r1 = []
# r2 = []
# for ind in xrange(fit['gamma'].shape[0]):
#   tmat = []
#   cs = []
#   for s in ['gamma','rho1']:
#     c = fit[s][ind,:]/((fit[s][ind,1]-fit[s][ind,2]))
#     cs.append(c)

#     y = numpy.array([numpy.dot(c,dAdAv),numpy.dot(c,dAdebv)])

#     ans = numpy.linalg.solve(a,y)

#     tmat.append(ans)

#   tmat = numpy.array(tmat)
#   r1.append(tmat[0,1]/tmat[0,0])
  

# print numpy.percentile(r1,(50,50-34,50+34))
# plt.hist(r1)
# plt.show()





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
