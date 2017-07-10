#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import f99_band
import matplotlib as mpl
mpl.rcParams['font.size'] = 16

f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

def gettmat():
  # Partial derivatives with respect to av and ebv
  av=0.1
  ebv=0.1/2.5
  A1= f99_band.A_X(r_v=av/ebv, ebv=ebv)
  A2= f99_band.A_X(r_v=(av+0.01)/ebv, ebv=ebv)
  dAdAv = (A2 - A1)/0.01

  A3= f99_band.A_X(r_v=av/(ebv+0.001), ebv=ebv+0.001)
  dAdebv = (A3 - A1)/0.001

  # print 'f99 derivatives'
  # print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdAv)
  # print '{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}, {0[3]:6.2f}, {0[4]:6.2f}'.format(dAdebv)

  # The equation of interest is
  # gammma0 = ans00 F0 + ans01 F1 + res
  # gammma1 = ans10 F0 + ans11 F1 + res
  # where F are the Fitzpatrick vectors (partial derivatives above) and
  # the residues are perpendicular to a and b
  # Note that the gammas are really gamma_X/(gamma_B-gamma_V)

  norm_dAdebv = numpy.dot(dAdebv, dAdebv)
  norm_dAdAv = numpy.dot(dAdAv, dAdAv)
  cross = numpy.dot(dAdebv, dAdAv)

  a = numpy.array([[norm_dAdAv,cross],[cross,norm_dAdebv]])

  tmat = []
  res = []
  # c_n = []
  cs = []
  for s in ['gamma','rho1']:
    c, cmin, cmax = numpy.percentile(fit[s]/((fit[s][:,1]-fit[s][:,2])[:,None]),(50,50-34,50+34),axis=0)
    # print s,
    # print "{:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}".format(c[0],c[1],c[2],c[3],c[4])
    cs.append(c)
    # c_norm = numpy.linalg.norm(c)
    # c_n.append(c_norm)

    y = numpy.array([numpy.dot(c,dAdAv),numpy.dot(c,dAdebv)])
    ans = numpy.linalg.solve(a,y)

    tmat.append(ans)
    ans = c-ans[1]*dAdebv - ans[0]*dAdAv
    res.append(ans)

  tmat = numpy.array(tmat)
  gvec = numpy.array([[4.8,3.89,2.89,2.22,1.59],[-3.83,-4.42,-5.42,-5.15,-4.71]])
  avec = numpy.array([[0.96,1,1,0.97,0.77],[1.77,0.98,0.12,-0.50,-0.53]])
  return tmat


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


#print the matrix and the residues
tmat=gettmat()
myavebv = numpy.dot(tmat.T,Egamma)

print tmat

print "my parameters"
print Egamma
print ugamma

print "F99 parameters"
print myavebv

newmat = numpy.dot(tmat.T, numpy.dot(ugamma, tmat))
print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in newmat])

print myavebv[0]/myavebv[1], numpy.sqrt(newmat[0,0]/myavebv[1]**2 - 2*myavebv[0]*newmat[0,1]/myavebv[1]*3 + newmat[1,1]*myavebv[1]**4)

dum = numpy.random.multivariate_normal(myavebv,newmat,1000)
dum =numpy.percentile(dum[:,0]/dum[:,1],(50,50-34,50+34))
print "{:5.2f}^{{+{:5.2f}}}_{{-{:5.2f}}}".format(dum[0],dum[2]-dum[0],dum[0]-dum[1])
# sanity check
# print "sanity check"
# print "my model"
# gvec = numpy.array([[4.8,3.89,2.89,2.22,1.59],[-3.83,-4.42,-5.42,-5.15,-4.71]])
# avec = numpy.array([[0.96,1,1,0.97,0.77],[1.77,0.98,0.12,-0.50,-0.53]])
# # print avec[0]
# # print avec[1]
# # print tmat[0,0],tmat[0,1]
# print tmat[0,0]*avec[0]+tmat[0,1]*avec[1]
# print tmat[1,0]*avec[0]+tmat[1,1]*avec[1]
# test= gvec[0]*Egamma[0]+ gvec[1]*Egamma[1]
# print  test
# test= avec[0]*myavebv[0]+ avec[1]*myavebv[1]
# print test
# wwefwe


draw = numpy.random.multivariate_normal(myavebv,newmat,100)

fidA1 = f99_band.A_X(r_v=myavebv[0]/myavebv[1], ebv=myavebv[1])
# fidA1 = numpy.delete(fidA1-fidA1[2],2)
temp=[]
for a in draw:
  A1= f99_band.A_X(r_v=a[0]/a[1], ebv=a[1])
  A1 = numpy.delete(A1-A1[2],2)
  # A1 = fidA1 + (A1-fidA1)*2
  temp.append(A1)
temp=numpy.array(temp)
tempperc = numpy.percentile(temp,(50,50-34,50+34),axis=0)

# wefwe
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

axes[0].errorbar(elam, gammam1*ans[0]+deltam1*ans[1],yerr=[fiterr,fiterr],label='Best-fit Model',color='red',fmt='.')
axes[0].errorbar(elam, tempperc[0],yerr=[tempperc[0]-tempperc[1], tempperc[2]-tempperc[0]],label='Projection onto F99',color='blue',fmt='.')
axes[0].errorbar(elam,exv,yerr=[e_exv,e_exv],fmt='o',label='SN2014J',color='black')
# shit=numpy.zeros((temp.shape[0],4))
# for i in xrange(shit.shape[0]):
#   shit[i]=elam
# axes[0].scatter(shit.flatten(), temp.flatten(),label='Best-Inferred F99',color='blue')

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

axes[1].errorbar(elam,-(exv-(gammam1*ans[0]+deltam1*ans[1])), \
  yerr=[numpy.sqrt(e_exv**2+fiterr**2),numpy.sqrt(e_exv**2+fiterr**2)],fmt='.',color='red')
axes[1].errorbar(elam,-(exv-tempperc[0]), \
  yerr=[numpy.sqrt(e_exv**2+((tempperc[2]-tempperc[1])/2)**2),numpy.sqrt(e_exv**2+((tempperc[2]-tempperc[1])/2)**2)],fmt='.', color='blue')
axes[1].axhline(0,linestyle=':')
axes[1].set_ylabel(r'$\Delta E_o(X-V)$')


plt.subplots_adjust(hspace=.07)
pp = PdfPages('output11/sn2014j.pdf')

plt.savefig(pp,format='pdf')
plt.tight_layout()

pp.close()
plt.close()



