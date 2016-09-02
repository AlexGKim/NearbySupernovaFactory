#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc

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

gammam1 = numpy.array([0.643, 0.336,-0.228,-.447])
deltam1 = numpy.array([-.393, -.224, -.056, -0.147])

fisher = numpy.array([[sum(gammam1**2/e_exv), sum(gammam1*deltam1/e_exv)],[sum(gammam1*deltam1/e_exv), sum(deltam1**2/e_exv)]])
u=numpy.linalg.inv(fisher)
g = numpy.array([sum(exv*gammam1/e_exv), sum(exv*deltam1/e_exv)])
ans = numpy.dot(u,g)
print ans, u

fiterr = numpy.sqrt(gammam1**2 * u[0,0] + 2*deltam1*gammam1*u[0,1] + deltam1**2 * u[1,1])

import sncosmo


plt.errorbar(elam,exv,yerr=[e_exv,e_exv],fmt='.',label='SN2014J',color='black')
plt.errorbar(elam, gammam1*ans[0]+deltam1*ans[1],yerr=[fiterr,fiterr],label='Best-fit Model',color='red',fmt='.')
print  gammam1*ans[0], deltam1*ans[1]

lambdas = numpy.arange(3500.,8000,100)
rvs=[1.4, 3.1]

for rv in rvs:
    f99 = sncosmo.F99Dust(r_v =rv)
    f99.set(ebv=1.37)
    A_ = f99.propagate(lambdas,1.) / f99.propagate(numpy.array([5477.]),1.)[0]
    A_=-2.5*numpy.log10(A_)
    # A_ = sncosmo._extinction.ccm89(lambdas, 1.37, rv) - sncosmo._extinction.ccm89(numpy.array([5477.]), 1.37, rv)[0]
    # norm  = sncosmo._extinction.ccm89(numpy.array([5477.]), 1., rv)
    # A_ = A_/norm[0]-1
    plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))

plt.legend()
plt.xlim((3200,8000))
plt.ylim((-2,3))
plt.ylabel(r'E(B-V)')
plt.xlabel(r'Wavelength (\AA)')
pp = PdfPages('output11/sn2014j.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
