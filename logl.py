#!/usr/bin/env python

import pickle
import csv

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy

def getloglik(filename):
    loglik11= []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        head = spamreader.next()
        index = len(head)-1
        for row in spamreader:
            loglik11.append(row[index])
    return numpy.array(loglik11,dtype='f')

loglik11 = getloglik('output11.csv')
loglik15 = getloglik('output15.csv')

lp2 = loglik11.max()
lp = loglik15.max()


print lp,lp2,lp-lp2


#   Delta   167
#   EWca    168
#   EWSi    168
#   EWl     168
#   D       166
#   k0      166
#   k1      166
#   c       5
#   alpha   5
#   beta    5
#   eta     5
#   delta   3
#   gamma1  5
#   gamma2  5
#   Cc      10

k1 = 167+3*168+3*166+6*5 + 3 + 10
k2 = 167+3*168+2*166+6*5 + 10
n=8.*168

AIC1 = -2*k1 - 2*lp
AIC2 = -2*k2 - 2*lp2
print AIC1-AIC2

c1 = 2*k1*(k1+1)/(n-k1-1)
c2 = 2*k2*(k2+1)/(n-k2-1)

print AIC1-AIC2+c1-c2

print (k1*numpy.log(n) - 2*lp) - (k2*numpy.log(n) - 2*lp2)

plt.hist(loglik11, normed=1, histtype='step', cumulative=-1,label='2 parameters',bins=20)
plt.hist(loglik15, normed=1, histtype='step', cumulative=-1,label='2+1 parameters',bins=20)

plt.legend(loc=9)
plt.xlabel(r'$\log \mathcal{L}$')
plt.ylabel(r'complementary CDF')
pp = PdfPages("output15/lp.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()




