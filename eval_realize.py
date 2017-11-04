import csv
import numpy
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc
rc('text', usetex = True)


shit = csv.reader(open('output.5.csv', 'rb'))

rows=[]
for row in shit:
    if len(row) != 1:
        rows.append(row)

# print rows[0][0:5]
# print rows[0][-25:-1]

rows=rows[-1]
rows=rows[2:-50]
rows = numpy.array(rows,dtype='f')

rows = numpy.reshape(rows,(len(rows)/5,5),order='F')

cum=[]
for r in rows:
    cum.append([r[0:i+1].sum() for i in xrange(5)])

cum=numpy.array(cum)


# plt.hist(cum,bins=numpy.arange(0,2.01,0.1))
# plt.show()

norm = cum/cum[:,4][:,None]

plt.hist(norm,bins=numpy.arange(0,1.01,0.05),weights=numpy.array([1/numpy.prod(norm,axis=1)**2.5 for i in xrange(5)]).T,label=['1','2','3','4','5'])
plt.xlim((0.25,1.05))
plt.xlabel(r'normalized cumulative variance in $n$ eignvectors')
plt.ylabel(r'histogram weighted by $\mbox{det}(C_c)^{-5/2}$')
plt.legend(loc=2)
pp = PdfPages('output25/prior.5_normvar.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.clf()

