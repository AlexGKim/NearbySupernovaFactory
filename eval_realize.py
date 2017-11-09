#!/usr/bin/env python
import csv
import numpy
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import corner
from matplotlib import rc
rc('text', usetex = True)

def cumeval():
    # bins = numpy.exp(numpy.arange(numpy.log(.25),numpy.log(1.01),numpy.log(1/.25)/20)) 
    # bins = 0.25+numpy.log(numpy.arange(numpy.exp(0.),numpy.exp(0.751),(numpy.exp(1.)-numpy.exp(.25))/40))  
    bins = numpy.arange(0.25,1.001,0.02) 
    # files = ['.1','1','4','20']
    # colors = ['red','green','blue','black']
    files = ['.5','20']
    colors = ['blue','red','green','orange']
    ans=[]
    fig, axes = plt.subplots(nrows=len(files)+1,sharex=True)

    ind = 0
    for f,c in zip(files,colors):
        shit = csv.reader(open('output{}.csv'.format(f), 'rb'))

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

        # det = numpy.exp(-2.5*numpy.sum(numpy.log(rows),axis=1))
        # print det.min(), det.max()
        # wefew

        # weight = numpy.array([1/numpy.prod(rows*1e2,axis=1)**2.5 for i in xrange(4)]).T

        shit = axes[ind].hist(norm[:,0:-1],bins=bins, \
            # weights=weight, \
               # label=('{}'.format('diag'),'_nolegend_', '_nolegend_','_nolegend_'), \
           label=['1','2','3','4'],
            histtype='step',\
            # linestyle = ('-', '--', '-.', ':'),\
            color=colors,normed=True,log=True)

        axes[ind].set_ylim(1e-1,1e2)
        axes[ind].set_title(r'LKJ $\eta={}$'.format(f),fontsize=9)
        ind = ind+1

    #diagonal
    cum = []
    for i in xrange(10000):
        cumon=[]
        for j in xrange(5):
            temp = -1
            while temp < 0:
                temp = numpy.random.standard_cauchy()
                temp = temp*.1
                temp = temp+.1
            cumon.append(temp)
        cumon = numpy.sort(cumon)
        cumon = [cumon[i] for i in xrange(4,-1,-1)]
        cum.append(cumon)
    cum=numpy.array(cum)
    rows= cum**2

    cum=[]
    for r in rows:
        cum.append([r[0:i+1].sum() for i in xrange(5)])

    cum=numpy.array(cum)
    norm = cum/cum[:,4][:,None]


    # plt.hist(norm[:,0:-1],bins=numpy.arange(0,1.01,0.02), \
    #     # weights=numpy.array([1/numpy.prod(norm,axis=1)**2.5 for i in xrange(4)]).T, \
    #     label=['{} 1'.format('diag'),'{} 2'.format('diag'),'{} 3'.format('diag'),'{} 4'.format('diag')],histtype='step',\
    #     color=['orange','orange','orange','orange'],normed=True)
    axes[-1].hist(norm[:,0:-1],bins=bins, \
        # weights=numpy.array([1/numpy.prod(norm,axis=1)**2.5 for i in xrange(4)]).T, \
        # label=('{}'.format('diag'),'_nolegend_', '_nolegend_','_nolegend_'), \
        label=['1','2','3','4'],
        histtype='step',\
        # linestyle = ('-', '--', '-.', ':'),\
        color=colors,normed=True,log=True)
    axes[-1].set_title('diag',fontsize=11)
    axes[-1].set_ylim(1e-1,1e2)

    plt.xlim((0.25-.01,1.01))
    plt.xlabel(r'normalized cumulative variance in $n$ eignvectors')
    # plt.ylabel(r'histogram weighted by $\mbox{det}(C_c)^{-5/2}$')
    plt.legend(loc=2,fontsize=9)
    # plt.xscale('logit')
    pp = PdfPages('output25/prior.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.clf()



def det(weight=False):
    # bins = numpy.exp(numpy.arange(numpy.log(.25),numpy.log(1.01),numpy.log(1/.25)/20)) 
    # bins = 0.25+numpy.log(numpy.arange(numpy.exp(0.),numpy.exp(0.751),(numpy.exp(1.)-numpy.exp(.25))/40))  
    bins = numpy.arange(-50,5,2) 
    # files = ['.1','1','4','20']
    # colors = ['red','green','blue','black']
    files = ['.5','1','4','20']
    colors = ['red','orange','green','blue','black']
    ans=[]

    ind = 0
    for f,c in zip(files,colors):
        shit = csv.reader(open('output{}.csv'.format(f), 'rb'))

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

        lndet = numpy.sum(numpy.log(rows),axis=1)

        if weight:
            w=numpy.prod(rows,axis=1)**(-2.5)
        else:
            w=None


        shit = plt.hist(lndet,bins=bins, \
            label=r'$\eta ={}$'.format(f), \
            histtype='step',\
            # linestyle = ('-', '--', '-.', ':'),\
            color=c,normed=True, weights=w)


    #diagonal
    cum = []
    for i in xrange(10000):
        cumon=[]
        for j in xrange(5):
            temp = -1
            while temp < 0:
                temp = numpy.random.standard_cauchy()
                temp = temp*.1
                temp = temp+.1
            cumon.append(temp)
        # cumon = numpy.sort(cumon)
        # cumon = [cumon[i] for i in xrange(4,-1,-1)]
        cum.append(cumon)
    cum=numpy.array(cum)
    rows= cum**2
    lndet = numpy.sum(numpy.log(rows),axis=1)

    if weight:
        w=numpy.prod(cum,axis=1)**(-2.5)
    else:
        w=None

    plt.hist(lndet,bins=bins, \
        # weights=numpy.array([1/numpy.prod(norm,axis=1)**2.5 for i in xrange(4)]).T, \
        # label=('{}'.format('diag'),'_nolegend_', '_nolegend_','_nolegend_'), \
        label='diagonal', \
        histtype='step',\
        # linestyle = ('-', '--', '-.', ':'),\
        color=colors[-1],normed=True,weights=w)


    # plt.xlim((0.25-.01,1.01))
    plt.xlabel(r'ln(det)')
    # plt.ylabel(r'histogram weighted by $\mbox{det}(C_c)^{-5/2}$')
    plt.legend(loc=2)
    ext=''
    if weight:
        ext='w'
    pp = PdfPages('output25/priordet{}.pdf'.format(ext))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.clf()

# det()
# det(True)

# wefwe
cumeval()

wefwe
def elements():
    shit = csv.reader(open('diag.csv', 'rb'))

    rows=[]
    for row in shit:
        if len(row) != 1:
            rows.append(row)

    # print rows[0][0:5]
    # print rows[0][-25:-1]

    rows=rows[-1]
    rows=rows[2:-50]
    rows4 = numpy.array(rows,dtype='f')


    shit = csv.reader(open('diag.05.csv', 'rb'))

    rows=[]
    for row in shit:
        if len(row) != 1:
            rows.append(row)

    # print rows[0][0:5]
    # print rows[0][-25:-1]

    rows=rows[-1]
    rows=rows[2:-50]
    rows5 = numpy.array(rows,dtype='f')

    shit = csv.reader(open('diag.20.csv', 'rb'))

    rows=[]
    for row in shit:
        if len(row) != 1:
            rows.append(row)

    # print rows[0][0:5]
    # print rows[0][-25:-1]

    rows=rows[-1]
    rows=rows[2:-50]
    rows20 = numpy.array(rows,dtype='f')

    plt.hist([rows5, rows4, rows20],histtype='step',label=['0.5','4','20'],normed=True,bins=50)
    plt.legend()
    plt.show()



    rows = numpy.reshape(rows,(len(rows)/10,10),order='F')
    corner.corner(rows)

    m, l, h = numpy.percentile(rows,(50,50-34,50+34),axis=0)

    print r'\hat{{U}} & ${:6.3f}^{{{:6.3f}}}_{{{:6.3f}}} & ${:6.3f}^{{{:6.3f}}}_{{{:6.3f}}} & ${:6.3f}^{{{:6.3f}}}_{{{:6.3f}}} & ${:6.3f}^{{{:6.3f}}}_{{{:6.3f}}}'.format(m[0],m[0]-l[0],h[0]-m[0], \
        m[1],m[1]-l[1],h[1]-m[1],m[2],m[2]-l[2],h[2]-m[2],m[3],m[3]-l[3],h[3]-m[3],m[4],m[4]-l[4],h[4]-m[4])
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


