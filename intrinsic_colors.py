#!/usr/bin/env python

import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo
import scipy
import cPickle
import matplotlib as mpl
import sivel
mpl.rcParams['font.size'] = 18
ext=''

f = open('temp18'+ext+'.pkl','rb')
(fit18,_) = pickle.load(f)



f = open('temp11'+ext+'.pkl','rb')
(fit11,_) = pickle.load(f)

# for key in fit.keys():
#     print key, fit[key].min(), fit[key].max()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel,sivel_err,x1,x1_err, _, _ = sivel.sivel(data)
use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
x1=x1[use]
x1_err = x1_err[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn

nsne, nmags = mag_obs.shape
color_obs = numpy.zeros((nsne,nmags-1))
color_obs[:,0] = mag_renorm[:,0]- mag_renorm[:,2]
color_obs[:,1] = mag_renorm[:,1]- mag_renorm[:,2]
color_obs[:,2] = mag_renorm[:,3]- mag_renorm[:,2]
color_obs[:,3] = mag_renorm[:,4]- mag_renorm[:,2]

pkl_file.close()

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))

colors = fit18['c']+mag_mn[None,:]
colors  = colors - colors[:,2][:,None]
colors18 = numpy.delete(colors, 2,axis=1)


colors = fit11['c']+mag_mn[None,:]
colors  = colors - colors[:,2][:,None]
colors11 = numpy.delete(colors, 2,axis=1)


from matplotlib.ticker import FuncFormatter, MaxNLocator
def format_fn(tick_val, tick_pos):
    if int(tick_val) in numpy.arange(4):
        return labels[int(tick_val)]
    else:
        return ''

colors18_perc = numpy.percentile(colors18,(50,50-34,50+34),axis=0)
colors11_perc = numpy.percentile(colors11,(50,50-34,50+34),axis=0)
print colors18_perc[0], colors11_perc[0]
ecolor = colors11_perc[0] - colors18_perc[0]
decolor = numpy.sqrt((colors18_perc[2]-colors18_perc[1])**2/4 + (colors11_perc[2]-colors11_perc[1])**2/4)
fig = plt.figure()
ax = fig.add_subplot(111)
labels = list('UBRI')
ax.errorbar(numpy.arange(4),ecolor,yerr=decolor,label='Model Difference')
ax.xaxis.set_major_formatter(FuncFormatter(format_fn))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlim((-0.5,3.5))


snmod='salt2'

synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))

model_nodust = sncosmo.Model(source=snmod)
flux_nodust = model_nodust.bandflux(synbands,0.)
rv=numpy.arange(2,4,.5)
for r in rv:
    dust = sncosmo.F99Dust(r_v=r)
    dust.set(ebv=0.3)
    model = sncosmo.Model(source=snmod, effects=[dust], effect_names=['host'], effect_frames=['rest'])
    AX=-2.5*numpy.log10(model.bandflux(synbands,0.)/flux_nodust)
    AX  = AX - AX[2]
    AX = numpy.delete(AX, 2)
    ax.plot(numpy.arange(4),AX,label=r'$R={}$'.format(r))
    ax.set_ylabel(r'$E(X-V)$')
    ax.set_xlabel(r'Band $X$')
plt.legend(fontsize=15)
pp = PdfPages('output18'+ext+'/model1118ebv.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


labels = list('UBVRI')
def format_fn2(tick_val, tick_pos):
    if int(tick_val) in numpy.arange(5):
        return labels[int(tick_val)]
    else:
        return ''


rs = numpy.percentile((fit11['c'][:16000,:]-fit18['c'])/(colors11[:16000,1]-colors18[:,1])[:,None],(50,50-34,50+34),axis=0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(numpy.arange(5), rs[0], yerr=[rs[0]-rs[1],rs[2]-rs[0]],fmt='o')
ax.xaxis.set_major_formatter(FuncFormatter(format_fn2))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlim((-0.5,4.5))
ax.set_xlabel(r'Band $X$')
ax.set_ylabel(r'$A_X / E(B-V)$ of F99 relative to Linear Model')
pp = PdfPages('output18'+ext+'/comptolinear.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wefwe
c11mn = numpy.median(colors11,axis=0)
c18mn = numpy.median(colors18,axis=0)
print c11mn-c18mn

wef


figure = corner.corner(colors,labels=[r"$U_0-V_0$",r"$B_0-V_0$",r"$R_0-V_0$",r"$I_0-V_0$"])
pp = PdfPages('output18'+ext+'/col_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


