#!/usr/bin/env python

# perfect standard candle

import pickle
import cPickle
import numpy
import pystan
import sivel as sivel2
import sncosmo
print sncosmo.__version__
# two color parameter model
numpy.random.seed(seed=1)

synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()


f = open('temp11.pkl','rb')
(fit,_) = pickle.load(f)
f.close()

sivel,sivel_err,x1,x1_err,zcmb,zerr = sivel2.sivel(data)

use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

#The values to use

EW = numpy.median(fit['EW'],axis=0)
sivel = numpy.median(fit['sivel'],axis=0)


nsne, nmags = mag_obs.shape

for seed in xrange(10):

    numpy.random.seed(seed=seed)


    intrinsic=[]
    for i in xrange(nsne):
        intrinsic.append(numpy.zeros(5)+numpy.random.normal(0,0.08))
    intrinsic = numpy.array(intrinsic)

    av = fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R']
    av=numpy.median(av,axis=0)
    # ab = fit['gamma'][:,1][:,None]*fit['k'] + fit['rho1'][:,1][:,None]*fit['R']
    # ab=numpy.median(ab,axis=0)

    ebv = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + \
        (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']
    ebv=numpy.median(ebv,axis=0)
    ebv=ebv-ebv.min()

    # wef
    # import matplotlib.pyplot as plt
    # plt.hist(ebv-ebv.min())
    # plt.show()

    rv = numpy.random.lognormal(numpy.log(2.5),0.05, size=nsne)
    # import matplotlib.pyplot as plt
    # plt.scatter(rv, rv*ebv)
    # plt.show()
    # wfew
    mag_obs=[]
    EW_obs=[]
    sivel_obs=[]
    for i in xrange(nsne):
        dust = sncosmo.F99Dust(r_v=rv[i])
        dust.set(ebv=ebv[i])
        # model0 = sncosmo.Model(source='hsiao')
        # model0.parameters[2]=x0
        # model0.parameters[3]=x1
        # model0.parameters[4]=c
        model = sncosmo.Model(source='salt2', effects=[dust], effect_names=['host'], effect_frames=['rest'])
        # model.parameters[2]=x0
        # model.parameters[3]=x1
        # model.parameters[4]=c
        # mag_obs.append(numpy.random.multivariate_normal(intrinsic[i,:]+2.5*numpy.log10(model0.bandflux(synbands,0.)/model.bandflux(synbands,0.)),mag_cov[i]))
        mag_obs.append(numpy.random.multivariate_normal(intrinsic[i,:]+model.bandmag(synbands,'vega',0.),mag_cov[i]))
        EW_obs.append(numpy.random.multivariate_normal(EW[i,:],EW_cov[i]))
        sivel_obs.append(numpy.random.normal(sivel[i],sivel_err[i]))

    mag_obs=numpy.array(mag_obs)
    EW_obs=numpy.array(EW_obs)
    sivel_obs = numpy.array(sivel_obs)

    # # renormalize the data
    EW_mn = EW_obs.mean(axis=0)
    EW_renorm = (EW_obs - EW_mn)

    mag_mn = mag_obs.mean(axis=0)
    mag_renorm  = mag_obs-mag_mn

    # import matplotlib.pyplot as plt
    # plt.scatter(mag_renorm[:,1]-mag_renorm[:,2],mag_renorm[:,2])
    # plt.show()

    sivel_mn = sivel.mean()
    sivel_renorm = sivel-sivel_mn

    data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
       'sivel_obs': sivel_renorm, 'sivel_err': sivel_err}


    output = open('simdata_gerard{}.pkl'.format(seed),'wb')
    pickle.dump(data, output, protocol=2)
    output.close()
