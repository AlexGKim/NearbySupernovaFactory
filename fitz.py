#!/usr/bin/env python
import pickle
import numpy
import pystan
import sncosmo

synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))

model_nodust = sncosmo.Model(source='hsiao')
flux_nodust = model_nodust.bandflux(synbands,0.)

av = numpy.arange(0,1.8,0.05)
rv = numpy.arange(1.,5.,0.05)

avs=[]
ebvs=[]
rvs=[]
AX = []

for a in av:
    for r in rv:
        dust = sncosmo.F99Dust(r_v=r)
        dust.set(ebv=a/r)
        model = sncosmo.Model(source='hsiao', effects=[dust], effect_names=['host'], effect_frames=['rest'])
        AX.append(-2.5*numpy.log10(model.bandflux(synbands,0.)/flux_nodust))
        avs.append(a)
        ebvs.append(a/r)
        rvs.append(r)

avs = numpy.array(avs)
ebvs = numpy.array(ebvs)
AX = numpy.array(AX)
rvs=numpy.array(rvs)

data = {'D': avs.size, 'AV': avs, 'EBV': ebvs, 'AX': AX}

av_ = numpy.array([0.01,0.01,0.01+1e-4])
ebv_ = numpy.array([0.01/2.5, 0.01/2.5+1e-4, 0.01/2.5])
rv_=av_/ebv_

ans_=[]
for av0, ebv0, rv0 in zip(av_,ebv_,rv_):
    dust = sncosmo.F99Dust(r_v=rv0)
    dust.set(ebv=ebv0)
    model = sncosmo.Model(source='hsiao', effects=[dust], effect_names=['host'], effect_frames=['rest'])
    ans_.append(2.5*numpy.log10(model.bandflux(synbands,0.)))

ans_ = numpy.array(ans_)

dum1 = (ans_[0]-ans_[2])/1e-4
dum2= (ans_[0]-ans_[1])/1e-4

init1 = {'a' : numpy.array([dum1,numpy.zeros(5),dum2,numpy.zeros(5),numpy.zeros(5),numpy.zeros(5),numpy.zeros(5),numpy.zeros(5),numpy.zeros(5)])}

sm = pystan.StanModel(file='fitz.stan')
fit = sm.sampling(data=data, iter=2000, chains=4,init=[init1,init1,init1,init1])

ans = fit.extract()

amed= numpy.median(fit['a'], axis=0)
diff = AX - (amed[0][None,:]*avs[:,None]+ amed[1][None,:] * avs[:,None]**2 \
    +amed[2][None,:]*ebvs[:,None]+ amed[3][None,:] * ebvs[:,None]**2 \
    +amed[4][None,:] * (avs*ebvs)[:,None] \
    +amed[5][None,:] * (avs**3)[:,None] \
    +amed[6][None,:] * (ebvs**3)[:,None] \
    +amed[7][None,:] * ((avs**2)*ebvs)[:,None] \
    +amed[8][None,:] * (avs*(ebvs**2))[:,None] \
    )

print numpy.max(numpy.abs(diff))
# 0.00589190110442

output = open('fitz.pkl','wb')
pickle.dump(amed, output, protocol=2)
output.close()
