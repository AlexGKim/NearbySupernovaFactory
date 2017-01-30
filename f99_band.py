import sncosmo
import numpy

#efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])
synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))

model_nodust = sncosmo.Model(source='salt2')
flux_nodust = model_nodust.bandflux(synbands,0.) 

def A_X(r_v=3.1, ebv=1.):

    dust = sncosmo.F99Dust(r_v=r_v)
    dust.set(ebv=ebv)
    model = sncosmo.Model(source='salt2', effects=[dust], effect_names=['host'], effect_frames=['rest'])
    return -2.5*numpy.log10(model.bandflux(synbands,0.)/flux_nodust)