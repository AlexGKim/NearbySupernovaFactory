#!/usr/bin/env python

import pickle
import cPickle
import numpy
import sivel as siv
import sncosmo
import fitz_band
import os

from astropy.io import fits

from specutils.io import read_fits

synlam = numpy.array([[3978.02,4795.35],[4795.35,5780.60]])

synname=['B','V']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))


names168 = [line.split()[0] for line in open('table.txt')]
dic_meta=cPickle.load(open("CABALLOv2/META.pkl"))
ans=[]
for sn in names168:
    meta = dic_meta[sn]
    for nm in meta['spectra']:
        spec=meta['spectra'][nm]
        if abs(spec['salt2.phase']) <2.5:
            name = "CABALLOv2/"+spec['idr.spec_merged']
            spectra = read_fits.read_fits_spectrum1d(name, dispersion_unit='angstrom')


            model0 = sncosmo.TimeSeriesSource(numpy.array([0.,1,2]), spectra.dispersion.value/(1+meta['salt2.Redshift']), \
                numpy.tile(spectra.flux.value,(3,1)))
            dust = sncosmo.F99Dust(r_v=2.5)
            dust.set(ebv=0.01)
            model = sncosmo.Model(source=model0, effects=[dust], effect_names=['host'], effect_frames=['rest'])
            try:

                A= -2.5*numpy.log10(model.bandflux(synbands,0.)/model0.bandflux(synbands,0.))
                ans.append(A[1]/(A[0]-A[1]))
            except:
                pass

ans = numpy.array(ans)
print ans.mean(), ans.std()