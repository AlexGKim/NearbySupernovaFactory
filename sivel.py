import pickle
import cPickle
import numpy


# two color parameter model
def sivel(data):
   dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

   dic_meta=cPickle.load(open("SNF-0203-CABALLO2/META.pkl"))

   sivel=[]
   sivel_err=[]
   x1 = []
   x1_err = []
   zcmb = []
   zerr= []
   for sn in data['snlist']:
      if sn in dic_meta.keys() and sn in dic_phreno.keys():
         meta = dic_meta[sn]
         vSiII_6355_lbd=0.
         vSiII_6355_lbd_err=0.
         counter  = 0
         for sp in dic_phreno[sn]["spectra"]:
            if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase']) < 2.5 and numpy.isfinite(dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]):
               vSiII_6355_lbd += dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
               vSiII_6355_lbd_err += 1/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
               counter +=1
         if counter !=0:
            sivel.append(vSiII_6355_lbd / vSiII_6355_lbd_err)
            sivel_err.append(1./numpy.sqrt(vSiII_6355_lbd_err))
         else:
            sivel.append(float('nan'))
            sivel_err.append(float('nan'))
         x1.append(meta['salt2.X1'])
         x1_err.append(numpy.sqrt(meta['salt2.CovX1X1']))
         zcmb.append(meta['host.zcmb'])
         zerr.append(meta['host.zhelio.err'])
      else:
         sivel.append(float('nan'))
         sivel_err.append(float('nan'))
         x1.append(float('nan'))
         x1_err.append(float('nan')) 
         zcmb.append(float('nan'))
         zerr.append(float('nan'))

   sivel = numpy.array(sivel)
   sivel_err = numpy.array(sivel_err)
   x1 = numpy.array(x1)
   x1_err = numpy.array(x1_err)
   zcmb = numpy.array(zcmb)
   zerr = numpy.array(zerr)

   return sivel,sivel_err,x1,x1_err,zcmb,zerr