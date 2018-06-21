#!/bin/env python

# code to read and shape data for Gerard's like analysis.

import numpy as N
import matplotlib as plt
import pylab as P
from scipy import stats
from ToolBox.Astro import Extinction
import scipy as S

# get ubvri at max for good SN

def spec_near_phase(meta,sn,origin=0, delta=2.5):
    allphases = [(N.abs( meta[sn]['spectra'][spec]['salt2.phase'])-origin,spec) for spec in meta[sn]['spectra'] ]
    closer = sorted(allphases)[0]
    if closer[0] < delta:
        return closer[1]
    else:
        return False
  
def build_ubvri_phreno(meta,ubvri,phreno=None):
    # meta needed for Caballov2 flags
    # good meta SN
    snlist=list()
    obsveclist=list()
    obscovlist=list()
    bands = ['U','B','V','R','I']
    SI = ['EWCaIIHK','EWSiII4000']
    bad_sn = ['SN2012cu','SNF20061108-001'] + ['SNF20080905-005','SNNGC7589']
    #bad_sn=[]

    for sn in ubvri:
        if sn in bad_sn:
            continue
        if ubvri[sn]['idr.subset'] in {'validation','training'}:
            # is sn in phreno ?
            if phreno:
                if sn not in phreno:
                    print sn+" is not in Phreno..."
                    continue
                spec = spec_near_phase(ubvri,sn)
                if not spec:
                    continue
                #print sn,spec
                try:
                    if meta[sn]['spectra'][spec]['procB.Quality'] !=1:
                        continue
                except: # key or sn not present
                    continue
                try:
                    if meta[sn]['spectra'][spec]['procR.Quality'] !=1:
                        continue
                except: # key or sn not present
                    continue

            snlist.append(sn)
            obsvec=list()
            if phreno:
                cov=N.zeros((len(bands)+len(SI),len(bands)+len(SI)))
                for i,si in enumerate(SI):
                    obsvec.append(phreno[sn]['spectra'][spec]['phrenology.'+si])
                    cov[i,i] = phreno[sn]['spectra'][spec]['phrenology.'+si+'.err']**2
                skip=len(SI)
            else:
                cov=N.zeros((len(bands),len(bands)))
                skip=0
            for i,b in enumerate(bands):
                obsvec.append(ubvri[sn]['AbsoluteMagMax.'+b])
                cov[i+skip,i+skip] = ubvri[sn]['AbsoluteMagMax.'+b+'.err']**2
                for j,b2 in enumerate(bands):
                    if j>i:
                        cov[i+skip,j+skip]=cov[j+skip,i+skip]=ubvri[sn]['AbsoluteMagMax.'+b+'_'+b2+'.cov']
            obscovlist.append(cov)
            obsveclist.append(N.array(obsvec))
    return snlist, obsveclist, obscovlist


def center_and_condition(obs, covs):
    # pre-condition the matrix in order to have smoother computations.
    # returns the center and the scale to be applied to find the original data from the conditioned ones
    # i.e. input = output x scale + center
    obs=N.array(obs)
    covs=N.array(covs)
    center = N.mean(obs,axis=0)
    scale = N.std(obs,axis=0)
    # data transform :
    transobs = (obs-center)/scale
    transcov = covs/N.outer(scale,scale)
    return center, scale, transobs, transcov

if __name__ == "__main__2":

    import cPickle
    ubvri = cPickle.load(open('ubvri.pkl'))
    phreno = cPickle.load(open('phrenology_2016_12_01_CABALLOv1.pkl'))
    meta=cPickle.load(open('SNF-0203-CABALLO2/META.pkl'))

    snlist,obs,cov=build_ubvri_phreno(meta,ubvri,phreno)
    obs=N.array(obs)
    cov=N.array(cov)
    center,scale,tobs,tcov=center_and_condition(obs,cov)
    tweight=N.array([N.linalg.inv(t) for t in tcov])

