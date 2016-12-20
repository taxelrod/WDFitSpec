#!/usr/bin/env python
"""
Take an .flm file and write a fits
table that ULySS can take as input
"""
import pyfits as pf
import numpy as np
import sys

def writeFitsTbl(fileName, wl, spec, var, mask):
    ivar = 1./var
    c1 = pf.Column(name='SPEC', format='E', array=spec)
    c2 = pf.Column(name='LAMBDA', format='E', array=wl)
    c3 = pf.Column(name='IVAR', format='E', array=ivar)
    c4 = pf.Column(name='ORMASK', format='I', array=mask)
    tblHdu = pf.BinTableHDU.from_columns([c1, c2, c3, c4])
    tblHdr = tblHdu.header
    tblHdu.writeto(fileName)

def readFLM(fileName):
    flm = np.loadtxt(fileName,skiprows=1)
    return flm

if __name__ == "__main__":
    inFileName = sys.argv[1]
    outFileName = sys.argv[2]
    flm = readFLM(inFileName)
    wl = flm[:,0]
    flux = flm[:,1]
    fluxerr = flm[:,2]
    mask = np.zeros(len(wl), int)
    
    writeFitsTbl(outFileName, wl, flux, fluxerr, mask)
