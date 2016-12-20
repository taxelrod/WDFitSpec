"""
Process set of Hubeny WD model spectra into Ulyss Type ? TGM file
"""
import numpy as np
import os
import re
from astropy.io import ascii
import astropy.io.fits as fits
from scikits.fitting import Polynomial2DFit

#
# parse a name like t110g60l.sp.11 into Teff=110000 logg=6.00
#
def HubenyNameToParams(name):
    m = re.search(r"t(\d+)g(\d+).*", name)
    if m:
        teff = int(m.group(1))*100.0
        logg = int(m.group(2))*0.10
        return(teff, logg)
    else:
        return(None, None)
#
# Go through directory of Hubeny models, and return array of (Teff, logg)
# It is assumed that the directory contains ONLY valid model files!
#
def GetGrid(dirName):
    fileList = os.listdir(dirName)
    nSpectra = len(fileList)
    grid = np.zeros((nSpectra,2))
    for (i, f) in enumerate(fileList):
        (teff, logg) = HubenyNameToParams(f)
        grid[i,0] = np.log(teff)
        grid[i,1] = logg
        fileList[i] = dirName + '/' + fileList[i]

    return (fileList, grid)

def ReadGridData(fileList):
    for (i, fileName) in enumerate(fileList):
        gridPtData = ascii.read(fileName)
        if i==0:
            wl = gridPtData.columns[0]
            nSpec = len(fileList)
            nWl = len(wl)
            gridSpectra = np.zeros((nSpec, nWl))
        gridSpectra[i,:] = gridPtData.columns[1]

    return(wl, gridSpectra)

#
# Fit 2D polynomial to gridSpectra at wavelength index i over the grid
#
def Fit2DPoly(grid, gridSpectra, i, degree=3, fullInfo=False):
    interp = Polynomial2DFit((degree,degree))
    interp.fit(np.transpose(grid), gridSpectra[:,i])
    checkVals = interp(np.transpose(grid))
    if fullInfo:
        return interp.poly, checkVals, interp
    else:
        return interp.poly

def Eval2DPoly(interp, poly, xlo, xhi, nx, ylo, yhi, ny):
    grid_x, grid_y = np.mgrid[xlo:xhi:(nx*1j), ylo:yhi:(ny*1j)]
    x = grid_x.ravel()
    y = grid_y.ravel()
    x2 = x*x
    x3 = x2*x
    y2 = y*y
    y3 = y2*y
    vals = interp(np.vstack((grid_x.ravel(), grid_y.ravel())))

    checkVals = poly[0]*x3*y3 + poly[1]*x3*y2 + poly[2]*x3*y + poly[3]*x3 + poly[4]*x2*y3 + poly[5]*x2*y2 + poly[6]*x2*y + poly[7]*x2 + poly[8]*x*y3 + poly[9]*x*y2 + poly[10]*x*y + poly[11]*x + poly[12]*y3 + poly[13]*y2 + poly[14]*y + poly[15]

    return vals.reshape((nx, ny)), checkVals.reshape((nx, ny))

def RemapPoly(poly):
    if len(poly) != 16:
        print 'Error: wrong length poly in RemapPoly ', len(poly)
        return None
    ulyPoly = np.zeros(16)
    ulyPoly[12] = poly[0]
    return ulyPoly

def GenUlyPoly(grid, gridSpectra):
    nWl = gridSpectra.shape[1]
    fitData = np.zeros((16, nWl))
    for n in range(nWl):
        fitData[:, n] = Fit2DPoly(grid, gridSpectra, n)

    return fitData
        

def WriteUlyTGM(fileName, wl, polyArray):
#
# check to be sure that polyArray is consistent with wl
#
    nWl = len(wl)
    if not polyArray.shape==(16, nWl):
        print 'Could not write TGM file: Incompatible sizes ', nWl, polyArray.shape
        return None

    # need accurate calculation of logarithmic step size
    deltaLogWl = (np.log(wl[nWl-1])-np.log(wl[0]))/float(nWl-1)
    priHdr = fits.Header()
    priHdr['CRPIX1'] = 1.
    priHdr['CRVAL1'] = np.log(wl[0])
    priHdr['CD1_1'] = deltaLogWl
    priHdr['CDELT1'] = deltaLogWl
    priHdr['H_AXIS1'] = 2
    priHdr['CTYPE1'] = 'AWAV-LOG'
    priHdr['INTRP_V'] = '4'
    priHdr['ULY_TYPE'] = 'TGM'

    hdu = fits.PrimaryHDU(polyArray, header=priHdr)
    hdu.writeto(fileName)
    

                    

    

    
