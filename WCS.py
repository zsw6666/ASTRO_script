import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle

import os


def CentralWCS(RA,DEC):
    ra,dec=Angle(RA+' hours'),Angle(DEC+' hours')

    degree_ra,degree_dec=ra.degree,dec.degree
    return degree_ra,degree_dec



def CentralPhysics(Lx,Ly):

    CentralX,CentralY=int(Lx/2.),int(Ly/2.)
    return CentralX,CentralY


def MoidifyWCS(path,fitslistname):

    os.chdir(path)

    fitslist=open(fitslistname).readlines()

    for i in range(len(fitslist)):
        wcsimg=fits.open(fitslist[i][:-1])
        RA,DEC=wcsimg[0].header['OBSRA'],wcsimg[0].header['OBSDEC']
        ref_ra,ref_dec=CentralWCS(RA,DEC)
        Lx,Ly=wcsimg[0].header['NAXIS1'],wcsimg[0].header['NAXIS2']
        ref_x,ref_y=CentralPhysics(Lx,Ly)
        fits.setval(fitslist[i][:-1], 'CTYPE1', value='RA---TAN')
        fits.setval(fitslist[i][:-1], 'CTYPE2', value='DEC---TAN')
        fits.setval(fitslist[i][:-1],'CRVAL1',value=ref_ra)
        fits.setval(fitslist[i][:-1], 'CRVAL2', value=ref_dec)
        fits.setval(fitslist[i][:-1], 'CRPIX1', value=ref_x)
        fits.setval(fitslist[i][:-1], 'CRPIX2', value=ref_y)
        wcsimg.close()






