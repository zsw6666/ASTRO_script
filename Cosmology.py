import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck13 as cosmos

def Angle_distance(z):
    return cosmos.angular_diameter_distance(z)

def Luminosity_distance(z):
    return cosmos.luminosity_distance(z)

def Comoving_distance(z):
    return cosmos.comoving_distance(z)

def Arcsec2rad(angle):
    arc=angle*u.arcsec
    ra=arc.to(u.rad)

    return ra.value