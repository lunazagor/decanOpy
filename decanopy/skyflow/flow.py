from __future__ import division

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, Angle, Longitude
from sunpy.coordinates import frames, sun
import star_chart_spherical_projection as scsp
from scipy.stats import gaussian_kde
import random
import pandas as pd
import csv
import argparse
import os
import errno

# supress dubious year warnings 
import warnings
warnings.simplefilter('ignore', UserWarning)


def compute_start_jd(year: str, month: str, location: EarthLocation, dhour: float, dS: float) -> float:
    """
    Compute the starting Julian Date for a sky run.

    Converts the first day of the given month/year BCE to JD, then applies
    a longitude correction for local sidereal time and an optional Stellarium
    offset (dS).

    Parameters
    ----------
    year : str
        Year BCE as a string, e.g. '1300'.
    month : str
        Month as a zero-padded string, e.g. '01'.
    location : EarthLocation
        Observer location on Earth.
    dhour : float
        JD equivalent of one hour (used for longitude correction).
    dS : float
        Stellarium JD offset. Pass 0 if not matching Stellarium.

    Returns
    -------
    float
        Starting Julian Date.
    """
    iso_str = f"-0{year}-{month}-01T00:00:00.000"
    lon_correction = (location.lon.deg / 15.0) * dhour
    return Time(iso_str, scale="local", location=location).jd - lon_correction + dS


# def _normalize_ra(h: int, m: int, s: float) -> tuple[int, int, float]:
#     """Carry overflow in RA seconds/minutes into the next unit."""
#     extra_m, s = divmod(s, 60)
#     extra_h, m = divmod(m + int(extra_m), 60)
#     h = (h + extra_h) % 24
#     return h, m, s


# def precessed_coords(declist: list, year: str) -> tuple[list, list]:
#     """
#     Find the RA and Dec of a list of stars accounting for precession of the equinoxes.
    
#     Uses the Vondrak algorithm via star_chart_spherical_projection. Note that this
#     does not exactly match Stellarium — small differences are expected. See issue #X.

#     Parameters
#     ----------
#     declist : list
#         List of star names to look up.
#     year : str
#         Year BCE as a string, e.g. '1300'.

#     Returns
#     -------
#     obj_list : list of SkyCoord
#         Precessed coordinates for each star.
#     hd_list : list of str
#         Column headers for the output file.
#     """
#     years_since = -2000 - int(year)
#     star_dict = scsp.final_position(declist, year_since_2000=years_since)
    
#     obj_list = []
#     hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
   
#     for name, pos in star_dict.items():
#         parts = pos["RA"].split(".")
#         h, m = int(parts[0]), int(parts[1])

#         # deal with seconds overflow manually 
#         s = float(parts[2][0:2] + "." + parts[2][2:]) if len(parts[2]) > 2 else float(parts[2])
#         if s >= 60:
#             h, m, s = _normalize_ra(h, m, s)
        
#         ra_str = f"{h}h{m}m{s}s" # make legible RA string
#         RA = Angle(ra_str).deg
#         Dec = Angle(pos["Declination"], unit="deg").deg
#         # update object and header lists
#         obj_list.append(SkyCoord(ra=RA, dec=Dec, unit="deg"))
#         hd_list.extend([f"{name} Azimuth", f"{name} Altitude"])

#     return obj_list, hd_list


def precessed_coords(declist, year):
    """
    Find the Dec and RA of a list of decans accounting for precession of the equinoxes.
    """
    years_since = -2000 - int(year)
    print(years_since)
    star_dict = scsp.final_position(declist, year_since_2000=years_since)

    obj_list = []
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]

    for name, pos in star_dict.items():
        ra_parts = pos["RA"].split(".")
        # Defensive: ensure RA has 3 parts
        if len(ra_parts) == 3:
            h, m = ra_parts[0], ra_parts[1]
            s = float(ra_parts[2][0:2] + "." + ra_parts[2][2:]) if len(ra_parts[2]) > 2 else float(ra_parts[2])
            if s >=60.0:
                print(name, pos["RA"], s)
            ra_str = f"{h}h{m}m{s}s"
        else:
            # fallback: treat as decimal hours
            ra_str = f"{pos['RA']}h"
        RA = Angle(ra_str).deg
        Dec = Angle(pos["Declination"], unit="deg").deg
        obj_list.append(SkyCoord(ra=RA, dec=Dec, unit="deg"))
        hd_list.extend([f"{name} Azimuth", f"{name} Altitude"])

    return obj_list, hd_list

def generate_time_grid(start, dhour, d4min):
    days = start + np.arange(0, 365)
    hours = dhour * np.arange(0, 24)
    minutes = d4min * np.arange(0, 15)
    return [
        (day, hour, min)
        for day in days
        for hour in hours
        for min in minutes
    ]

def calc_sun_altaz(temptime, Luxor):
    c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=Time(temptime, format='jd'), observer="earth", frame=frames.Helioprojective)
    frame_altaz = AltAz(obstime=Time(temptime, format='jd'), location=Luxor)
    sun_altaz = c.transform_to(frame_altaz)
    return sun_altaz.T.az, sun_altaz.T.alt

def calc_altaz(ra, dec, loc, time):
    '''
    Find the Altitude and Azimuth of a given star at a given place and time.
    Inputs: 
        ra = right ascension (in decimal hours)
        dec = declination (in degrees)
        loc = location (astropy EarthLocation object)
        time = time (Julian date)
    Outputs:
        alt = altitude of star
        az = azimuth of star
    '''
    # calculations taken from
    # https://www.cloudynights.com/topic/587586-azimuth-altitude-calculation-script/
    #
    # time and location
    obs_time =  Time(time, format = 'jd', location=loc)
    lat = loc.lat
    lon = loc.lon
    # local time
    lst = obs_time.sidereal_time('mean').hour
    gst = lst - (lon.deg / 15.0)
    lha = Angle((gst - ra) * 15 + lon.deg, unit = "deg")
    # alt and az
    alt = altitude(lha, dec, lat)
    az = azimuth(lha, dec, lat)

    return alt, az

def altitude(lha, dec, lat):
    '''
    Find altitude given local hour angle, declination, 
    and latitude of Earth location.
    '''
    a = np.cos(lha.rad)
    b = np.cos(dec.rad)
    c = np.cos(Angle(lat, unit = "deg").rad)
    d = np.sin(dec.rad)
    e = np.sin(Angle(lat, unit = "deg").rad)

    ret = np.arcsin(a*b*c + d*e)
    return Angle(ret, unit="rad").deg

def azimuth(lha, dec, lat):
    '''
    Find azimuth given local hour angle, declination, 
    and latitude of Earth location.
    '''
    a = -1 * np.sin(lha.rad)
    b = np.tan(dec.rad)
    c = np.cos(Angle(lat, unit = "deg").rad)
    d = np.sin(Angle(lat, unit = "deg").rad)
    e = np.cos(lha.rad)

    ret = np.arctan2(a, (b * c - d*e))
    if ret < 0:
        ret += 2 * np.pi
    return Angle(ret, unit="rad").deg

def dS_offset(year):
    '''
    Introduce offset number of days to align Stellarium and Astropy JD. 
    Tested for 1600 to 1100 BCE. 
    '''
    if int(year) > 1499:
        dS_off = -14 
    elif int(year) < 1201:    
        dS_off = -11
    else:
        dS_off = 1 - (int(year) / 100.0)
    return dS_off
