"""Functions for procedurally generating skies."""

from __future__ import division

import random
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from scipy.stats import gaussian_kde

# supress dubious year warnings 
import warnings
warnings.simplefilter('ignore', UserWarning)

# paths
from decanopy.config import DEFAULT_STAR_DATA, OBS_LAT


def generate_mags_from_Hdist(num, data_path: DEFAULT_STAR_DATA):
    '''Generate a number of magnitudes from a kernel density estimator based on the distribution of star magnitudes in
    the Hipparcos data stored in decanopy/data/star_data.csv. 

    The input num is the number of magnitudes to generate.
    The output is an array of num magnitudes.
    '''
    # read file 
    try:
        name_df = pd.read_csv(data_path, index_col=None, header=0, names=['Name', 'RA', 'Dec', 'Mag'])
    except Exception as e:
        raise RuntimeError(f"Could not read star data from {data_path}: {e}")
    data = name_df['Mag'].values #np.random.normal(size=1000)
    # create kernel density estimator pdf
    x_grid = np.linspace(min(data), max(data), num)
    kdepdf = kde(data, x_grid, bandwidth=0.1)
    random_from_kde = generate_rand_from_pdf(kdepdf, x_grid, num)
    # return an array of num magnitudes
    return random_from_kde

def generate_star(lat: OBS_LAT):
    '''
    Generate a single star at a random position in the sky.
    The star is uniformly distributed across a hemisphere, with a minimum declination visible at some optional latitude. 
    Currently defaults to the latitude specified in config file (Luxor = 25.6989 degrees). 
    Returns:
        (RA, Dec): a tuple containing the right ascension and declination of the star in hours and degrees, respectively. 
    '''
    # minimum visible Dec 
    maxnum = 0.5 * (np.cos(np.pi * lat / 180) + 1) # for generating stars visible at Luxor below the hemisphere
    # select a random point uniformly across the visible hemisphere
    x = random.random()
    y = random.uniform(0, maxnum) # 0-0.5: above horizon, 0.5-1 below horizon
    phi = np.arccos(2*y-1) - np.pi/2
    theta = 2 * np.pi * x
    # transform coordinates into RA and Dec
    RA = Angle(theta * u.radian).hour
    Dec = Angle(phi * u.radian).deg
    return RA, Dec

def generate_star_field(num):
    '''
    Generate a number = num of stars with random Dec and RA values.
    It is expected that num < 1000, otherwise the naming convention will be wrong. 
    To add more, change "{:04.0f}" to "{:0X.0f}", where X is the number of digits. 
    '''
    star_names = [f"R{i:04.0f}" for i in range(num)]
    ra_dec_list = [generate_star() for _ in range(num)]
    ra_list, dec_list = zip(*ra_dec_list)
    mag_list = np.round(generate_mags_from_Hdist(num), 2)
    hd_list = (
        ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
        + [f"{name} Azimuth" for name in star_names]
        + [f"{name} Altitude" for name in star_names]
    )
    obj_list = SkyCoord(list(ra_list) * u.hour, list(dec_list) * u.deg)
    df = pd.DataFrame({
        "Name": star_names,
        "RA": ra_list,
        "Dec": dec_list,
        "Mag": mag_list
    })
    return obj_list, hd_list, df

def generate_star_belt(star, num, year, dec_off: 0.0):
    ''' 
    Generate a belt of stars around a reference star (e.g., Sirius), 
    precessed to a given year and offset in right ascension.

    Args:
        star: Reference SkyCoord object (e.g., Sirius).
        num: Number of stars to generate.
        year: Year (BC) for precession.
        dec_off: Declination offset in degrees [optional].

    Returns:
        obj_list: List of SkyCoord objects for the belt stars.
        star_names: List of star names.
        hd_list: List of header strings for output tables.
    '''
    # Precess reference star to given year
    precessed_objs, _ = precessed_coords([star], year)
    ref_obj = precessed_objs[0]
    ra0 = Angle(ref_obj.ra, unit="deg").hour
    dec = Angle(ref_obj.dec, unit="deg") + Angle(dec_off, unit="deg")

    obj_list = []
    star_names = []
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]

    for i in range(num):
        name = f"S{i:02d}"
        star_names.append(name)
        ra = (ra0 + 360 / num * i) % 360
        obj = SkyCoord(ra=ra, dec=dec, unit="deg")
        obj_list.append(obj)
        hd_list.append(f"{name} Azimuth")
        hd_list.append(f"{name} Altitude")
    # TO DO: return star names!
    return obj_list, hd_list

# Helper functions for generate_mags_from_Hdist and generate_rand_from_pdf

def kde(x, x_grid, bandwidth=0.2, **kwargs):
    '''
    Kernel Density Estimation with Scipy.
    Adapted from: https://stackoverflow.com/questions/17821458/random-number-from-histogram
    Inputs:
        x: data points to estimate the density from
        x_grid: points at which to evaluate the density
        bandwidth: bandwidth for the Gaussian kernel (default is 0.2)
        **kwargs: additional keyword arguments for gaussian_kde
    Returns:
        kde.evaluate(x_grid): estimated density at the points in x_grid 
    '''
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

def generate_rand_from_pdf(pdf, x_grid, num):
    '''
    Generate random samples from a probability density function (PDF) using the cumulative distribution function (CDF).
    Adapted from: https://stackoverflow.com/questions/17821458/random-number-from-histogram 
    Inputs:
        pdf: probability density function values at x_grid points
        x_grid: points at which the PDF is defined
        num: number of random samples to generate       
    Returns:
        random_from_cdf: an array of num random samples drawn from the PDF  
    '''
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(num)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf
    




# def mockCoords_randomStar():
#     '''
#     Make a single random star in Dec and RA.
#     '''
#     #minimum dec visible 
#     lat = 25.6989 # lat of Luxor in degrees
#     maxnum = 0.5 * (np.cos(np.pi * lat / 180) + 1) # for generating stars visible at Luxor below the hemisphere
#     #select a random point uniformly across a hemisphere +
#     x = random.random()
#     y = random.uniform(0, maxnum) #0-0.5: above horizon, 0.5-1 below horizon
#     phi = np.arccos(2*y-1) - np.pi/2
#     theta = 2 * np.pi * x
#     #mag = random.uniform(1,6) # magnitude of star (if we want it)
#     # make them into RA and Dec
#     RA = Angle(theta * u.radian).hour
#     Dec = Angle(phi * u.radian).deg
#     # obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
#     return (RA, Dec)

# def mockCoords_randomStarField(num):
#     '''
#     Make a number = num of random stars in Dec and RA to analyze with decanOpy.
#     It is expected that num < 100, otherwise the naming convention will be wrong. 
#     To add more, change "{:02.0f}" to "{:0X.0f}", where X is the number of digits. 
#     '''
#     # initalize empty key lists
#     star_names = []
#     obj_list = []
#     RA_list = []
#     Dec_list = []
#     mag_list = np.round(generate_mags_from_Hdist(num), 2)
#     hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
#     for i in range(0, num):
#         name = "R" + "{:04.0f}".format(i)
#         star_names.append(name)
#         (RA, Dec) = mockCoords_randomStar()
#         RA_list.append(RA)
#         Dec_list.append(Dec)
#         #mag_list.append(round(random.uniform(-1.5, 6), 2))
#         #obj = randomStar()
#         #obj_list.append(obj)
#         hd_list.append(name + " Azimuth")
#         hd_list.append(name + " Altitude")
#     obj_list = SkyCoord(RA_list * u.hour, Dec_list * u.deg)
#     df = pd.DataFrame({"Name" : star_names, 
#                 "RA" : RA_list,
#                 "Dec": Dec_list,
#                 "Mag": mag_list})  
#     return(obj_list, hd_list, df)  

# def mockCoords_StarLike(star, num, year, dec_off):
#     '''
#     A function to create Dec and RA structures of fake stars for testing Sirius-like behavior. 
#     The input "num" refers to the number of stars created and must be an integer. 
#     The input "year" designates the year BC for performing precession on Sirius coords. 
#     The input "dec_off" is an optional offset in Dec for Sirius, assumed to be in degrees.
#     '''
#     # initalize empty key lists
#     star_names = []
#     obj_list = []
#     hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
#     # Sirius Dec
#     (obj_list, hd_list) = precessedCoords([star], year) # get Sirius data for given year BC
#     obj = obj_list[0] #extract Sirius data from list strucure
#     RA0 = Angle(obj.ra, unit="deg").hour
#     Dec = Angle(obj.dec, unit="deg") + Angle(dec_off, unit = "deg")
#     #Dec = Angle(-17.849335700373032, unit="deg").deg
#     # populate both
#     for i in range(0, num):
#         name = "S" + "{:02.0f}".format(i)
#         star_names.append(name)
#         # Sirius RA
#         RA = (RA0 + 360/num * i) % 360 # stepping by 1 hr = 15 deg mod 360
#         obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
#         obj_list.append(obj)
#         hd_list.append(name + " Azimuth")
#         hd_list.append(name + " Altitude")
#     return(obj_list, hd_list)
