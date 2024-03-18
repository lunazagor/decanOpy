### Import Statements

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

def precessedCoords(declist, year):
    '''
    A function to find the Dec and RA of list of decans accounting for the precession of the equinoxes.
    Inputs: 
        declist = list of named decans to querry
        year = 
    Outputs:
        obj_list = list of decans' RAs and Decs as an astropy SkyCoords object
        year = year BC for which to return data
        hd_list = a variable to set the header for the csv file
    NOTE: I should change this so it uses SkyCoord(List) instead of List(SkyCoord--it'll be faster)
    '''
    # querry available stars from scsp
    # note: need to fix this for stars not listed by name in scsp!
    years_since = -2000 - int(year)
    star_dict = scsp.finalPositionOfStars(declist, yearSince2000=years_since)
    star_names = star_dict.keys()

    # set up lists 
    obj_list = []
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]

    # calculate objects
    for name in star_names:
        ra_temp = star_dict[name]["RA"].split(".")
        RA = Angle(ra_temp[0] + "h" + ra_temp[1] + "m" + ra_temp[2] + "s").deg
        dec_temp = star_dict[name]["Declination"]
        Dec = Angle(dec_temp, unit="deg").deg
        obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
        obj_list.append(obj)
        hd_list.append(name + " Azimuth")
        hd_list.append(name + " Altitude")
    return(obj_list, hd_list)

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

    return(alt, az)


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



# Import Decan Data

def ImportDecanData(direct, filename):
    
    '''
    A function to import data from a decanOpy-generated .txt file. 
    Inputs: 
        direct = string with the directory where the .txt file is located
        filename = string with name of file (name + month + year)
    Outputs:
        jd = Julian date
        date = human readable date
        DecAz = the azimuth of the decan
        DecAlt = the altitude of the decan
        SunAz = the azimuth of the Sun
        SunAlt = the altitude of the Sun
    '''
    
    jd = []
    date = []
    DecAz = []
    DecAlt = []
    SunAz = []
    SunAlt = []
    # Import Single Object
    with open(direct + filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='|')
        decan = next(csv_reader)[0]
        location = next(csv_reader)[0]
        trash = next(csv_reader)
        headers = next(csv_reader)
        for row in csv_reader:
            # time info
            jd.append(float(row[0]))
            date.append(row[1])
            # decan info
            DecAz.append(float(row[2][0:-4]))
            DecAlt.append(float(row[3][0:-4]))
            # solar info
            SunAz.append(float(row[4][0:-4]))
            SunAlt.append(float(row[5][0:-4]))
    return(jd, date, DecAz, DecAlt, SunAz, SunAlt)

def JustDecanData(direct, filename):
    
    '''
    A function to import just the can data data from a decanOpy-generated .txt file. 
    Used for the MaxMinAltAz function.
    Inputs: 
        direct = string with the directory where the .txt file is located
        filename = string with name of file (name + month + year)
    Outputs:
        DecAz = the azimuth of the decan
        DecAlt = the altitude of the decan
    '''
    
    DecAz = []
    DecAlt = []
    # Import Single Object
    with open(direct + filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='|')
        decan = next(csv_reader)[0]
        location = next(csv_reader)[0]
        trash = next(csv_reader)
        headers = next(csv_reader)
        for row in csv_reader:
            # decan info
            DecAz.append(float(row[2][0:-4]))
            DecAlt.append(float(row[3][0:-4]))
            # solar info
    return(DecAz, DecAlt)

# Get Sunset and Sunrise Times

def SunRiseSet(jd, SunAlt):
    
    '''
    A function to create a list of indices where the Sun rises and sets in a given year. 
    This is useful for making sure we're tracking nightly, visible motion of the decans.
    Inputs: 
        jd = Julian date
        SunAlt = the altitude of the Sun
    Outputs:
        sunriseset = indices of sunrize and sunset in the jd & date columns
    '''
    
    sunriseset = []
    for i in range(360, len(jd), 360):
        temp = []
        for j in range(i - 360, i):
            if SunAlt[j] <= 0.4 and SunAlt[j] >= -0.4:
                if len(temp) == 0: 
                    temp.append(j)
                elif temp[-1] != j - 1:
                    temp.append(j)
        sunriseset.append(temp)
    return sunriseset

# Maximum and Minimum Nightly Altitude of Object

def MaxMinAltAz(direct, filename, jd, sunriseset):
    
    '''
    A function to create lists of minimum and maximum azimuths and altitudes of the decan. 
    This is useful for making sure we're tracking nightly, visible motion of the decans.
    Inputs: 
        direct = string with the directory where the .txt file is located
        filename = string with name of file (name + month + year)
        jv = Julian date
        sunriseset = indices of sunrize and sunset in the jd & date columns
    Outputs:
        sunriseset = indices of sunrize and sunset in the jd & date columns
        days = list of indices when it's daylight 
        minaz, maxaz = minimum and maximum azimuths of the decan per night
        minalt, maxalt = minimum and maximum altitudes of the decan per night
        riseaz, setaz = azimuth of decan at rise & set
        risealt, setalt = altitude of decan at rise & set
    '''
    
    (DecAz, DecAlt) = JustDecanData(direct, filename)
    maxalt = []
    minalt = []
    maxaz = []
    minaz = []
    riseaz = []
    setaz = []
    risealt = []
    setalt = []
    days = []
    for i in range(0, int(len(jd)/360) - 2):
        sset = sunriseset[i][1]
        srise = sunriseset[i + 1][0]
        maxalt.append(max(DecAlt[sset:srise]))
        minalt.append(min(DecAlt[sset:srise]))
        maxaz.append(max(DecAz[sset:srise]))
        minaz.append(min(DecAz[sset:srise]))
        riseaz.append(DecAz[srise])
        setaz.append(DecAz[sset])
        risealt.append(DecAlt[srise])
        setalt.append(DecAlt[sset])
        days.append(DecAlt[srise:sset])
    return(days, minaz, maxaz, minalt, maxalt, riseaz, setaz, risealt, setalt)


def kde(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)


def generate_rand_from_pdf(pdf, x_grid, num):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(num)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf
    

def generate_mags_from_Hdist(num):
    # get data 
    direct = os.getcwd()
    name_df = pd.read_csv(direct + '/StarLists/RealSky/ICs/star_data_names.csv', index_col=None, header=0, names=['Name', 'RA', 'Dec', 'Mag'])
    data = name_df['Mag'].values #np.random.normal(size=1000)
    # create kernel density estimator pdf
    x_grid = np.linspace(min(data), max(data), num)
    kdepdf = kde(data, x_grid, bandwidth=0.1)
    random_from_kde = generate_rand_from_pdf(kdepdf, x_grid, num)
    # return an array of num magnitudes
    return random_from_kde

def mockCoords_randomStar():
    '''
    Make a single random star in Dec and RA.
    '''
    #minimum dec visible 
    lat = 25.6989 # lat of Luxor in degrees
    maxnum = 0.5 * (np.cos(np.pi * lat / 180) + 1) # for generating stars visible at Luxor below the hemisphere
    #select a random point uniformly across a hemisphere +
    x = random.random()
    y = random.uniform(0, maxnum) #0-0.5: above horizon, 0.5-1 below horizon
    phi = np.arccos(2*y-1) - np.pi/2
    theta = 2 * np.pi * x
    #mag = random.uniform(1,6) # magnitude of star (if we want it)
    # make them into RA and Dec
    RA = Angle(theta * u.radian).hour
    Dec = Angle(phi * u.radian).deg
    # obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
    return (RA, Dec)

def mockCoords_randomStarField(num):
    '''
    Make a number = num of random stars in Dec and RA to analyze with decanOpy.
    It is expected that num < 100, otherwise the naming convention will be wrong. 
    To add more, change "{:02.0f}" to "{:0X.0f}", where X is the number of digits. 
    '''
    # initalize empty key lists
    star_names = []
    obj_list = []
    RA_list = []
    Dec_list = []
    mag_list = np.round(generate_mags_from_Hdist(num), 2)
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
    for i in range(0, num):
        name = "R" + "{:04.0f}".format(i)
        star_names.append(name)
        (RA, Dec) = mockCoords_randomStar()
        RA_list.append(RA)
        Dec_list.append(Dec)
        #mag_list.append(round(random.uniform(-1.5, 6), 2))
        #obj = randomStar()
        #obj_list.append(obj)
        hd_list.append(name + " Azimuth")
        hd_list.append(name + " Altitude")
    obj_list = SkyCoord(RA_list * u.hour, Dec_list * u.deg)
    df = pd.DataFrame({"Name" : star_names, 
                "RA" : RA_list,
                "Dec": Dec_list,
                "Mag": mag_list})  
    return(obj_list, hd_list, df)  

def mockCoords_StarLike(star, num, year, dec_off):
    '''
    A function to create Dec and RA structures of fake stars for testing Sirius-like behavior. 
    The input "num" refers to the number of stars created and must be an integer. 
    The input "year" designates the year BC for performing precession on Sirius coords. 
    The input "dec_off" is an optional offset in Dec for Sirius, assumed to be in degrees.
    '''
    # initalize empty key lists
    star_names = []
    obj_list = []
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
    # Sirius Dec
    (obj_list, hd_list) = precessedCoords([star], year) # get Sirius data for given year BC
    obj = obj_list[0] #extract Sirius data from list strucure
    RA0 = Angle(obj.ra, unit="deg").hour
    Dec = Angle(obj.dec, unit="deg") + Angle(dec_off, unit = "deg")
    #Dec = Angle(-17.849335700373032, unit="deg").deg
    # populate both
    for i in range(0, num):
        name = "S" + "{:02.0f}".format(i)
        star_names.append(name)
        # Sirius RA
        RA = (RA0 + 360/num * i) % 360 # stepping by 1 hr = 15 deg mod 360
        obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
        obj_list.append(obj)
        hd_list.append(name + " Azimuth")
        hd_list.append(name + " Altitude")
    return(obj_list, hd_list)

def randomStar():
    '''
    Make a single random star in Dec and RA.
    '''
    #select a random point uniformly across a hemisphere
    x = random.random()
    y = random.uniform(0, 0.5) #0-0.5: above horizon, 0.5-1 below horizon
    phi = np.arccos(2*y-1) - np.pi/2
    theta = 2 * np.pi * x
    #mag = random.uniform(1,6) # magnitude of star (if we want it)
    # make them into RA and Dec
    RA = Angle(theta * u.radian).hour
    Dec = Angle(phi * u.radian).deg
    obj = SkyCoord(ra=RA, dec=Dec, unit="deg")
    return obj

def randomStarField(num):
    '''
    Make a number = num of random stars in Dec and RA to analyze with decanOpy.
    It is expected that num < 100, otherwise the naming convention will be wrong. 
    To add more, change "{:02.0f}" to "{:0X.0f}", where X is the number of digits. 
    '''
    # initalize empty key lists
    star_names = []
    obj_list = []
    hd_list = ["Julian Date", "Local Date and Time", "Sun Azimuth", "Sun Altitude"]
    for i in range(0, num):
        name = "R" + "{:02.0f}".format(i)
        star_names.append(name)
        obj = randomStar()
        obj_list.append(obj)
        hd_list.append(name + " Azimuth")
        hd_list.append(name + " Altitude")
    return(obj_list, hd_list)    


# ##
# #### Parse the arguments
# ##


# parser = argparse.ArgumentParser('My program')
# parser.add_argument('-d','--decan', nargs='+', type=str, required=True)
# parser.add_argument('-yBC', '--yearBC', required=True)
# parser.add_argument('-m', '--month', required=False, default = "01")
# parser.add_argument('-mS', '--matchStellariumJD', required=False, default=True)

# args = parser.parse_args()

# decans = list(args.decan)
# year = str(args.yearBC)
# month = str(args.month)
# matchStellariumJD = bool(args.matchStellariumJD)


# # ##
# # #### Determining the when and where
# # ##

# # # Set Location on Earth 

# Luxor = EarthLocation(lat=25.6989*u.deg, lon=32.6421*u.deg, height=89*u.m) # data matched to Stellarium

# # # Length of month 

# len31 = ['01', '03', '05', '07', '08', '10', '12']
# len30 = ['02', '04', '06', '09', '11']

# if month in len31:
# 	end = 31
# else:
# 	end = 30

# # # Times and dates

# # hour and minute steps
# dhour = 0.04166666674427688
# d4min = 0.00277777784503996


# # introduce offset of a few days to get jd in line with Stellarium:
# dS = 0
# if matchStellariumJD:
#     dS = dS_offset(year)

# # star time to jd with offset for local sidereal time
# start = (Time('-0' + year + '-' + month + '-01T00:00:00.000', scale="local", location = Luxor).jd) - (Luxor.lon.deg/15.0) * dhour + dS

# # values to iterate over
# #days = start + 1 * np.arange(0, 365) # iterate for a year 
# days = start + 1 * np.arange(0, 1) # iterate for a day to debug 
# hours = dhour * np.arange(0, 24)
# minutes = d4min * np.arange(0, 15)

# ###
# ##### Writing the .txt file
# ###
# direct = os.getcwd() # current working directory
# direct = direct + '/DecanLists' # directory where the .txt files go

# if not os.path.exists(os.path.dirname(direct)):
#     try:
#         os.makedirs(os.path.dirname(direct))
#     except OSError as exc: # Guard against race condition
#         if exc.errno != errno.EEXIST:
#             raise

# # name the txt file
# # note to self: should I just make it a csv or something?

# filename = direct + "/" + "testdata" + month + year + "BC.txt"

# ## Get coordinates of decans while accounting for precession of the equinoxes

# #obj = SkyCoord.from_name(decan)
# (obj_list, hd_list) = precessedCoords(decans, year)


# # start writing the file

# with open(filename, "w", newline='') as file:
#     writer = csv.writer(file, delimiter='|')
#     #writer.writerow(["Object: " + str(decans)])
#     #writer.writerow(["Luxor = EarthLocation(lat=25.68*u.deg, lon=31.55*u.deg, height=89*u.m)"])
#     #writer.writerow(["\n"])
#     writer.writerow(hd_list) # write headers
#     for day in days:
#         for hour in hours:
#             for mins in minutes:
#                 temptime = day + hour + mins
#                 # Sun coords
#                 c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=Time(temptime, format = 'jd'), observer="earth", frame=frames.Helioprojective)
#                 frame_altaz = AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor)
#                 sun_altaz = c.transform_to(frame_altaz)
#                 # decan coords
#                 info = ['{0:.16f}'.format(np.round(temptime, 10)),
#                                  str(Time(temptime - dS + (Luxor.lon.deg/15.0) * dhour, format = 'jd').fits), # local time and date
#                                  '{0:.1f}'.format(sun_altaz.T.az)[0:-4], # [0:-4] get rid of trailing " deg" for later analysis
#                                  '{0:.1f}'.format(sun_altaz.T.alt)[0:-4]]
#                 for obj in obj_list:
#                     (alt, az) = calc_altaz(Angle(obj.ra, unit="deg").hour, obj.dec, Luxor, temptime)
#                     info.append('{0:.3f}'.format(az))
#                     info.append('{0:.3f}'.format(alt))
#                     # info_temp = obj.transform_to(AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor))
#                     # info.append('{0.az:.1f}'.format(info_temp)[0:-4]) 
#                     # info.append('{0.alt:.1f}'.format(info_temp)[0:-4]) 
#                 # Write to file
#                 writer.writerow(info)



# # ##
# # #### Functions for data processing
# # ##

