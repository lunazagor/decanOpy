### Import Statements
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, Angle, Longitude
from sunpy.coordinates import frames, sun
import star_chart_spherical_projection as scsp
import csv
import argparse
import time, os, fnmatch, shutil
import random
import pandas as pd
import errno
from decanO import *

# supress dubious year warnings 
import warnings
warnings.simplefilter('ignore', UserWarning)
 

# ##
# #### Parse the arguments
# ##
parser = argparse.ArgumentParser('Restart Unfinished Run from Last Save')
parser.add_argument('-p', '--path', required=True) 
parser.add_argument('-p2ICs', '--path2ICs', required=True) 
parser.add_argument('-matchS', '--matchStellariumJD', required=False, default=True) #JD offset to match stellarium


args = parser.parse_args()

matchStellariumJD = bool(args.matchStellariumJD)
path = str(args.path)
p2ICs = str(args.path2ICs)

# check if paths to files exist
direct = os.getcwd() # current working directory
abspath = direct +  path
absp2ICs = direct + p2ICs


if not os.path.exists(abspath):
    error = "ERROR: file does not exist at " + abspath
    raise Exception(error)

if not os.path.exists(absp2ICs):
    error = "ERROR: initial conditions file does not exist at " + absp2ICs
    raise Exception(error)    

# read files & check that number of stars lines up
df = pd.read_csv(abspath, sep='|')
df2 = pd.read_csv(absp2ICs, sep=',')

(trash, nhead) = df.shape
(nstar, trash) = df2.shape

if nhead != 2 * nstar + 4:
    raise Exception("Initial conditions and aborted run do not match in number of stars!")

# get last saved iteration
dat = df.iloc[-1]["Local Date and Time"]
year = dat[2:6] # assuming 4 digits in year
jdinit = df.iloc[0]["Julian Date"]
jdsave = df.iloc[-1]["Julian Date"]

# determine the last saved hour and minute
# savhour = int(dat[13:15])
# savmin = int(dat[16:18])//4

# # Set Location on Earth (currently hardcoded to Luxor, Egypt)
Luxor = EarthLocation(lat=25.6989*u.deg, lon=32.6421*u.deg, height=89*u.m) # data matched to Stellarium

# # Times and dates
# hour and minute steps
dhour = 0.04166666674427688 #iterate every hour
d4min = 0.00277777784503996 # iterate every 4 minutes


# start time to jd with offset for local sidereal time
start = jdinit #(Time('-0' + year + '-' + month + '-01T00:00:00.000', scale="local", location = Luxor).jd) - (Luxor.lon.deg/15.0) * dhour + dS

# values to iterate over
days = start + 1 * np.arange(0, 365) # iterate for a year 
hours = dhour * np.arange(0, 24)
minutes = d4min * np.arange(0, 15)

endtime = days[-1] + hours[-1] + minutes[-1] 
numiter = round((endtime - jdsave) / d4min) #number of iterations left to finish the year 

# introduce offset of a number of days to get JD in line with Stellarium:
dS = 0
if matchStellariumJD:
    dS = dS_offset(year)

## Get RA/Dec of decans while accounting for precession of the equinoxes
# NOTE: scsp uses the Vondrak precession algorithm which doesn't EXACTLY match Stellarium, so some differences are to be expected!
# Other algorithms may be introduced in the future

# get headers 
hd_list = list(df)

# get object coords 
RA_list = list(df2['RA'])
Dec_list = list(df2['Dec'])
obj_list = SkyCoord(RA_list * u.hour, Dec_list * u.deg)

###
##### Writing the .txt file
###

with open(abspath, "a", newline='') as file:
    writer = csv.writer(file, delimiter='|')
    for i in range (0, numiter):
        temptime = jdsave + d4min * i + d4min
        # Sun coords
        c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=Time(temptime, format = 'jd'), observer="earth", frame=frames.Helioprojective)
        frame_altaz = AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor)
        sun_altaz = c.transform_to(frame_altaz)
        # decan coords
        info = ['{0:.16f}'.format(np.round(temptime, 10)),
                            str(Time(temptime - dS + (Luxor.lon.deg/15.0) * dhour, format = 'jd').fits), # local time and date
                            '{0:.3f}'.format(sun_altaz.T.az)[0:-4], # [0:-4] get rid of trailing " deg" for later analysis
                            '{0:.3f}'.format(sun_altaz.T.alt)[0:-4]]
        for obj in obj_list:
            (alt, az) = calc_altaz(Angle(obj.ra, unit="deg").hour, obj.dec, Luxor, temptime)
            info.append('{0:.3f}'.format(az))
            info.append('{0:.3f}'.format(alt))
        # Write row to file
        writer.writerow(info)


