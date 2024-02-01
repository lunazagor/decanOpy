### Import Statements
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, Angle, Longitude
from sunpy.coordinates import frames, sun
import star_chart_spherical_projection as scsp
import csv
import argparse
import os
import errno
from decanO import *

# supress dubious year warnings 
import warnings
warnings.simplefilter('ignore', UserWarning)


# ##
# #### Parse the arguments
# ##
parser = argparse.ArgumentParser('My program')
parser.add_argument('-d','--decan', nargs='+', type=str, required=True)
parser.add_argument('-yBC', '--yearBC', required=True)
parser.add_argument('-m', '--month', required=False, default = "01")
parser.add_argument('-matchS', '--matchStellariumJD', required=False, default=True)
parser.add_argument('-n', '--name', required=False, default="data")

args = parser.parse_args()

decans = list(args.decan)
year = str(args.yearBC)
month = str(args.month)
matchStellariumJD = bool(args.matchStellariumJD)
name = str(args.name)


# ##
# #### Determining the when and where
# ##

# # Set Location on Earth (currently hardcoded to Luxor, Egypt)
Luxor = EarthLocation(lat=25.6989*u.deg, lon=32.6421*u.deg, height=89*u.m) # data matched to Stellarium

# # Times and dates
# hour and minute steps
dhour = 0.04166666674427688 #iterate every hour
d4min = 0.00277777784503996 # iterate every 4 minutes

# introduce offset of a number of days to get JD in line with Stellarium:
dS = 0
if matchStellariumJD:
    dS = dS_offset(year)

# star time to jd with offset for local sidereal time
start = (Time('-0' + year + '-' + month + '-01T00:00:00.000', scale="local", location = Luxor).jd) - (Luxor.lon.deg/15.0) * dhour + dS

# values to iterate over
#days = start + 1 * np.arange(0, 365) # iterate for a year 
days = start + 1 * np.arange(0, 1) # iterate for a day (for debugging)
hours = dhour * np.arange(0, 24)
minutes = d4min * np.arange(0, 15)

## Get RA/Dec of decans while accounting for precession of the equinoxes
# NOTE: scsp uses the Vondrak precession algorithm which doesn't EXACTLY match Stellarium, so some differences are to be expected!
# Other algorithms may be introduced in the future
(obj_list, hd_list) = precessedCoords(decans, year)

###
##### Writing the .txt file
###
direct = os.getcwd() # current working directory
direct = direct + '/DecanLists' # directory where the .txt files go

if not os.path.exists(os.path.dirname(direct)):
    try:
        os.makedirs(os.path.dirname(direct))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

# name the txt file (default is "data")
filename = direct + "/" + name + year + "BC.txt"


# start writing the file

with open(filename, "w", newline='') as file:
    writer = csv.writer(file, delimiter='|')
    writer.writerow(hd_list) # write headers
    for day in days:
        for hour in hours:
            for mins in minutes:
                temptime = day + hour + mins
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
                    # info_temp = obj.transform_to(AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor))
                    # info.append('{0.az:.1f}'.format(info_temp)[0:-4]) 
                    # info.append('{0.alt:.1f}'.format(info_temp)[0:-4]) 
                # Write to file
                writer.writerow(info)