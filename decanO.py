### Import Statements
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from sunpy.coordinates import frames, sun
import csv
import argparse
import os
import errno


###
##### Parse the arguments
###


parser = argparse.ArgumentParser('My program')
parser.add_argument('-decan', '--one')
parser.add_argument('-yearBC', '--two')
parser.add_argument('-month', '--three')

args = parser.parse_args()


decan = str(args.one)
year = str(args.two)
month = str(args.three)

###
##### Determining the when and where
###

## Set Location on Earth 

Luxor = EarthLocation(lat=25.68*u.deg, lon=31.55*u.deg, height=89*u.m)

## Get coordinates of Object

obj = SkyCoord.from_name(decan)

## Length of month 

len31 = ['01', '03', '05', '07', '08', '10', '12']
len30 = ['02', '04', '06', '09', '11']

if month in len31:
	end = 31
else:
	end = 30

## Times and dates

start = (Time('-0' + year + '-' + month + '-01T00:00:00.000').jd)
days = start + 1 * np.arange(0, 365) # iterate for a year 
dhour = 0.04166666674427688
d4min = 0.00277777784503996
hours = dhour * np.arange(0, 24)
minutes = d4min * np.arange(0, 15)


###
##### Writing the .txt file
###
direct = os.getcwd() # current working directory
ditect = direct + '/DecanLists' # directory where the .txt files go

if not os.path.exists(os.path.dirname(direct)):
    try:
        os.makedirs(os.path.dirname(direct))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise


filename = direct + "/" + decan + month + year + "BC.txt"

with open(filename, "w", newline='') as file:
    writer = csv.writer(file, delimiter='|')
    writer.writerow(["Object: " + decan])
    writer.writerow(["Luxor = EarthLocation(lat=25.68*u.deg, lon=31.55*u.deg, height=89*u.m)"])
    writer.writerow(["\n"])
    writer.writerow(["Julian Date", " Human Readable Date",
                     decan + " Azimuth", decan + " Altitude", "Sun Azimuth", "Sun Altitude"])
    for day in days:
        for hour in hours:
            for mins in minutes:
                temptime = day + hour + mins
                # decan coords
                info = obj.transform_to(AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor))
                # Sun coords
                c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=Time(temptime, format = 'jd'), observer="earth", frame=frames.Helioprojective)
                frame_altaz = AltAz(obstime=Time(temptime, format = 'jd'), location=Luxor)
                sun_altaz = c.transform_to(frame_altaz)
                # Write to file
                writer.writerow(['{0:.16f}'.format(np.round(temptime, 10)),
                                 str(Time(temptime, format = 'jd').fits),
                                 '{0.az:.3}'.format(info),
                                 '{0.alt:.3}'.format(info),
                                 '{0:.1f}'.format(sun_altaz.T.az),
                                 '{0:.1f}'.format(sun_altaz.T.alt)])



###
##### Functions for data processing
###


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
                        
