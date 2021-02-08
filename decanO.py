### Import Statements
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from sunpy.coordinates import frames, sun
import csv
import argparse

parser = argparse.ArgumentParser('My program')
parser.add_argument('-decan', '--one')
parser.add_argument('-yearBC', '--two')
parser.add_argument('-month', '--three')

args = parser.parse_args()


decan = str(args.one)
year = str(args.two)
month = str(args.three)

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


with open(decan + month + year + "BC.txt", "w", newline='') as file:
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


                        
