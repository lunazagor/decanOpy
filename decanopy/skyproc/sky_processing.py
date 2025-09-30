# functions for processing sky data from skyflow 

# NOTE: this is the some of the oldest decanOpy code, could use some refactoring
# may particularly benefit from renaming!
# pivot to object-oriented design?

import numpy as np
from pathlib import Path
import pandas as pd
from collections import Counter


def StarRiseSet(jd, starAlt, deg):
    # create return arrays
    star_rise = np.zeros(int(len(jd)/360), dtype=int)
    star_set = np.zeros(int(len(jd)/360), dtype=int)
    # loop through each day
    for i in range(360, len(jd) + 1, 360):
        # isolate data for one day
        dailyStarAlt = starAlt[i - 360 : i]
        # is the star visible at the beginning of the day
        starVis = dailyStarAlt[0] > deg
        # make array of booleans
        bool_arr = dailyStarAlt > deg   
        # find two indices of change
        ind1 = np.argwhere(bool_arr != starVis)[0][0]
        if len(np.argwhere(bool_arr[ind1:]== starVis)) != 0: # edge case of next index is in next day
            ind2 = np.argwhere(bool_arr[ind1:]== starVis)[0][0] + ind1
        else:
            ind2 = 360    
        # assign indices
        if starVis:
            star_set[int((i-360)/360)] = int(-360 + i  + ind1)
            star_rise[int((i-360)/360)] = int(-360 + i  + ind2)
        else:
            star_set[int((i-360)/360)] = int(-360 + i  + ind2)
            star_rise[int((i-360)/360)] = int(-360 + i  + ind1)
    return (star_rise, star_set)

def SunRiseSet(jd, SunAlt, deg):
    
    '''
    A function to create a list of indices where the Sun rises and sets in a given year.
    This is generalized to find the sun above/below soem altitude given in degrees deg.  
    This is useful for making sure we're tracking nightly, visible motion of the decans.
    NOTE: as written, this code assumes that data is collected every 4 minutes. 
    To change this, change number to number of collection intervals per day! (360 = 24 * 60/4)
    Inputs: 
        jd = Julian date
        SunAlt = the altitude of the Sun
    Outputs:
        sunriseset = indices of sunrise and sunset in the jd & date columns
    '''
    
    sunriseset = []
    for i in range(360, len(jd), 360):
        temp = []
        for j in range(i - 360, i):
            if SunAlt[j] <= deg + 0.4 and SunAlt[j] >= deg - 0.4:
                if len(temp) == 0: 
                    temp.append(j)
                elif temp[-1] != j - 1:
                    temp.append(j)
        sunriseset.append(temp)
    return sunriseset

def isStarVisible(sunSet, sunRise, starAlt):
    vis_arr = np.full(365, True)
    max_alt_arr = np.zeros(364)
    for i in range(0, 364):
        maxalt = max(starAlt[sunSet[i]:sunRise[i + 1]])
        max_alt_arr[i] = maxalt
        #print(maxalt)
        if maxalt < 0:
            vis_arr[i] = False
    return(max_alt_arr, vis_arr)

def MaxMinAltAz(jd, sunriseset, DecAz, DecAlt):
    
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
    maxalt = []
    minalt = []
    maxaz = []
    minaz = []
    riseaz = []
    setaz = []
    risealt = []
    setalt = []
    days = []
    for i in range(0, int(len(jd)/360) - 1):
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

def initialize_sky1(
    filepath,
    star_rise_deg: float = 10,
    sun_rise_deg: float = -12
):
    """
    Initialize the sky model by reading star position data and calculating rise/set times and visibility.

    Args:
        filepath (str or Path): Path to the star position data file.
        star_rise_deg (float): Altitude (deg) above horizon to define star rise. The default is 10.
        sun_rise_deg (float): Sun altitude (deg) to define sunrise/set. The default is -12 (nautical twilight). 

    Returns:
        skydict: Dictionary containing all relevant arrays and lists.
    """

    # Read data
    decan_output = pd.read_csv(filepath, sep="|")
    header = decan_output.columns

    # Star names
    starlist = [name[0:-8] for name in header[4:-1:2]]

    # Standard data
    jd = decan_output[header[0]].to_numpy()
    hrd = decan_output[header[1]]
    sunAz = decan_output[header[2]].to_numpy()
    sunAlt = decan_output[header[3]].to_numpy()

    # Star data
    num_decs = (len(header) - 4) // 2
    starsAz = np.stack([decan_output[header[4 + 2 * i]].to_numpy() for i in range(num_decs)])
    starsAlt = np.stack([decan_output[header[5 + 2 * i]].to_numpy() for i in range(num_decs)])

    # Sunrise and sunset times
    sunRise, sunSet = StarRiseSet(jd, sunAlt, sun_rise_deg)
    sunAzSet = sunAz[sunSet]

    # Star rise/set/visibility/max altitude
    starAzRiseList = np.zeros((num_decs, len(sunRise)))
    starVisList = np.full((num_decs, len(sunRise)), True)
    starMaxAltList = np.zeros((num_decs, len(sunRise)-1))

    for i in range(num_decs):
        min_alt = np.min(starsAlt[i])
        max_alt = np.max(starsAlt[i])
        if min_alt >= star_rise_deg:
            # Circumpolar
            starVisList[i, :] = True
            starMaxAltList[i, :] = max_alt
        elif max_alt < star_rise_deg:
            # Never rises
            starVisList[i, :] = False
            starMaxAltList[i, :] = max_alt
        else:
            # Sometimes visible
            starRise, starSet = StarRiseSet(jd, starsAlt[i], star_rise_deg)
            starAzRise = starsAz[i, starRise]
            starAzRiseList[i, :] = starAzRise
            maxAlt, starVis = isStarVisible(sunSet, sunRise, starsAlt[i])
            starVisList[i, :] = starVis
            starMaxAltList[i, :] = maxAlt

    skydict =  {
        "jd": jd,
        "hrd": hrd,
        "sunAz": sunAz,
        "sunAlt": sunAlt,
        "starlist": starlist,
        "starsAz": starsAz,
        "starsAlt": starsAlt,
        "sunRise": sunRise,
        "sunSet": sunSet,
        "sunAzSet": sunAzSet,
        "starAzRiseList": starAzRiseList,
        "starVisList": starVisList,
        "starMaxAltList": starMaxAltList
    }

    return skydict

def initialize_sky(
        skyname: str,     
        skytype: str, 
        star_rise_deg: float = 10,
        sun_rise_deg: float = -12):
    """
    Calls on initialize_paths and initialize_sky1 to set up paths and read in sky data. 

    Args:
        filename (str): Name of the input data file.
        skytype (str): One of 'real_sky', 'rand_sky', or 'star_like'.
        star_rise_deg (float): Altitude (deg) above horizon to define star rise. The default is 10.
        sun_rise_deg (float): Sun altitude (deg) to define sunrise/set. The default is -12 (nautical twilight).
    
    Returns:
        skydict: Dictionary containing all relevant arrays and lists.
    """
    from decanopy.io.fileops import initialize_paths
    # set paths and read in magnitude data
    params = initialize_paths(skyname, skytype)
    filepath = params["filepath"]
    writepath = params["writepath"]
    mag_dict = params["mag_dict"]
    # create sky dictionary
    skydict = initialize_sky1(filepath)

    # add path and magnitude data to skydict
    skydict["mag_dict"] = mag_dict
    skydict["filepath"] = filepath
    skydict["writepath"] = writepath
    skydict["skyname"] = skyname
    skydict["skytype"] = skytype    

    return skydict
