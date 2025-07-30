## Some of these functions should ultimately go elsewhere!! 

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


def horizon_altitude_transform(horizon_window, alt_window):
    '''
    Given a horizon and altitude window, return an effective horizon window
    calculated via the midpoint of the altitude window. 
    '''
    a = (alt_window[0] + alt_window[1])/2 # average altitude of window (degrees)
    # a = alt_window[1] # maximum altitude check
    w = horizon_window[1] - horizon_window[0] # width of horizon window 
    ib = horizon_window[0] + w/2 # central position of the horizon window
    # calculate the transformation 
    sina = np.sin(np.radians(a)) # sin(a) with a in radians 
    cosa = np.cos(np.radians(a)) # cos(a) with a in radians 
    cosw = np.cos(np.radians(w)) # cos(w) with a in radians 
    #w1 = np.arccos(sina**2 + (cosa**2)*cosw) # transformed horizon width in radians
    w1 = np.arccos( (cosw - sina**2) /cosa**2 ) # transformed horizon width in radians
    w1 = np.degrees(w1) # transform back to degrees
    return(tuple(np.round((ib - w1/2, ib + w1/2),1)))

def horizonBins(horizon, bsize, gsize):
    # if gap size is zero, return as usual
    if gsize == 0:
        return((np.linspace(horizon[0], horizon[1], 8), 1))
    # otherwise, divide into smallest common denominator
    num = 7 * bsize + 6 * gsize + 1
    horizon_bins = np.linspace(horizon[0], horizon[1], num)
    #print(len(horizon_bins))
    # select indices to merge by bin/gap size
    ## first gaps
    ind_gap = np.zeros(14)
    ind_gap[1::2] = np.arange(0,7)
    ind_gap[::2] = np.arange(0,7)
    #print(ind_gap)
    ## then bins
    ind_bin = np.zeros(14)
    ind_bin[0:13] = ind_gap[1:]
    ind_bin[-1] = 7 
    #print(ind_bin)
    # add together with bin/gap sizes
    inds = ind_bin * bsize + ind_gap * gsize
    inds = inds.astype(int)
    #print(inds)
    # return only those bin and gap indices
    return((horizon_bins[inds], 2))


def isStarInWindow(alt_window, starAlt):
    '''
    Function to check whether star is in a given altitude window.
    '''
    (alt_min, alt_max) = alt_window
    if starAlt >= alt_min and starAlt <=alt_max:
        return(True)
    else:
        return(False)


def synRSC1star(date_ind, alt_window, horizon, bsize, gsize, sunSet, sunRise, starName, starAz, starAlt, starVis):
    '''
    Function which, when given a date, size of horizon, and single star's data*, 
    makes a synthetic Ramesside Star Clock.
    Also returns a dictionary of dbin_center for each star and time bin.
    '''
    # transform horizon limits according to the altitude window
    horizon = horizon_altitude_transform(horizon, alt_window)
    # define horizon bin limits (assuming equal binsize)
    (horizon_bins, skip) = horizonBins(horizon, bsize, gsize)
    horizon_bin_centers = (horizon_bins[0:-1] + horizon_bins[1:]) / 2 # centers of bins for calculating dbc
    # define time bin limits (assuming equal binsize)
    sset = sunSet[date_ind]
    sris = sunRise[date_ind + 1]
    time_ind_bins = np.round(np.linspace(sset, sris, 14, endpoint=False))
    rsc_table = []
    dbc_dict = {}  # distance from bin center dictionary, will look like {i: {starName: dbc}}
    for i in range(0, 13):
        ind = int(time_ind_bins[i])
        inWindow = isStarInWindow(alt_window, starAlt[ind])
        az_list = starAz[ind]
        if starVis[date_ind] and inWindow and az_list > horizon[0] and az_list < horizon[1]:
            row = np.histogram(az_list, horizon_bins)[0][::skip]
            dbc = np.round(min(np.abs(horizon_bin_centers - az_list)), 4)
            dbc_dict[i] = {starName: dbc}
        else:
            row = np.zeros(7, dtype=int)
        rsc_table.append(row)
    return np.array(rsc_table), dbc_dict

def synRSC(date_ind, alt_window, horizon, bsize, gsize, sunSet, sunRise, starlist, starAzlist, starAltlist, starVislist):
    '''
    Function to create synthetic Ramesside Star Clocks given data from several stars.
    Also returns a dictionary: {i: {starName: dbc, ...}, ...}
    '''
    df = pd.DataFrame(data=np.empty((13,7), dtype=str))
    dbc_dict = {i: {} for i in range(13)}  # initialize dbc_dict

    for i in range(len(starlist)):
        temp_table, dbc_temp = synRSC1star(
            date_ind, alt_window, horizon, bsize, gsize,
            sunSet, sunRise, starlist[i], starAzlist[i], starAltlist[i], starVislist[i]
        )
        inds = np.argwhere(temp_table == 1)
        for ind in inds:
            df.at[ind[0], ind[1]] += (starlist[i] + " ")
        # Merge dbc_temp into dbc_dict
        for k, v in dbc_temp.items():
            dbc_dict[k].update(v)
    df.columns = [-3, -2, -1, 0, 1, 2, 3]
    return df, dbc_dict



def mag_data(df, mag_dict):
    '''
    Given a data frame made with synRSC and a name-to-magnitude value dictionary, 
    this function will select the brightest star in each row (aka horizon bin) to 
    create a magnitude-selected Ramesside Star Clock.  
    '''
    # data frame to save magnitude-selected data
    df_mag = pd.DataFrame(data=np.empty((13,7), dtype=str))
    # iterate through df of all possible stars and select for magnitude
    for i in range(0, 13):
        sname = ""
        min_mag = 10 # all human visible magnitudes should be higher than this 
        # iterate through columns in row ( = horizon bins)
        for j in range(-3, 4):
            j *= -1 # testing something 
            dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
            for item in dlist:
                if mag_dict[item] < min_mag:
                    min_mag = mag_dict[item] # update brightest available star
                    cind = j + 3 # column index
                    sname = item # star name 
        if len(sname) > 1:                     
            df_mag.at[i, cind] = sname
    df_mag.columns = [-3, -2, -1, 0, 1, 2, 3]
    return(df_mag)


def name_or_mag_data(df, mag_dict, known_stars):
    '''
    Given a data frame made with synRSC, a name-to-magnitude value dictionary, and a known-star dictionary, 
    this function will select the known brightest star in each row (aka horizon bin) to 
    create a magnitude-selected Ramesside Star Clock;
    if no known stars, it defaults to brightest. 

    BEWARE: call this the alpha version of this function; it is the FARTHEST thing from elegant or optimized.   
    '''
    # data frame to save magnitude-selected data
    df_magname = pd.DataFrame(data=np.empty((13,7), dtype=str))
    # iterate through df of all possible stars and select for "known stars", then magnitude
    for i in range(0, 13): # for each row
        row_list = [] #list of stars in row
        for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
            dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
            row_list += dlist
        scand_list = [] # list of known candidate stars 
        for scand in row_list:
            if scand in known_stars:
                scand_list.append(scand)     
        # CASE 1: no previously known stars, choose by magnitude and add to known star list 
        if len(scand_list) == 0: 
            sname=''
            min_mag = 10 # all human visible magnitudes should be higher than this 
            for j in range(-3, 4):
                dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
                for item in dlist:
                    if mag_dict[item] < min_mag:
                        min_mag = mag_dict[item] # update brightest available star
                        cind = j + 3 # column index
                        sname = item # star name 
            if len(sname) > 1: # if it's found *no* stars, leave blank                
                df_magname.at[i, cind] = sname  
                known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary           
        # CASE 2: one known star, choose that one
        elif len(scand_list) == 1:
            sname = scand_list[0]
            for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
                dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
                if sname in dlist:
                    cind = j + 3 # column index
                    df_magname.at[i, cind] = sname 
                    #known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
        # CASE 3: several known stars, choose the brightest one 
        else:
            # first find the brightest available star
            min_mag = 10 # all human visible magnitudes should be higher than this 
            for scand in scand_list:
                if mag_dict[scand] < min_mag:
                    min_mag = mag_dict[scand] # update brightest available star
                    sname = scand # star name 
                    #known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
            for j in range(-3, 4): # now find position of star
                dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
                if sname in dlist:
                    cind = j + 3 # column index
                    df_magname.at[i, cind] = sname 
    # return                 
    df_magname.columns = [-3, -2, -1, 0, 1, 2, 3]
    return(df_magname, known_stars)    

# main and helper functions for full-choice algorithm

def cosorted_by_dict(list1, list2, sort_dict, ascending=True):
    """
    Co-sorts list1 and list2 based on the values of list1's elements in dict_.
    Returns the sorted lists.
    
    Args:
        list1: List of keys to sort by dict_ values.
        list2: List to be co-sorted with list1.
        dict_: Dictionary mapping elements of list1 to values.
        ascending: Sort order (default True for ascending).
    """
    zipped = sorted(zip(list1, list2), key=lambda x: sort_dict[x[0]], reverse=not ascending)
    if zipped:
        list1_sorted, list2_sorted = zip(*zipped)
        return list(list1_sorted), list(list2_sorted)
    else:
        return [], []


def sorted_magnitude_filter(row_list, ind_list, row_mag_list, threshold=0.5):
    ''' Function to filter out stars that are not within 0.5 mag of the first star.
        IMPORTANT: this function assumes that row_list and ind_list are already sorted by magnitude!
    '''    
    # set threshold for magnitude relative to the first star
    threshold += row_mag_list[0] 
    # filter out stars that are not within the threshold
    keep = [idx for idx, mag in enumerate(row_mag_list) if mag <= threshold]
    row_list = [row_list[k] for k in keep]
    ind_list = [ind_list[k] for k in keep]
    return (row_list, ind_list)

def full_choice_data_row(df, i, mag_dict, dbc_dict):
    '''
    Optimized implementation of the full-choice algorithm for one row.
    '''
    df_all = pd.DataFrame(data=np.empty((13, 10), dtype=str))  # 13 rows, 7 bins + 3 columns for separation, choice type, and comments
    #df_all.columns = [-3, -2, -1, 0, 1, 2, 3, "", "Code", " d_dbc"]

    row_list = []
    ind_list = []
    for j in range(-3, 4):
        dlist = list(filter(None, df[j][i].split(' ')))
        row_list += dlist
        ind_list += [j] * len(dlist)
    # STEP 0: trivial cases
    if len(row_list) == 1:      # Only one star in the row
        df_all.at[i, ind_list[0] + 3] = row_list[0]
        df_all.at[i, 8] = "S"
        if abs(ind_list[0] + 3) == 3:   # overwrite if also B3
            df_all.at[i, 8] = "B3"    
    elif len(row_list) == 0:    # No stars in the row
        df_all.at[i, 8] = "D"
    else:                       # Sort by magnitude
        #(row_listM, ind_listM) = cosorted_by_magnitude(row_list, ind_list, mag_dict)
        (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, mag_dict)
        #if row_listM != row_list or ind_listM != ind_list:
        #    print("Warning: row_listM and row_list are not equal after sorting by magnitude.")
        row_mag_list = [mag_dict[x] for x in row_list]  # Get magnitudes of sorted stars

        # STEP A: magnitude cut
        if row_mag_list[1] - row_mag_list[0] > 0.5: 
            df_all.at[i, ind_list[0] + 3] = row_list[0]
            df_all.at[i, 8] = "A"
            return df_all

        # Filter to only stars within 0.5 mag of the brightest
        row_list, ind_list = sorted_magnitude_filter(row_list, ind_list, row_mag_list)

        # Count occurrences by absolute bin index
        abs_counts = Counter(abs(x) for x in ind_list)

        # STEPS B0-B3: bin with exactly one star in order of priority (0, -1/+1, -2/+2,-3/+3)
        for bin_idx in range(4):
            if abs_counts[bin_idx] == 1:
                chosen_idx = [k for k, v in enumerate(ind_list) if abs(v) == bin_idx][0]
                df_all.at[i, ind_list[chosen_idx] + 3] = row_list[chosen_idx]
                df_all.at[i, 8] = f"B{bin_idx}"
                return df_all
            elif abs_counts[bin_idx] > 1:
                # If more than one star in this bin, go to STEP C
                break
        # STEP C: closest to bin center
        (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, dbc_dict)
        df_all.at[i, ind_list[0] + 3] = row_list[0]  # Assign the first star in the sorted list
        df_all.at[i, 8] = "C"  # Closest to bin center
        df_all.at[i, 9] = str(np.round(dbc_dict[row_list[1]] - dbc_dict[row_list[0]],2))
    return df_all

def full_choice_data(df, mag_dict, dbc_dict):
    '''
    Calls on full_choice_data_row to process all the rows of a given table. 
    '''
    df_out = pd.DataFrame(data=np.empty((13, 10), dtype=str))  # 13 rows, 7 bins + 3 columns for separation, choice type, and comments
    choices_dict = {'A': 0, 'B0': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'C': 0, 'D': 0, 'S': 0} # dictionary to count choice types
    # iterate through each row and apply full_choice_data_row
    for i in range(0, 13): # for each row
        df_row = full_choice_data_row(df, i, mag_dict, dbc_dict[i])
        df_out += df_row
        choices_dict[df_out.at[i, 8]] += 1
    df_out.columns = [-3, -2, -1, 0, 1, 2, 3, "", "Choice", " d_dbc"]  
    return df_out, choices_dict  


def dbc_data(df, dbc_dict):
    '''
    Implementation of the algorithm to minimize distance from center of bin. 
    '''
    df_dbc = pd.DataFrame(data=np.empty((13, 10), dtype=str))  # 13 rows, 7 bins + 3 columns for separation, choice type, and comments
    for i in range(0, 13):
        row_list = []
        ind_list = []
        dbc_dict_row = dbc_dict[i]
        # make list of candidates
        for j in range(-3, 4):
            # list of available stars
            dlist = list(filter(None, df[j][i].split(' ')))
            row_list += dlist
            ind_list += [j] * len(dlist)
        # make selections
        if len(row_list)==0:
            # mark if no star candidates in row
            df_dbc.at[i, 8] = "D"
        elif len(row_list)==1:   
            # only one option 
            df_dbc.at[i, ind_list[0] + 3] = row_list[0] 
            df_dbc.at[i, 8] = "S" # note single star
            df_dbc.at[i, 9] = str(np.round(dbc_dict_row[row_list[0]]))
        else: # if 2 or more stars
            # cosort by dbc
            (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, dbc_dict_row)
            df_dbc.at[i, ind_list[0] + 3] = row_list[0]  # Assign the first star in the sorted list
            df_dbc.at[i, 9] = str(np.round(dbc_dict_row[row_list[1]] - dbc_dict_row[row_list[0]],2))
    df_dbc.columns = [-3, -2, -1, 0, 1, 2, 3, "", "Code", "d_dbc"]
    return df_dbc