# functions for making synthetic Ramesside Star Clocks
# NOTE: cosorted_by_dict() should maybe go into utils.py or something at some point 

import numpy as np
from pathlib import Path
import pandas as pd
from collections import Counter
from decanopy.io.fileops import clobberCheck

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

def mag_select_distinct(row_list, ind_list, mag_dict):
    '''
    Sort stars by magnitude and return only the lowest magnitude star in the row.  
    If there are multiple stars with the same lowest magnitude, pick one randomly.
    '''
    # sort and filter candidates by magnitude
    (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, mag_dict) # sort
    row_mag_list = [mag_dict[x] for x in row_list]  # get magnitudes of sorted stars
    row_list, ind_list = sorted_magnitude_filter(row_list, ind_list, row_mag_list, 0.0) # filter to only lowest mag value
    # how many stars?
    if len(row_list) == 1: # only one star in row
        return (row_list[0], ind_list[0])
    elif len(row_list) > 1: # more than one star in row
        # pick randomly
        pick = np.random.randint(0, len(row_list))
        return (row_list[pick], ind_list[pick])

def mag_select_distinct_dbc(row_list, ind_list, mag_dict, dbc_dict):
    '''
    Sort stars by magnitude and return only the lowest magnitude star in the row.  
    If there are multiple stars with the same lowest magnitude, pick the one closest to its bin centre (lowest dbc).
    '''
    # sort and filter candidates by magnitude
    (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, mag_dict) # sort
    row_mag_list = [mag_dict[x] for x in row_list]  # get magnitudes of sorted stars
    row_list, ind_list = sorted_magnitude_filter(row_list, ind_list, row_mag_list, 0.0) # filter to only lowest mag value
    # how many stars?
    if len(row_list) == 1: # only one star in row
        return (row_list[0], ind_list[0])
    elif len(row_list) > 1: # more than one star in row
        # pick lowest dbc
        min_idx = np.argmin([dbc_dict[star] for star in row_list])
        return (row_list[min_idx], ind_list[min_idx])
        # NOTE: this can theoretically still be degenerate if multiple stars have the same dbc value

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
        row_list = [] 
        ind_list = []
        # extract all stars in one row
        for j in range(-3, 4):
            dlist = list(filter(None, df[j][i].split(' ')))
            row_list += dlist
            ind_list += [j] * len(dlist)
        # select by magnitude
        (star, ind) = mag_select_distinct(row_list, ind_list, mag_dict)   
        df_mag.at[i, ind + 3] = star
    # name columns
    df_mag.columns = [-3, -2, -1, 0, 1, 2, 3]
    return(df_mag)

def name_or_mag_data(df, mag_dict, known_stars):
    '''
    Given a data frame made with synRSC, a name-to-magnitude value dictionary, and a known-star dictionary, 
    this function will select the known brightest star in each row (aka horizon bin) to 
    create a name-selected Ramesside Star Clock WITHOUT table duplicates ;
    if no known stars, it defaults to brightest.  
    '''
    # data frame to save magnitude-selected data
    df_magname = pd.DataFrame(data=np.empty((13,7), dtype=str))
    tablist = [] # track stars in this table to avoid duplicates
    # iterate through df of all possible stars and select for "known stars", then magnitude
    for i in range(0, 13): # for each row
        row_list = [] 
        ind_list = [] 
        for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
            dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
            row_list += dlist
            ind_list += [j] * len(dlist)
        ## filter out stars from THIS table and their indices
        indices_to_remove = [i for i, item in enumerate(row_list) if item in tablist]
        row_list = [item for i, item in enumerate(row_list) if i not in indices_to_remove]
        ind_list = [item for i, item in enumerate(ind_list) if i not in indices_to_remove]
        # check if there are any known stars in row
        scand_list = [] # list of known candidate stars 
        scand_list_ind = [] # list of known candidate star indices
        for scind in range(len(row_list)):
            if row_list[scind] in known_stars:
                scand_list.append(row_list[scind])
                scand_list_ind.append(ind_list[scind])
        if len(row_list) == 0:
            df_magname.at[i, 3] = "" # if no stars in row, leave blank
            #print("Blank row found!")
        elif len(scand_list) == 0:
            # select by magnitude
            (star, ind) = mag_select_distinct(row_list, ind_list, mag_dict)   
            df_magname.at[i, ind + 3] = star
            known_stars[star] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
        # CASE 2: one known star, choose that one
        elif len(scand_list) == 1:
            star = scand_list[0]
            df_magname.at[i,  scand_list_ind[0] + 3] = star
        # CASE 3: several known stars, choose the brightest one 
        else:
            # find the brightest available star and add to known star list
            (star, ind) = mag_select_distinct(scand_list, scand_list_ind, mag_dict)   
            df_magname.at[i, ind + 3] = star
            # known_stars[star] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
        # add chosen start to tablist 
        tablist.append(star)
    # return                 
    df_magname.columns = [-3, -2, -1, 0, 1, 2, 3]
    return(df_magname, known_stars)   

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
            df_dbc.at[i, 9] = str(np.round(dbc_dict_row[row_list[0]], 2))
        else: # if 2 or more stars
            # cosort by dbc
            (row_list, ind_list) = cosorted_by_dict(row_list, ind_list, dbc_dict_row)
            df_dbc.at[i, ind_list[0] + 3] = row_list[0]  # Assign the first star in the sorted list
            df_dbc.at[i, 9] = str(np.round(dbc_dict_row[row_list[1]] - dbc_dict_row[row_list[0]],2))
    df_dbc.columns = [-3, -2, -1, 0, 1, 2, 3, "", "Code", "d_dbc"]
    return df_dbc

def cbin_data(df, dbc_dict, mag_dict):
    '''
    Algorithm to prioritize centre bin, breaks degeneracy via dbc. 
    If no centre bin option, use mag select on other bins, break degeneracy via dbc.
    '''
    df_cb = pd.DataFrame(data=np.empty((13, 10), dtype=str))  # 13 rows, 7 bins + 3 columns for separation, choice type, and comments
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
            df_cb.at[i, 8] = "D"
        elif len(row_list)==1:   
            # only one option 
            df_cb.at[i, ind_list[0] + 3] = row_list[0] 
            df_cb.at[i, 8] = "S" # note single star
            df_cb.at[i, 9] = str(np.round(dbc_dict_row[row_list[0]], 2))
        else: # if 2 or more stars
            # check if centre bin has candidates
            if 0 in ind_list:
                # if so, select from centre bin candidates
                cb_list = [row_list[k] for k in range(len(row_list)) if ind_list[k]==0]
                (row_list, ind_list) = cosorted_by_dict(cb_list, [0]*len(cb_list), dbc_dict_row)
                df_cb.at[i, 3] = row_list[0]  # Assign the first star in the sorted list
                if len(row_list)>1:
                    df_cb.at[i, 9] = str(np.round(dbc_dict_row[row_list[1]] - dbc_dict_row[row_list[0]],2))
                    df_cb.at[i, 8] = "DBC0" # central chosen by dbc
                else:
                    df_cb.at[i, 9] = str(np.round(dbc_dict_row[row_list[0]],2))
                    df_cb.at[i, 8] = "CB0" # only one centre bin 
            # else, mag select from other bins
            else:
                (star, ind) = mag_select_distinct_dbc(row_list, ind_list, mag_dict, dbc_dict_row)   
                df_cb.at[i, ind + 3] = star
                df_cb.at[i, 8] = "mag" # chosen by mag (only dbc, not random))
    df_cb.columns = [-3, -2, -1, 0, 1, 2, 3, "", "Code", "d_dbc"]
    return df_cb

def init_synRSC_excel_writer(writepath, writename, horizon, alt_window, bsize, gsize):
    writer = pd.ExcelWriter(writepath / writename, engine='xlsxwriter')
    workbook = writer.book
    sheets = {}
    for name in ['RSCs', 'Mag Select', 'Name Select', 'DBC Select', 'CBin Select', 'Full Choice']:
        ws = workbook.add_worksheet(name)
        writer.sheets[name] = ws
        sheets[name] = ws
    cell_format = workbook.add_format()
    cell_format.set_font_size(11)
    # Metadata
    sheets['RSCs'].write(0, 0, f"horizon is {horizon}", cell_format)
    sheets['RSCs'].write(1, 0, f"alt window is {alt_window}", cell_format)
    sheets['RSCs'].write(2, 0, f"bsize = {bsize}", cell_format)
    sheets['RSCs'].write(3, 0, f"gsize = {gsize}", cell_format)
    return writer, workbook, sheets, cell_format

def write_rsc_table(i, df, sheets, writer, cell_format):
    df.to_excel(writer, sheet_name='RSCs', startrow=i * 15 + 5, startcol=0)
    sheets['RSCs'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)

def write_mag_select(i, df, mag_dict, sheets, writer, cell_format):
    df_mag = mag_data(df, mag_dict)
    df_mag.to_excel(writer, sheet_name='Mag Select', startrow=i * 15 + 5, startcol=0)
    sheets['Mag Select'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)

def write_name_select(i, df, mag_dict, known_stars_dict, sheets, writer, cell_format):
    df_name, known_stars_dict = name_or_mag_data(df, mag_dict, known_stars_dict)
    df_name.to_excel(writer, sheet_name='Name Select', startrow=i * 15 + 5, startcol=0)
    sheets['Name Select'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)
    sheets['Name Select'].write(4, 10, f"Number of known stars = {len(known_stars_dict)}", cell_format)
    df_dict = pd.DataFrame(list(known_stars_dict.items()), columns=["H-index", "'Known' index"])
    df_dict.to_excel(writer, sheet_name='Name Select', startrow=5, startcol=10, index=False)
    return known_stars_dict

def write_dbc_select(i, df, dbc_table, sheets, writer, cell_format):
    df_dbc = dbc_data(df, dbc_table)
    df_dbc.to_excel(writer, sheet_name='DBC Select', startrow=i * 15 + 5, startcol=0)
    sheets['DBC Select'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)

def write_cbin_select(i, df, mag_dict, dbc_table, sheets, writer, cell_format):
    df_cb = cbin_data(df, dbc_table, mag_dict)
    df_cb.to_excel(writer, sheet_name='CBin Select', startrow=i * 15 + 5, startcol=0)
    sheets['CBin Select'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)

def write_full_choice(i, df, mag_dict, dbc_table, sheets, writer, cell_format):
    df_choices, choices_dict = full_choice_data(df, mag_dict, dbc_table)
    df_choices.to_excel(writer, sheet_name='Full Choice', startrow=i * 15 + 5, startcol=0)
    sheets['Full Choice'].write(i * 15 + 5, 0, f"Table {i + 1}", cell_format)
    return choices_dict

def write_choices_summary(all_choices_dict, writer):
    df_dict = pd.DataFrame(list(all_choices_dict.items()), columns=["Code", "Count"])
    df_dict.to_excel(writer, sheet_name='Full Choice', startrow=5, startcol=12, index=False)

def write_synRSC_to_excel(writename, horizon, alt_window, bsize, gsize, skydict, clobberSave=False):

    # check if file exists and clobber if neede
    clobberCheck(skydict["writepath"], writename, clobberSave)

    # initalize excel writer
    writer, workbook, sheets, cell_format = init_synRSC_excel_writer(
        skydict["writepath"], writename, horizon, alt_window, bsize, gsize
    )
    # initialize helper dictionaries
    known_stars_dict = {}
    all_choices_dict = {}
    dbc_dict = {i: {} for i in range(24)}

    # loop over 24 tables
    for i in range(24):
        date = i * 15 # every 15 days
        # main synRSC table & update dbc dict
        df, dbc_table = synRSC(date, alt_window, horizon, bsize, gsize,
                            skydict["sunSet"], skydict["sunRise"], skydict["starlist"], skydict["starsAz"], skydict["starsAlt"], skydict["starVisList"])
        dbc_dict[i] = dbc_table # add to main dbc dictionary
        write_rsc_table(i, df, sheets, writer, cell_format)
        
        # write mag select 
        write_mag_select(i, df, skydict["mag_dict"], sheets, writer, cell_format)
        
        # write name select & update known stars dict
        known_stars_dict = write_name_select(i, df, skydict["mag_dict"], known_stars_dict, sheets, writer, cell_format)
        
        # write dbc select 
        write_dbc_select(i, df, dbc_table, sheets, writer, cell_format)

        # write cbin select 
        write_cbin_select(i, df, skydict["mag_dict"], dbc_table, sheets, writer, cell_format)
        
        # write full choice & update all choices dict
        choices_dict = write_full_choice(i, df, skydict["mag_dict"], dbc_table, sheets, writer, cell_format)
        for k, v in choices_dict.items():
            all_choices_dict[k] = all_choices_dict.get(k, 0) + v
    # write summary of all choices
    write_choices_summary(all_choices_dict, writer)
    writer.close()


### DEPRECATED VERSIONS BELOW

# def write_synRSC_to_excel(writepath, writename, horizon, alt_window, bsize, gsize, sunSet, sunRise, starlist, starsAz, starsAlt, starVisList, mag_dict):
    
#     # # Initialize Excel writer and sheets
#     # ## TODO: do I really need cell_format?
#     # writer, workbook, cell_format, rsc_wsheet, mag_wsheet, name_wsheet, dbc_wsheet, fc_wsheet = initialize_synRSC_excel(
#     #     writepath, writename, horizon, alt_window, bsize, gsize
#     # )

#     # Create Excel Writer Object from Pandas  
#     writer = pd.ExcelWriter(writepath / writename, engine='xlsxwriter')
#     workbook = writer.book

#     # Format
#     cell_format = workbook.add_format()
#     cell_format.set_font_size(11)

#     # Create worksheets with all possible stars 
#     rsc_wsheet = workbook.add_worksheet('RSCs')
#     writer.sheets['RSCs'] = rsc_wsheet
#     # Write metadata
#     rsc_wsheet.write(0, 0, "horizon is " + str(horizon), cell_format)
#     rsc_wsheet.write(1, 0, "alt window is " + str(alt_window), cell_format)
#     rsc_wsheet.write(2, 0, "bsize = " + str(bsize), cell_format)
#     rsc_wsheet.write(3, 0, "gsize = " + str(gsize), cell_format)

#     # add choices sheets
#     mag_wsheet = workbook.add_worksheet('Mag Select')
#     writer.sheets['Mag Select'] = mag_wsheet

#     name_wsheet = workbook.add_worksheet('Name Select')
#     writer.sheets['Name Select'] = name_wsheet

#     dbc_wsheet = workbook.add_worksheet('DBC Select')
#     writer.sheets['DBC Select'] = dbc_wsheet

#     fc_wsheet = workbook.add_worksheet('Full Choice')
#     writer.sheets['Full Choice'] = fc_wsheet

#     # create dictionaries to store data
#     known_stars_dict = {}
#     all_choices_dict = {}
#     dbc_dict = {i: {} for i in range(24)}  # initialize dbc dict for each table

#     # loop over 24 tables
#     for i in range(0, 24):
#         date = i * 15 # days from first day in decan data
#         # all star candidates
#         df, dbc_table = synRSC(date, alt_window, horizon, bsize, gsize, sunSet, sunRise, starlist, starsAz, starsAlt, starVisList)
#         dbc_dict[i] = dbc_table # store dbc_table for each table date
#         df.to_excel(writer, sheet_name='RSCs',startrow= i * 15 + 5, startcol=0)   
#         rsc_wsheet.write(i * 15 + 5,  0, "Table " + str(i + 1), format)

#         #add mag data
#         df_mag = mag_data(df, mag_dict)
#         df_mag.to_excel(writer, sheet_name='Mag Select',startrow= i * 15 + 5, startcol=0) 
#         mag_wsheet.write(i * 15 + 5,  0, "Table " + str(i + 1), format)

#         # add name or mag data
#         (df_name, known_stars_dict) = name_or_mag_data(df, mag_dict, known_stars_dict)
#         df_name.to_excel(writer, sheet_name='Name Select',startrow= i * 15 + 5, startcol=0) 
#         name_wsheet.write(i * 15 + 5,  0, "Table " + str(i + 1), format)
#         name_wsheet.write(4,  10, "Number of known stars = " + str(len((known_stars_dict))), format)
#         df_dict3 = pd.DataFrame(list(known_stars_dict.items()), columns=["H-index", "'Known' index"])
#         df_dict3.to_excel(writer, sheet_name='Name Select', startrow=5, startcol=10, index=False)

#         # add dbc data 
#         df_dbc = dbc_data(df, dbc_dict[i])
#         df_dbc.to_excel(writer, sheet_name='DBC Select', startrow= i * 15 + 5, startcol=0) 
#         dbc_wsheet.write(i * 15 + 5,  0, "Table " + str(i + 1), format)

#         # add choices data 
#         df_choices, choices_dict  = full_choice_data(df, mag_dict, dbc_table)
#         df_choices.to_excel(writer, sheet_name='Full Choice', startrow= i * 15 + 5, startcol=0) 
#         fc_wsheet.write(i * 15 + 5,  0, "Table " + str(i + 1), format)
        
#         # dynamically update all_choices_dict
#         for k, v in choices_dict.items():
#             all_choices_dict[k] = all_choices_dict.get(k, 0) + v
#     # write final all_choices_dict and close    
#     df_dict5 = pd.DataFrame(list(all_choices_dict.items()), columns=["Code", "Count"])
#     df_dict5.to_excel(writer, sheet_name='Full Choice', startrow=5, startcol=12, index=False)    
#     writer.close()

# def name_or_mag_data(df, mag_dict, known_stars):
#     '''
#     Given a data frame made with synRSC, a name-to-magnitude value dictionary, and a known-star dictionary, 
#     this function will select the known brightest star in each row (aka horizon bin) to 
#     create a name-selected Ramesside Star Clock ALLOWING table duplicates;
#     if no known stars, it defaults to brightest. 
#     '''
#     # data frame to save magnitude-selected data
#     df_magname = pd.DataFrame(data=np.empty((13,7), dtype=str))
#     # iterate through df of all possible stars and select for "known stars", then magnitude
#     for i in range(0, 13): # for each row
#         row_list = [] 
#         ind_list = [] 
#         for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
#             dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#             row_list += dlist
#             ind_list += [j] * len(dlist)
#         # check if there are any known stars in row
#         scand_list = [] # list of known candidate stars 
#         scand_list_ind = [] # list of known candidate star indices
#         for scind in range(len(row_list)):
#             if row_list[scind] in known_stars:
#                 scand_list.append(row_list[scind])
#                 scand_list_ind.append(ind_list[scind])
#         # CASE 1: no previously known stars, choose by magnitude and add to known star list 
#         if len(scand_list) == 0:
#             # select by magnitude
#             (star, ind) = mag_select_distinct(row_list, ind_list, mag_dict, df_magname)   
#             df_magname.at[i, ind + 3] = star
#             known_stars[star] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
#         # CASE 2: one known star, choose that one
#         elif len(scand_list) == 1:
#             df_magname.at[i,  scand_list_ind[0] + 3] = scand_list[0] 
#         # CASE 3: several known stars, choose the brightest one 
#         else:
#             # find the brightest available star and add to known star list
#             (star, ind) = mag_select_distinct(scand_list, scand_list_ind, mag_dict, df_magname)   
#             df_magname.at[i, ind + 3] = star
#             # known_stars[star] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
#     # return                 
#     df_magname.columns = [-3, -2, -1, 0, 1, 2, 3]
#     return(df_magname, known_stars)    

# def mag_data(df, mag_dict):
#     '''
#     Given a data frame made with synRSC and a name-to-magnitude value dictionary, 
#     this function will select the brightest star in each row (aka horizon bin) to 
#     create a magnitude-selected Ramesside Star Clock.  
#     '''
#     # data frame to save magnitude-selected data
#     df_mag = pd.DataFrame(data=np.empty((13,7), dtype=str))
#     # iterate through df of all possible stars and select for magnitude
#     for i in range(0, 13):
#         sname = ""
#         min_mag = 10 # all human visible magnitudes should be higher than this 
#         # iterate through columns in row ( = horizon bins)
#         for j in range(-3, 4):
#             #j *= -1 # testing something 
#             dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#             for item in dlist:
#                 if mag_dict[item] < min_mag:
#                     min_mag = mag_dict[item] # update brightest available star
#                     cind = j + 3 # column index
#                     sname = item # star name 
#         if len(sname) > 1:                     
#             df_mag.at[i, cind] = sname
#     df_mag.columns = [-3, -2, -1, 0, 1, 2, 3]
#     return(df_mag)


# def name_or_mag_data(df, mag_dict, known_stars):
#     '''
#     Given a data frame made with synRSC, a name-to-magnitude value dictionary, and a known-star dictionary, 
#     this function will select the known brightest star in each row (aka horizon bin) to 
#     create a magnitude-selected Ramesside Star Clock;
#     if no known stars, it defaults to brightest. 

#     BEWARE: call this the alpha version of this function; it is the FARTHEST thing from elegant or optimized.   
#     '''
#     # data frame to save magnitude-selected data
#     df_magname = pd.DataFrame(data=np.empty((13,7), dtype=str))
#     # iterate through df of all possible stars and select for "known stars", then magnitude
#     for i in range(0, 13): # for each row
#         row_list = [] #list of stars in row
#         for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
#             dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#             row_list += dlist
#         scand_list = [] # list of known candidate stars 
#         for scand in row_list:
#             if scand in known_stars:
#                 scand_list.append(scand)     
#         # CASE 1: no previously known stars, choose by magnitude and add to known star list 
#         if len(scand_list) == 0: 
#             sname=''
#             min_mag = 10 # all human visible magnitudes should be higher than this 
#             for j in range(-3, 4):
#                 dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#                 for item in dlist:
#                     if mag_dict[item] < min_mag:
#                         min_mag = mag_dict[item] # update brightest available star
#                         cind = j + 3 # column index
#                         sname = item # star name 
#             if len(sname) > 1: # if it's found *no* stars, leave blank                
#                 df_magname.at[i, cind] = sname  
#                 known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary           
#         # CASE 2: one known star, choose that one
#         elif len(scand_list) == 1:
#             sname = scand_list[0]
#             for j in range(-3, 4): # iterate through columns in row ( = horizon bins)
#                 dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#                 if sname in dlist:
#                     cind = j + 3 # column index
#                     df_magname.at[i, cind] = sname 
#                     #known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
#         # CASE 3: several known stars, choose the brightest one 
#         else:
#             # first find the brightest available star
#             min_mag = 10 # all human visible magnitudes should be higher than this 
#             for scand in scand_list:
#                 if mag_dict[scand] < min_mag:
#                     min_mag = mag_dict[scand] # update brightest available star
#                     sname = scand # star name 
#                     #known_stars[sname] = "K" + str(len(known_stars)).zfill(2) # update known star dictionary 
#             for j in range(-3, 4): # now find position of star
#                 dlist = list(filter(None, df[j][i].split(' '))) # split into star names and filter out empty strings 
#                 if sname in dlist:
#                     cind = j + 3 # column index
#                     df_magname.at[i, cind] = sname 
#     # return                 
#     df_magname.columns = [-3, -2, -1, 0, 1, 2, 3]
#     return(df_magname, known_stars)    

# main and helper functions for full-choice algorithm

# def initialize_synRSC_excel(writepath, writename, horizon, alt_window, bsize, gsize):
#     """
#     Helper to initialize Excel writer, sheets, and formatting for synRSC output.
#     Returns: writer, workbook, worksheet, worksheet2, worksheet3, worksheet4, worksheet5, format
#     """
#     # Create Excel Writer Object from Pandas  
#     writer = pd.ExcelWriter(writepath / writename, engine='xlsxwriter')
#     workbook = writer.book

#     # Format
#     cell_format = workbook.add_format()
#     cell_format.set_font_size(11)

#     # Create worksheets with all possible stars 
#     rsc_wsheet = workbook.add_worksheet('RSCs')
#     writer.sheets['RSCs'] = rsc_wsheet
#     # Write metadata
#     rsc_wsheet.write(0, 0, "horizon is " + str(horizon), cell_format)
#     rsc_wsheet.write(1, 0, "alt window is " + str(alt_window), cell_format)
#     rsc_wsheet.write(2, 0, "bsize = " + str(bsize), cell_format)
#     rsc_wsheet.write(3, 0, "gsize = " + str(gsize), cell_format)

#     # add choices sheets
#     mag_wsheet = workbook.add_worksheet('Mag Select')
#     writer.sheets['Mag Select'] = mag_wsheet

#     name_wsheet = workbook.add_worksheet('Name Select')
#     writer.sheets['Name Select'] = name_wsheet

#     dbc_wsheet = workbook.add_worksheet('DBC Select')
#     writer.sheets['DBC Select'] = dbc_wsheet

#     cbin_wsheet = workbook.add_worksheet('Cbin Select')
#     writer.sheets['Cbin Select'] = cbin_wsheet

#     fc_wsheet = workbook.add_worksheet('Full Choice')
#     writer.sheets['Full Choice'] = fc_wsheet

#     # add compare sheet here

#     # return all objects
#     return writer, workbook, cell_format, rsc_wsheet, mag_wsheet, name_wsheet, dbc_wsheet, cbin_wsheet, fc_wsheet