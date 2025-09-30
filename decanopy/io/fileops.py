# functions for i/o operations
# NOTE: Why the hell is _generate_time_grid here?

import sys
from pathlib import Path
import pandas as pd # for reading CSVs in initialize_paths 
from tqdm import tqdm # for progress bar in write_skyflow_file
from decanopy.skyflow.flow import calc_altaz, calc_sun_altaz # used in write_skyflow_file
import csv # for writing CSVs in write_skyflow_file # maybe use pandas instead?


def parse_args():
    parser = argparse.ArgumentParser('Sky Run')
    parser.add_argument('-yBC', '--yearBC', required=True)
    parser.add_argument('-d', '--decan', nargs='+', type=str, required=False, default=[])
    parser.add_argument('-m', '--month', required=False, default="01")
    parser.add_argument('-matchS', '--matchStellariumJD', required=False, default=True)
    parser.add_argument('-n', '--name', required=False, default="data")
    return parser.parse_args()

def _generate_time_grid(start, dhour, d4min):
    days = start + np.arange(0, 365)
    hours = dhour * np.arange(0, 24)
    minutes = d4min * np.arange(0, 15)
    return [
        (day, hour, mins)
        for day in days
        for hour in hours
        for mins in minutes
    ]

def write_skyflow_file(filename, hd_list, all_times, obj_list, checkpoint_file, start_idx, dS, obs_locale, dhour):
    with open(filename, "a", newline='') as file:
        writer = csv.writer(file, delimiter='|')
        if start_idx == 0:
            writer.writerow(hd_list)
        for idx, (day, hour, mins) in enumerate(tqdm(all_times[start_idx:], initial=start_idx, total=len(all_times))):
            temptime = day + hour + mins
            # Sun coordinates
            sun_az, sun_alt = calc_sun_altaz(temptime, obs_locale)
            local_time = Time(temptime - dS + (obs_locale.lon.deg/15.0) * dhour, format='jd').fits
            info = [
                f"{np.round(temptime, 10):.16f}",
                str(local_time),
                f"{sun_az:.3f}"[:-4],  # Remove trailing " deg"
                f"{sun_alt:.3f}"[:-4]
            ]
            # Decan coordinates
            altaz_pairs = [
                calc_altaz(Angle(obj.ra, unit="deg").hour, obj.dec, obs_locale, temptime)
                for obj in obj_list
            ]
            info.extend([f"{az:.3f}" for alt, az in altaz_pairs])
            info.extend([f"{alt:.3f}" for alt, az in altaz_pairs])
            writer.writerow(info)
            # Save checkpoint every N rows
            if idx % 100 == 0:
                with open(checkpoint_file, "w") as f:
                    f.write(str(idx + start_idx))

def initialize_paths(filename: str, skytype: str):
    """
    Initialize input/output paths for skyflow and RSC data based on skytype.
    Checks if the input file exists and aborts with an error if not.

    Args:
        filename (str): Name of the input data file.
        skytype (str): One of 'real_sky', 'rand_sky', or 'star_like'.

    Returns:
        dict: Dictionary with keys 'filepath', 'writepath', and 'mag_dict'.
    """
    from decanopy.config import (
        SKYFLOW_OUTPUT_REAL_SKY, RSC_OUTPUT_REAL_SKY, DEFAULT_STAR_NAMES,
        SKYFLOW_OUTPUT_RAND_SKY, RSC_OUTPUT_RAND_SKY, USER_INPUT_RAND_SKY,
        SKYFLOW_OUTPUT_STAR_LIKE, RSC_OUTPUT_STAR_LIKE, USER_INPUT_STAR_LIKE
    )
    import pandas as pd

    if skytype == "real_sky":
        filepath = SKYFLOW_OUTPUT_REAL_SKY / filename
        writepath = RSC_OUTPUT_REAL_SKY
        name_df = pd.read_csv(DEFAULT_STAR_NAMES, index_col=None, header=0, names=['Name', 'RA', 'Dec', 'Mag'])
    elif skytype == "rand_sky":
        filepath = SKYFLOW_OUTPUT_RAND_SKY / filename
        writepath = RSC_OUTPUT_RAND_SKY
        ic_filename = 'star_data_Mar-18-2024_1059.csv'
        name_df = pd.read_csv(USER_INPUT_RAND_SKY / ic_filename, index_col=None, header=0, names=['Name', 'RA', 'Dec', 'Mag'])
    elif skytype == "star_like":
        filepath = SKYFLOW_OUTPUT_STAR_LIKE / filename
        writepath = RSC_OUTPUT_STAR_LIKE
        ic_filename = 'star_data_other.csv'
        name_df = pd.read_csv(USER_INPUT_STAR_LIKE / ic_filename, index_col=None, header=0, names=['Name', 'RA', 'Dec', 'Mag'])
    else:
        raise ValueError(f"Unknown skytype: {skytype}. Please choose one of the following options: real_sky, rand_sky, or star_like.")

    # Check if file exists
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Input file not found: {filepath}\nPlease verify that the appropriate sky data exists in /data/skyflow/outpout/<skytype>/")

    mag_dict = name_df.set_index('Name')['Mag'].to_dict()

    return {
        "filepath": filepath,
        "writepath": writepath,
        "mag_dict": mag_dict
    }

def clobberCheck(filepath: Path, filename: str, clobberSave: bool):
    """
    Check if a file exists and handle according to clobberSave flag.

    Args:
        filepath (Path): Directory path where the file is located.
        filename (str): Name of the file to check.
        clobberSave (bool): If True, overwrite existing file. If False, abort if file exists.
    """
    full_path = filepath / filename
    if full_path.exists():
        if clobberSave:
            print(f"Warning: Overwriting existing file {full_path}")
        else:
            raise FileExistsError(f"File {full_path} already exists. To overwrite, set clobberSave=True.")     