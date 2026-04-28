"""
Run a decanOpy skyflow simulation for a real sky (precessed star coordinates) over a full year.

Computes and writes alt/az positions of the Sun and a set of stars (identified by
Hipparcos number) at 4-minute intervals over a full year, for a given year BCE and
observer location. Output is a pipe-delimited .txt file written to OUTPUT_SKYFLOW_DIR.

Usage
-----
    python run_sky_test.py -yBC 1300 [-n myrun]

Arguments
---------
    -yBC, --yearBC        Year BCE (required)
    -d,   --decan         One or more Hipparcos numbers (optional)
    -m,   --month         Starting month, zero-padded (default: '01')
    -matchS               Match Julian Date to Stellarium's convention (default: True)
    -n,   --name          Output file name prefix (default: 'data')

Notes
-----
    Precession uses the Vondrak algorithm via star_chart_spherical_projection.
    This may not exactly match Stellarium--small differences are to be expected!
"""

from pathlib import Path

import astropy.units as u
from astropy.coordinates import EarthLocation

from decanopy.config import OUTPUT_SKYFLOW_DIR, OBS_LAT, OBS_LON, OBS_HEIGHT, dhour, d4min
from decanopy.io.fileops import str2bool, parse_args, write_skyflow_file, load_checkpoint
from decanopy.skyflow.flow import precessed_coords, dS_offset, compute_start_jd, generate_time_grid

import warnings
warnings.simplefilter('ignore', UserWarning)

def main() -> None:
    args = parse_args()
    decans = list(args.decan)
    year = str(args.yearBC)
    month = str(args.month)
    matchStellariumJD = str2bool(args.matchStellariumJD)
    name = str(args.name)

    # skytype determination should go here (?)

    # Argument checks 
    if not (1 <= int(month) <= 12):
        raise ValueError(f"month must be between 01 and 12, got '{month}'")
    if int(year) <= 0:
        raise ValueError(f"yearBC must be a positive integer, got '{year}'")

    # Set place from config file
    obs_locale = EarthLocation(lat=OBS_LAT*u.deg, lon=OBS_LON*u.deg, height=OBS_HEIGHT*u.m)
    
    # Set Stellarium offset
    dS = dS_offset(year) if matchStellariumJD else 0 # Offset for JD matching with Stellarium

    # Set time
    start = compute_start_jd(year, month, obs_locale, dhour, dS)
    all_times = generate_time_grid(start, dhour, d4min)

    # Get star coordinates and headers
    obj_list, hd_list = precessed_coords(decans, year)

    # Output file
    OUTPUT_SKYFLOW_DIR.mkdir(parents=True, exist_ok=True) # make sure output directory exists
    filename = f"test_{name}{year}BC.txt"
    filename = OUTPUT_SKYFLOW_DIR / filename

    # Checkpointing
    checkpoint_file = filename.with_suffix(".checkpoint")
    start_idx = load_checkpoint(checkpoint_file)

    # Write to file 
    write_skyflow_file(filename, hd_list, all_times, obj_list, checkpoint_file, start_idx, dS, obs_locale, dhour)

if __name__ == "__main__":
    main()