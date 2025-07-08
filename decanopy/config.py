from pathlib import Path

# Project root (you may eventually want this to be dynamic)
ROOT_DIR = Path(__file__).resolve().parents[1]

# Static packaged data
PACKAGE_DATA = ROOT_DIR / "decanopy" / "data"
DEFAULT_STAR_DATA = PACKAGE_DATA / "star_data.csv"
DEFAULT_STAR_NAMES = PACKAGE_DATA / "star_data_names.csv"

# User input/output
INPUT_DIR = ROOT_DIR / "data" / "input" / "skyflow"
OUTPUT_SKYFLOW = ROOT_DIR / "data" / "output" / "skyflow"
OUTPUT_RSC = ROOT_DIR / "data" / "output" / "rsc"

# Output by sky type
REAL_SKY = OUTPUT_SKYFLOW / "real_sky"
RAND_SKY = OUTPUT_SKYFLOW / "rand_sky"
STAR_LIKE = OUTPUT_SKYFLOW / "star_like"