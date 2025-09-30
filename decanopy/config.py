from pathlib import Path

# === Project Root ===
ROOT_DIR = Path(__file__).resolve().parents[1]

# === Packaged Static Data (within decanopy/) ===
PACKAGE_DATA = ROOT_DIR / "decanopy" / "data"
DEFAULT_STAR_DATA = PACKAGE_DATA / "star_data.csv"
DEFAULT_STAR_NAMES = PACKAGE_DATA / "star_data_names.csv"

# === User Input Data ===
USER_INPUT_SKYFLOW_DIR = ROOT_DIR / "data" / "input" / "skyflow"
USER_INPUT_REAL_SKY = USER_INPUT_SKYFLOW_DIR / "real_sky"
USER_INPUT_RAND_SKY = USER_INPUT_SKYFLOW_DIR / "rand_sky"
USER_INPUT_STAR_LIKE = USER_INPUT_SKYFLOW_DIR / "star_like"

# === Skyflow Outputs (generated star motions) ===
OUTPUT_SKYFLOW_DIR = ROOT_DIR / "data" / "output" / "skyflow"
SKYFLOW_OUTPUT_REAL_SKY = OUTPUT_SKYFLOW_DIR / "real_sky"
SKYFLOW_OUTPUT_RAND_SKY = OUTPUT_SKYFLOW_DIR / "rand_sky"
SKYFLOW_OUTPUT_STAR_LIKE = OUTPUT_SKYFLOW_DIR / "star_like"

# === RSC Outputs (analysis results based on skyflow outputs) ===
OUTPUT_RSC_DIR = ROOT_DIR / "data" / "output" / "rsc"
RSC_OUTPUT_REAL_SKY = OUTPUT_RSC_DIR / "real_sky"
RSC_OUTPUT_RAND_SKY = OUTPUT_RSC_DIR / "rand_sky"
RSC_OUTPUT_STAR_LIKE = OUTPUT_RSC_DIR / "star_like"

# === Default observing locale latitude===
# Luxor, Egypt
OBS_LAT = 25.6989   # degrees
OBS_LON = 32.6421   # degrees
OBS_HEIGHT = 89     # meters