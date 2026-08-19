"""
schedule_analysis.py
====================
Utilities for reading scheduling Excel files and building heatmap-ready data structures.

Expected filename convention:
    results_h<hlo>-<hhi>_a<alo>-<ahi>.xlsx
    e.g.  results_h165-195_a0-30.xlsx

Each file has 5 data sheets plus "RSCs" (used only for fail-condition checking):
    Mag Select, Name Select, DBC Select, CBin Select, Full Choice

scan_folder() returns a dict of long-form DataFrames keyed by short sheet name:
    {
        "Mag":   DataFrame,
        "Name":  DataFrame,
        "DBC":   DataFrame,
        "CBin":  DataFrame,
        "Full":  DataFrame,
    }

Each DataFrame has columns:
    horizon_lo, horizon_hi, alt_lo, alt_hi,
    horizon_width,            ← hi - lo  (used for tick labels)
    alt_hi_label,             ← just the max altitude (used for tick labels)
    fail1, fail2, total_empty_rows,
    value                     ← unique record count minus reference, or FAIL1/FAIL2 code
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 1.  Filename helpers
# ---------------------------------------------------------------------------

_FNAME_RE = re.compile(
    r"results_h(?P<hlo>[\d.]+)-(?P<hhi>[\d.]+)_a(?P<alo>[\d.]+)-(?P<ahi>[\d.]+)\.xlsx",
    re.IGNORECASE,
)


def parse_filename(path: str | Path) -> dict | None:
    """Return {'horizon': (hlo,hhi), 'alt': (alo,ahi)} or None if no match."""
    m = _FNAME_RE.search(Path(path).name)
    if not m:
        return None
    return {
        "horizon": (int(float(m["hlo"])), int(float(m["hhi"]))),
        "alt":     (int(float(m["alo"])), int(float(m["ahi"]))),
    }


def build_filename(horizon: tuple, alt: tuple, prefix: str = "results") -> str:
    return f"{prefix}_h{horizon[0]}-{horizon[1]}_a{alt[0]}-{alt[1]}.xlsx"


# ---------------------------------------------------------------------------
# 2.  Single-sheet parser
# ---------------------------------------------------------------------------

_RCELL_RE = re.compile(r"R\d{4}")    # default: R#### format
_HCELL_RE = re.compile(r"H\d{6}")    # alternative: H###### format

#_RCELL_RE = re.compile(r"R\d{4}")

# def _parse_cell(val) -> list:
#     if pd.isna(val):
#         return []
#     return _RCELL_RE.findall(str(val))

def _parse_cell(val, pattern=None) -> list:
    if pd.isna(val):
        return []
    return (pattern or _RCELL_RE).findall(str(val))


def read_sheet(filepath: str | Path, sheet_name: str, has_metadata: bool = False,
               cell_re=None) -> dict:
    """
    Parse one sheet of a scheduling Excel file.

    Parameters
    ----------
    has_metadata : if True, parse horizon/alt/bsize/gsize from the first 4 rows
                   (only "RSCs" has these). For all other sheets, pass False
                   (default) and the parser skips straight to the first Table header.
    cell_re      : compiled regex for parsing cell values (default: R#### pattern).
                   Pass sa._HCELL_RE for H###### format, or your own re.compile().

    Returns dict with keys:
        'horizon', 'alt', 'bsize', 'gsize', 'tables'
    where 'tables' is {table_number: DataFrame} with list-of-str cells.
    """
    pattern = cell_re or _RCELL_RE
    raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None, dtype=str,
                        usecols=range(8))   # columns A-H only; ignore extra columns

    def _nums(s):
        return tuple(int(float(x)) for x in re.findall(r"[\d.]+", str(s)))

    if has_metadata:
        horizon = _nums(raw.iloc[0, 0])[:2]
        alt     = _nums(raw.iloc[1, 0])[:2]
        bsize   = int(float(re.search(r"[\d.]+", str(raw.iloc[2, 0])).group()))
        gsize   = int(float(re.search(r"[\d.]+", str(raw.iloc[3, 0])).group()))
    else:
        horizon, alt, bsize, gsize = None, None, None, None

    COL_LABELS = [-3, -2, -1, 0, 1, 2, 3]
    tables    = {}
    table_num = None
    rows_buf  = {}

    for ri in range(5, len(raw)):
        row   = raw.iloc[ri]
        cell0 = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""

        if cell0.startswith("Table"):
            if table_num is not None and rows_buf:
                tables[table_num] = _buf_to_df(rows_buf, COL_LABELS)
            table_num = int(re.search(r"\d+", cell0).group())
            rows_buf  = {}
            continue

        if cell0.strip() == "" and all(pd.isna(row.iloc[1:8])):
            continue

        if table_num is not None:
            try:
                row_idx = int(float(cell0))
            except (ValueError, TypeError):
                continue
            rows_buf[row_idx] = [_parse_cell(row.iloc[c], pattern) for c in range(1, 8)]

    if table_num is not None and rows_buf:
        tables[table_num] = _buf_to_df(rows_buf, COL_LABELS)

    return {"horizon": horizon, "alt": alt, "bsize": bsize, "gsize": gsize, "tables": tables}


def _buf_to_df(rows_buf: dict, col_labels: list) -> pd.DataFrame:
    df = pd.DataFrame.from_dict(
        {ri: vals for ri, vals in sorted(rows_buf.items())},
        orient="index",
        columns=col_labels,
    )
    df.index.name = "row"
    return df



# def read_sheet(filepath: str | Path, sheet_name: str, has_metadata: bool = False) -> dict:
#     """
#     Parse one sheet of a scheduling Excel file.

#     Parameters
#     ----------
#     has_metadata : if True, parse horizon/alt/bsize/gsize from the first 4 rows
#                    (only "RSCs" has these). For all other sheets, pass False
#                    (default) and the parser skips straight to the first Table header.

#     Returns dict with keys:
#         'horizon', 'alt', 'bsize', 'gsize', 'tables'
#     where 'tables' is {table_number: DataFrame} with list-of-str cells.
#     """
#     raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None, dtype=str,
#                         usecols=range(8))   # columns A-H only; ignore extra columns

#     def _nums(s):
#         return tuple(int(float(x)) for x in re.findall(r"[\d.]+", str(s)))

#     if has_metadata:
#         horizon = _nums(raw.iloc[0, 0])[:2]
#         alt     = _nums(raw.iloc[1, 0])[:2]
#         bsize   = int(float(re.search(r"[\d.]+", str(raw.iloc[2, 0])).group()))
#         gsize   = int(float(re.search(r"[\d.]+", str(raw.iloc[3, 0])).group()))
#     else:
#         horizon, alt, bsize, gsize = None, None, None, None

#     COL_LABELS = [-3, -2, -1, 0, 1, 2, 3]
#     tables    = {}
#     table_num = None
#     rows_buf  = {}

#     for ri in range(5, len(raw)):
#         row   = raw.iloc[ri]
#         cell0 = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""

#         if cell0.startswith("Table"):
#             if table_num is not None and rows_buf:
#                 tables[table_num] = _buf_to_df(rows_buf, COL_LABELS)
#             table_num = int(re.search(r"\d+", cell0).group())
#             rows_buf  = {}
#             continue

#         if cell0.strip() == "" and all(pd.isna(row.iloc[1:8])):
#             continue

#         if table_num is not None:
#             try:
#                 row_idx = int(float(cell0))
#             except (ValueError, TypeError):
#                 continue
#             rows_buf[row_idx] = [_parse_cell(row.iloc[c]) for c in range(1, 8)]

#     if table_num is not None and rows_buf:
#         tables[table_num] = _buf_to_df(rows_buf, COL_LABELS)

#     return {"horizon": horizon, "alt": alt, "bsize": bsize, "gsize": gsize, "tables": tables}


# def _buf_to_df(rows_buf: dict, col_labels: list) -> pd.DataFrame:
#     df = pd.DataFrame.from_dict(
#         {ri: vals for ri, vals in sorted(rows_buf.items())},
#         orient="index",
#         columns=col_labels,
#     )
#     df.index.name = "row"
#     return df


# ---------------------------------------------------------------------------
# 3.  Fail-condition checks  (always run on "RSCs" sheet)
# ---------------------------------------------------------------------------

def _is_row_empty(row_series) -> bool:
    return all(len(cell) == 0 for cell in row_series)


def check_fail_conditions(tables: dict) -> dict:
    """
    Fail 1 : total empty rows across all tables >= 10
    Fail 2 : any single table has >= 2 empty rows
    """
    empty_per_table = {
        tnum: sum(_is_row_empty(df.loc[ri]) for ri in df.index)
        for tnum, df in tables.items()
    }
    total_empty = sum(empty_per_table.values())
    fail1 = total_empty >= 10
    fail2 = any(v >= 2 for v in empty_per_table.values())
    return {
        "fail1":                fail1,
        "fail2":                fail2,
        "total_empty_rows":     total_empty,
        "empty_rows_per_table": empty_per_table,
    }


# ---------------------------------------------------------------------------
# 4.  Metrics
# ---------------------------------------------------------------------------

COL_OFFSETS = [-3, -2, -1, 0, 1, 2, 3]

def count_unique_records(tables: dict) -> int:
    """Count unique R#### IDs appearing anywhere in all tables."""
    seen = set()
    for df in tables.values():
        for col in df.columns:
            for cell in df[col]:
                seen.update(cell)
    return len(seen)


def compute_histogram(tables: dict) -> dict:
    """
    Count total number of R#### entries (stars, duplicates included) in each
    column offset k in {-3,-2,-1,0,1,2,3} across all 24 tables.

    Returns {k: count} dict.
    """
    hist = {k: 0 for k in COL_OFFSETS}
    for df in tables.values():
        for k in COL_OFFSETS:
            if k in df.columns:
                hist[k] += sum(len(cell) for cell in df[k])
    return hist


def compute_total_stars(tables: dict) -> int:
    """Total number of populated rows across all tables (24 * 13 - empty rows)."""
    empty = sum(
        sum(
            all(len(cell) == 0 for cell in df.loc[ri])
            for ri in df.index
        )
        for df in tables.values()
    )
    return 24 * 13 - empty


def compute_ratios(hist: dict):
    """
    Compute R0 and R1 from a histogram dict {k: count}.

    R0 = C(0) / sum(C(k) for k != 0)
    R1 = (C(-1) + C(+1)) / (C(-3) + C(-2) + C(+2) + C(+3))

    Returns (R0, R1) — either may be NaN if denominator is zero.
    """
    denom_r0 = sum(hist[k] for k in COL_OFFSETS if k != 0)
    R0 = hist[0] / denom_r0 if denom_r0 > 0 else float("nan")

    denom_r1 = hist[-3] + hist[-2] + hist[2] + hist[3]
    R1 = (hist[-1] + hist[1]) / denom_r1 if denom_r1 > 0 else float("nan")

    denom_s = hist[-1] + hist[-2]+ hist[-3]
    S  = (hist[1] + hist[2]+ hist[3])/denom_s  if denom_s > 0 else float("nan")
    return R0, R1, S


def compute_all_metrics(tables: dict) -> dict:
    """
    Compute all metrics for a passing sheet. Returns a dict with keys:
        unicount, hist, total_stars, R0, R1
    """
    hist       = compute_histogram(tables)
    R0, R1, S     = compute_ratios(hist)
    return {
        "unicount":    count_unique_records(tables),
        "hist":        hist,
        "total_stars": compute_total_stars(tables),
        "R0":          R0,
        "R1":          R1,
        "S":           S,
    }


def plot_histogram(
    df: pd.DataFrame,
    horizon_width: int,
    alt_hi: int,
    title: str = None,
    figsize: tuple = (7, 4),
    ax=None,
):
    """
    Plot the bin histogram for a specific (horizon_width, alt_hi) cell.

    Parameters
    ----------
    df            : one sheet's DataFrame from scan_folder() output
    horizon_width : horizon width value (hi - lo) to select
    alt_hi        : maximum altitude value to select
    title         : plot title (auto-generated if None)
    figsize       : figure size
    ax            : existing Axes to draw on

    Returns (fig, ax).
    """
    import matplotlib.pyplot as plt

    row = df[(df["horizon_width"] == horizon_width) & (df["alt_hi_label"] == alt_hi)]
    if row.empty:
        raise ValueError(f"No data found for horizon_width={horizon_width}, alt_hi={alt_hi}")

    hist = row.iloc[0]["hist"]
    if hist is None:
        raise ValueError(f"Cell (horizon_width={horizon_width}, alt_hi={alt_hi}) is a failed cell — no histogram.")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    keys = COL_OFFSETS
    vals = [hist[k] for k in keys]
    colors = ["#4477AA" if k != 0 else "#EE6677" for k in keys]
    ax.bar(keys, vals, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(keys)
    ax.set_xticklabels([str(k) for k in keys])
    ax.set_xlabel("Bin offset k")
    ax.set_ylabel("Star count C(k)")
    ax.set_title(title or f"Histogram — horizon width {horizon_width}°, max alt {alt_hi}")
    fig.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# 5.  Scan a folder → dict of long-form DataFrames, one per data sheet
# ---------------------------------------------------------------------------

FAIL1_CODE = -1001   # >= 10 total empty rows across all tables
FAIL2_CODE = -1002   # any single table has >= 2 empty rows
FAILB_CODE = -1003   # both conditions true simultaneously

# Maps short key → full Excel sheet name
DATA_SHEETS = {
    "Mag":  "Mag Select",
    "Name": "Name Select",
    "DBC":  "DBC Select",
    "CBin": "CBin Select",
    "Full": "Full Choice",
}


# def scan_folder(
#     folder: str | Path,
#     glob_pattern: str = "results_h*_a*.xlsx",
# ) -> dict:
#     """
#     Walk *folder* for all matching Excel files and return a dict of DataFrames,
#     one per data sheet (keys: "Mag", "Name", "DBC", "CBin", "Full").

#     Fail conditions are evaluated from the "RSCs" sheet and propagate to all
#     metric columns. Metrics are only computed for passing files.

#     Returns
#     -------
#     dict of DataFrames keyed by short sheet name, each with columns:
#         horizon_lo, horizon_hi, alt_lo, alt_hi,
#         horizon_width, alt_hi_label,
#         fail1, fail2, total_empty_rows,
#         unicount, hist, total_stars, R0, R1
#     Failing cells store the appropriate FAIL*_CODE in every metric column.
#     """
#     folder  = Path(folder)
#     buffers = {key: [] for key in DATA_SHEETS}

#     for fpath in sorted(folder.glob(glob_pattern)):
#         params = parse_filename(fpath)
#         if params is None:
#             continue

#         # --- fail conditions from RSCs ----------------------------------------
#         try:
#             rsc_data = read_sheet(fpath, "RSCs", has_metadata=True)
#         except Exception as e:
#             print(f"[WARN] Could not read RSCs in {fpath.name}: {e}")
#             continue

#         fc = check_fail_conditions(rsc_data["tables"])

#         if fc["fail1"] and fc["fail2"]:
#             fail_code = FAILB_CODE
#         elif fc["fail1"]:
#             fail_code = FAIL1_CODE
#         elif fc["fail2"]:
#             fail_code = FAIL2_CODE
#         else:
#             fail_code = None

#         h = params["horizon"]
#         a = params["alt"]
#         base = {
#             "horizon_lo":       h[0],
#             "horizon_hi":       h[1],
#             "alt_lo":           a[0],
#             "alt_hi":           a[1],
#             "horizon_width":    h[1] - h[0],
#             "alt_hi_label":     a[1],
#             "fail1":            fc["fail1"],
#             "fail2":            fc["fail2"],
#             "total_empty_rows": fc["total_empty_rows"],
#         }

#         # --- compute metrics for each data sheet ------------------------------
#         for key, sheet_name in DATA_SHEETS.items():
#             if fail_code is not None:
#                 metrics = {
#                     "unicount":    fail_code,
#                     "hist":        None,
#                     "total_stars": fail_code,
#                     "R0":          fail_code,
#                     "R1":          fail_code,
#                     "S":           fail_code,
#                 }
#             else:
#                 try:
#                     data    = read_sheet(fpath, sheet_name, has_metadata=False)
#                     metrics = compute_all_metrics(data["tables"])
#                 except Exception as e:
#                     print(f"[WARN] Could not read '{sheet_name}' in {fpath.name}: {e}")
#                     metrics = {
#                         "unicount":    np.nan,
#                         "hist":        None,
#                         "total_stars": np.nan,
#                         "R0":          np.nan,
#                         "R1":          np.nan,
#                         "S":           np.nan,
#                     }

#             buffers[key].append({**base, **metrics})

#     return {key: pd.DataFrame(rows) for key, rows in buffers.items()}


def scan_folder(
    folder: str | Path,
    glob_pattern: str = "results_h*_a*.xlsx",
    cell_re=None,
) -> dict:
    """
    Walk *folder* for all matching Excel files and return a dict of DataFrames,
    one per data sheet (keys: "Mag", "Name", "DBC", "CBin", "Full").

    Fail conditions are evaluated from the "RSCs" sheet and propagate to all
    metric columns. Metrics are only computed for passing files.

    Parameters
    ----------
    cell_re : compiled regex for parsing cell values (default: R#### pattern).
              Pass sa._HCELL_RE for H###### format, or your own re.compile().

    Returns
    -------
    dict of DataFrames keyed by short sheet name, each with columns:
        horizon_lo, horizon_hi, alt_lo, alt_hi,
        horizon_width, alt_hi_label,
        fail1, fail2, total_empty_rows,
        unicount, hist, total_stars, R0, R1, S
    Failing cells store the appropriate FAIL*_CODE in every metric column.
    """
    folder  = Path(folder)
    buffers = {key: [] for key in DATA_SHEETS}

    for fpath in sorted(folder.glob(glob_pattern)):
        params = parse_filename(fpath)
        if params is None:
            continue

        # --- fail conditions from RSCs ----------------------------------------
        try:
            rsc_data = read_sheet(fpath, "RSCs", has_metadata=True, cell_re=cell_re)
        except Exception as e:
            print(f"[WARN] Could not read RSCs in {fpath.name}: {e}")
            continue

        fc = check_fail_conditions(rsc_data["tables"])

        if fc["fail1"] and fc["fail2"]:
            fail_code = FAILB_CODE
        elif fc["fail1"]:
            fail_code = FAIL1_CODE
        elif fc["fail2"]:
            fail_code = FAIL2_CODE
        else:
            fail_code = None

        h = params["horizon"]
        a = params["alt"]
        base = {
            "horizon_lo":       h[0],
            "horizon_hi":       h[1],
            "alt_lo":           a[0],
            "alt_hi":           a[1],
            "horizon_width":    h[1] - h[0],
            "alt_hi_label":     a[1],
            "fail1":            fc["fail1"],
            "fail2":            fc["fail2"],
            "total_empty_rows": fc["total_empty_rows"],
        }

        # --- compute metrics for each data sheet ------------------------------
        for key, sheet_name in DATA_SHEETS.items():
            if fail_code is not None:
                metrics = {
                    "unicount":    fail_code,
                    "hist":        None,
                    "total_stars": fail_code,
                    "R0":          fail_code,
                    "R1":          fail_code,
                    "S":           fail_code,
                }
            else:
                try:
                    data    = read_sheet(fpath, sheet_name, has_metadata=False, cell_re=cell_re)
                    metrics = compute_all_metrics(data["tables"])
                except Exception as e:
                    print(f"[WARN] Could not read '{sheet_name}' in {fpath.name}: {e}")
                    metrics = {
                        "unicount":    np.nan,
                        "hist":        None,
                        "total_stars": np.nan,
                        "R0":          np.nan,
                        "R1":          np.nan,
                        "S":           np.nan,
                    }

            buffers[key].append({**base, **metrics})

    return {key: pd.DataFrame(rows) for key, rows in buffers.items()}




# ---------------------------------------------------------------------------
# 6.  Pivot to heatmap grid
# ---------------------------------------------------------------------------

def make_heatmap_grid(df: pd.DataFrame, metric: str = "unicount") -> pd.DataFrame:
    """
    Pivot a long-form sheet DataFrame into a 2-D grid:
        rows    = alt_hi_label    (sorted ascending)  → Y axis
        columns = horizon_width   (sorted ascending)  → X axis
        values  = metric column (default: "unicount")

    Available metrics: unicount, total_stars, R0, R1
    (hist is not plottable directly as a heatmap)

    Returns a DataFrame ready for plot_heatmap().
    """
    if metric not in df.columns:
        raise ValueError(f"metric '{metric}' not found. Available: {list(df.columns)}")
    h_order = sorted(df["horizon_width"].unique())
    a_order = sorted(df["alt_hi_label"].unique())
    grid = df.pivot(index="alt_hi_label", columns="horizon_width", values=metric)
    return grid.loc[a_order, h_order]


# ---------------------------------------------------------------------------
# 7.  Heatmap plot
# ---------------------------------------------------------------------------

def plot_heatmap(
    grid: pd.DataFrame,
    title: str = "Unique records",
    # fail1_color: str = "#444444",#"#bbbbbb",
    # fail2_color: str = "#242424",#"#888888",
    # failb_color: str = "#0D0D0D",#"#444444",
    fail1_color = "#969696",
    fail2_color = "#808080",
    failb_color = "#696969",
    cmap: str = "BrBG",
    figsize: tuple = (15, 5),
    ax=None,
    annotate: bool = False,
    clims: tuple = None,
    clabel: str = "",
    aspect: float = None,
):
    """
    Plot *grid* (output of make_heatmap_grid) as a colour heatmap.

    Axes:
        X = Maximum altitude       (alt_hi_label)
        Y = Horizon width [degrees] (horizon_width)

    Fail cells:
        F1 (fail1 only)   → fail1_color  (medium gray)
        F2 (fail2 only)   → fail2_color  (light gray)
        FB (both)         → failb_color  (dark gray)
    Passing cells use a diverging colormap centred on 0 (the reference value).

    Returns (fig, ax).
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches

    FAIL_CODES = {FAIL1_CODE, FAIL2_CODE, FAILB_CODE}
    is_fail = grid.map(lambda v: float(v) in FAIL_CODES)

    plot_data = grid.copy().astype(float)
    plot_data[is_fail] = np.nan

    auto_vmax = np.nanmax(np.abs(plot_data.values)) if not np.all(np.isnan(plot_data.values)) else 1
    
    (vmin, vmax) = (clims[0], clims[-1]) if clims is not None else (-auto_vmax, auto_vmax)
    if clims is not None and len(clims)==3:
        vcenter = clims[1]
    else:
        vcenter = (vmin + vmax)/2 # check if vcenter is specified     
          
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax) # always want divergence around a a number
    #norm = mcolors.Normalize(vmin=vmin vmax=vmax)

    # if clims is not None and clims[0] >= 0:
    #     norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    # elif clims is not None and clims[0] < 0:
    #     norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # else:
    #     norm = mcolors.Normalize(0, 1)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    cmap_obj = plt.get_cmap(cmap)
    #norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=(vmin+vmax)/2, vmax=vmax) if vmax > 0 else mcolors.Normalize(0, 1)

    im = ax.imshow(plot_data.values, aspect=aspect if aspect is not None else "auto",
                   cmap=cmap_obj, norm=norm, interpolation="nearest")

    _fail_styles = {
        FAIL1_CODE: (fail1_color, "F1"),
        FAIL2_CODE: (fail2_color, "F2"),
        FAILB_CODE: (failb_color, "F1+F2"),
    }

    if annotate:
        for (ri, ci), val in np.ndenumerate(grid.values):
            fval = float(val)
            if fval in _fail_styles:
                color, label = _fail_styles[fval]
                rect = mpatches.FancyBboxPatch(
                    (ci - 0.5, ri - 0.5), 1, 1,
                    boxstyle="square,pad=0", linewidth=0, facecolor=color,
                )
                ax.add_patch(rect)
                ax.text(ci, ri, label, ha="center", va="center",
                        fontsize=8, color="white", fontweight="bold")
            elif not np.isnan(fval):
                ax.text(ci, ri, str(int(fval)), ha="center", va="center",
                        fontsize=7, color="white" if abs(norm(fval) - 0.5) > 0.3 else "black")
    else:
        # still draw fail cell colors even without text
        for (ri, ci), val in np.ndenumerate(grid.values):
            fval = float(val)
            if fval in _fail_styles:
                color, _ = _fail_styles[fval]
                rect = mpatches.FancyBboxPatch(
                    (ci - 0.5, ri - 0.5), 1, 1,
                    boxstyle="square,pad=0", linewidth=0, facecolor=color,
                )
                ax.add_patch(rect)

    # X = horizon width, Y = maximum altitude
    ax.set_xticks(range(grid.shape[1]))
    ax.set_xticklabels([f"{c}°" for c in grid.columns], fontsize=9)
    ax.set_xlabel("Horizon width [degrees]", fontsize=10)

    # Y: invert so low altitude is at bottom, tick every 5
    alt_vals = sorted(grid.index)   # ascending
    ax.set_yticks(range(grid.shape[0]))
    ax.set_yticklabels([f"{r}°" if r % 5 == 0 else "" for r in alt_vals], fontsize=9)
    ax.set_ylabel("Maximum altitude [degrees]", fontsize=10)
    ax.invert_yaxis()

    ax.set_title(title, fontsize=12)
    fig.colorbar(im, ax=ax, label=clabel)

    f1_patch = mpatches.Patch(color=fail1_color, label="F1 = ≥10 empty rows total")
    f2_patch = mpatches.Patch(color=fail2_color, label="F2 = ≥2 empty rows in one table")
    fb_patch = mpatches.Patch(color=failb_color, label="F1+F2 = both conditions")
    ax.legend(handles=[f1_patch, f2_patch, fb_patch], loc="upper left",
              bbox_to_anchor=(0, -0.18), fontsize=8, frameon=False)

    fig.tight_layout()
    #plt.close(fig)
    return fig, ax



def plot_heatmap_grid(
    gridlist: list[pd.DataFrame],
    titles: list,
    fail1_color = "#969696",
    fail2_color = "#808080",
    failb_color = "#696969",
    cmap: str = "BrBG",
    figsize: tuple = (7, 7),
    annotate: bool = False,
    clims: tuple = None,
    clabel: str = "",
    aspect: float = None,
):
    """
    Plot *grid* (output of make_heatmap_grid) as a colour heatmap.

    Axes:
        X = Maximum altitude       (alt_hi_label)
        Y = Horizon width [degrees] (horizon_width)

    Fail cells:
        F1 (fail1 only)   → fail1_color  (medium gray)
        F2 (fail2 only)   → fail2_color  (light gray)
        FB (both)         → failb_color  (dark gray)
    Passing cells use a diverging colormap centred on 0 (the reference value).

    Returns (fig, ax).
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches

    # create subplot
    fig, axs = plt.subplots(nrows=1, ncols=len(gridlist), figsize=figsize, layout='constrained') #subplots

    cmap_obj = plt.get_cmap(cmap)
    #norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=(vmin+vmax)/2, vmax=vmax) if vmax > 0 else mcolors.Normalize(0, 1)

    # for marking failure cases
    _fail_styles = {
        FAIL1_CODE: (fail1_color, "F1"),
        FAIL2_CODE: (fail2_color, "F2"),
        FAILB_CODE: (failb_color, "F1+F2"),
    }

    FAIL_CODES = {FAIL1_CODE, FAIL2_CODE, FAILB_CODE}

    for i in range(0, len(gridlist)):
        # choose grid and axis
        grid = gridlist[i]
        ax = axs[i]
        
        # select data 
        plot_data = grid.copy().astype(float)
        is_fail = grid.map(lambda v: float(v) in FAIL_CODES)
        plot_data[is_fail] = np.nan
        
        # still draw fail cell colors even without text
        for (ri, ci), val in np.ndenumerate(grid.values):
            fval = float(val)
            if fval in _fail_styles:
                color, _ = _fail_styles[fval]
                rect = mpatches.FancyBboxPatch(
                    (ci - 0.5, ri - 0.5), 1, 1,
                    boxstyle="square,pad=0", linewidth=0, facecolor=color,
                )
                ax.add_patch(rect)

        # everything on same vmin-vmax norm
        auto_vmax = np.nanmax(np.abs(plot_data.values)) if not np.all(np.isnan(plot_data.values)) else 1
        (vmin, vmax) = (clims[0], clims[-1]) if clims is not None else (-auto_vmax, auto_vmax)
        if clims is not None and len(clims)==3:
            vcenter = clims[1]
        else:
            vcenter = (vmin + vmax)/2 # check if vcenter is specified      
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax) # always want divergence around a a number

        # start plotting 
        im = ax.imshow(plot_data.values, aspect=aspect if aspect is not None else "auto",
                    cmap=cmap_obj, norm=norm, interpolation="nearest")

        # X = horizon width, Y = maximum altitude
        ax.set_xticks(range(grid.shape[1]))
        ax.set_xticklabels([f"{c}°" for c in grid.columns], fontsize=9)
        ax.set_xlabel("Horizon width [degrees]", fontsize=10)

        # Y: invert so low altitude is at bottom, tick every 5
        alt_vals = sorted(grid.index)   # ascending
        ax.set_yticks(range(grid.shape[0]))
        ax.set_yticklabels([f"{r}°" if r % 5 == 0 else "" for r in alt_vals], fontsize=9)
        #ax.set_ylabel("Maximum altitude [degrees]", fontsize=10)
        ax.invert_yaxis()
        # force to be 1:1 aspect ratio
        ax.set_box_aspect(1)
        ax.set_title(titles[i], fontsize=12)
    # whole figure niceties
    axs[0].set_ylabel("Maximum altitude [degrees]", fontsize=10)
    fig.colorbar(im, ax=axs, label=clabel, shrink=0.84)

        # f1_patch = mpatches.Patch(color=fail1_color, label="F1 = ≥10 empty rows total")
        # f2_patch = mpatches.Patch(color=fail2_color, label="F2 = ≥2 empty rows in one table")
        # fb_patch = mpatches.Patch(color=failb_color, label="F1+F2 = both conditions")
        # ax.legend(handles=[f1_patch, f2_patch, fb_patch], loc="upper left",
        #         bbox_to_anchor=(0, -0.18), fontsize=8, frameon=False)

    #fig.tight_layout()
    #plt.close(fig)
    return fig, axs




# ---------------------------------------------------------------------------
# 8.  Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys, matplotlib.pyplot as plt

    if len(sys.argv) < 2:
        print("Usage: python schedule_analysis.py <folder>  [reference_value]")
        sys.exit(0)

    folder    = sys.argv[1]
    reference = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    sheets = scan_folder(folder, reference=reference)

    for key, df in sheets.items():
        grid    = make_heatmap_grid(df)
        fig, ax = plot_heatmap(grid, title=f"{DATA_SHEETS[key]} — unique records vs. ref ({reference})")
        out = f"heatmap_{key}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out}")
        plt.close(fig)
