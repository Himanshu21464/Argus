#!/usr/bin/env python3
"""
jwst_to_csv.py — convert MAST / Zenodo / pipeline-output JWST exoplanet
transmission spectra to the 3-column CSV the Argus loader eats.

USAGE
  jwst_to_csv.py INPUT [-o OUTPUT.csv] [--source 'citation tag']
  jwst_to_csv.py --zenodo 7388032 --member transit_spectra/FIREFLy_transit_spec.h5 -o w39.csv
  jwst_to_csv.py --mast-fits jw01366001001_x1dints.fits -o w39.csv

INPUT FORMATS RECOGNISED
  *.h5 / *.hdf5    — HDF5. Tries the FIREFLy / Tiberius / Eureka / Tshirt
                     column conventions; falls back to whatever
                     wavelength + transit_depth + uncertainty arrays
                     are present at the root.
  *.fits / *.fit   — FITS BINTABLE. Tries MAST x1d (WAVELENGTH / FLUX /
                     ERROR) and PHOENIX-style transit-spectrum tables.
  *.csv / *.ecsv   — already columnar; just normalises units.
  *.dat / *.txt    — whitespace-delimited; assumes 3 columns in
                     (wavelength, depth, sigma) order, with optional
                     `#` comment header.

UNIT NORMALISATION
  wavelength  : auto-detects μm / nm / Å. Default emit is μm.
  depth       : auto-detects fraction / ppm / percent.
                Default emit is dimensionless fraction (Rp/R*)^2.
  sigma_depth : same units as depth.

OUTPUT
  3-column CSV: `wavelength_um, transit_depth, sigma_depth`
  with `#`-prefixed header comments preserving the source citation
  + provenance (input filename, format, column mapping decisions).

NO RUNTIME ARGUS DEPENDENCY — pure-Python sidecar.
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import Optional, Tuple


# ─── Schema dictionaries — known reduction-pipeline column names. ─────
HDF5_WAVELENGTH_KEYS = (
    "wavelength", "wavelength_um", "wavelength_micron", "wave",
    "wl", "lambda", "WAVELENGTH",
)
HDF5_DEPTH_KEYS = (
    "transit_depth", "depth", "rprs2", "rp_rs2", "rprs_squared",
    "transit_depth_ppm", "TRANSIT_DEPTH",
)
HDF5_SIGMA_KEYS = (
    "transit_depth_uncertainty", "depth_uncertainty", "sigma",
    "depth_err", "transit_depth_err", "rprs2_err",
    "TRANSIT_DEPTH_UNCERTAINTY",
)

FITS_WAVELENGTH_KEYS = ("WAVELENGTH", "wavelength", "WAVE", "WL", "LAMBDA")
FITS_DEPTH_KEYS = (
    "TRANSIT_DEPTH", "DEPTH", "RPRS2", "TDEPTH", "FLUX",
)
FITS_SIGMA_KEYS = (
    "TRANSIT_DEPTH_UNCERTAINTY", "DEPTH_ERR", "ERROR", "SIGMA", "ERR",
)


def detect_format(path: Path) -> str:
    """Return one of: hdf5, fits, csv, dat."""
    sfx = path.suffix.lower()
    if sfx in (".h5", ".hdf5"):
        return "hdf5"
    if sfx in (".fits", ".fit"):
        return "fits"
    if sfx in (".csv", ".ecsv"):
        return "csv"
    if sfx in (".dat", ".txt"):
        return "dat"
    # Fall back to magic bytes.
    with open(path, "rb") as f:
        head = f.read(8)
    if head.startswith(b"\x89HDF"):
        return "hdf5"
    if head.startswith(b"SIMPLE  ="):
        return "fits"
    return "csv"


def first_match(names, keys):
    """Return the first name in `names` (case-insensitive) that matches any of `keys`."""
    lower = {n.lower(): n for n in names}
    for k in keys:
        if k.lower() in lower:
            return lower[k.lower()]
    return None


def detect_wavelength_unit(arr) -> Tuple[float, str]:
    """Return (factor_to_um, label). Heuristic: typical exoplanet spectra
    span ~0.3-30 μm. Same data in nm would be 300-30000; in Å, 3000-3e5."""
    import numpy as np
    median = float(np.median(arr))
    if 0.1 <= median <= 100.0:
        return 1.0, "μm"
    if 100.0 < median <= 1.0e5:
        return 1.0e-3, "nm"
    if median > 1.0e5:
        return 1.0e-4, "Å"
    return 1.0, "μm (assumed)"


def detect_depth_unit(arr) -> Tuple[float, str]:
    """Return (factor_to_fraction, label). Transit depth is typically
    1e-4 to 1e-1 as a fraction; 1e2 to 1e5 as ppm; 1e-2 to 1e1 as %."""
    import numpy as np
    median = float(np.median(arr))
    if 1.0e-5 <= median <= 0.1:
        return 1.0, "fraction"
    if median > 100.0:
        return 1.0e-6, "ppm"
    if 0.1 < median <= 100.0:
        return 1.0e-2, "percent"
    return 1.0, "fraction (assumed)"


# ─── HDF5 loader ──────────────────────────────────────────────────────
def load_hdf5(path: Path):
    try:
        import h5py
    except ImportError as e:
        sys.exit("jwst_to_csv.py: h5py is required to read HDF5 input — `pip install h5py`")
    with h5py.File(path, "r") as f:
        # Walk all datasets at root + one level deep so e.g. /spectrum/wavelength is found.
        candidates = {}
        def walk(name, obj):
            if hasattr(obj, "shape") and len(obj.shape) == 1:
                candidates[name] = obj[...]
        f.visititems(walk)
        # Also include root-level datasets (visititems skips them)
        for name in f.keys():
            obj = f[name]
            if hasattr(obj, "shape") and len(obj.shape) == 1:
                candidates[name] = obj[...]
        names = list(candidates.keys())
        wl_key  = first_match(names, HDF5_WAVELENGTH_KEYS)
        d_key   = first_match(names, HDF5_DEPTH_KEYS)
        sig_key = first_match(names, HDF5_SIGMA_KEYS)
        if not (wl_key and d_key and sig_key):
            sys.exit(
                "jwst_to_csv.py: could not find wavelength/depth/sigma datasets in HDF5.\n"
                f"  available 1-D datasets: {names}\n"
                f"  matched: wavelength={wl_key}, depth={d_key}, sigma={sig_key}"
            )
        return candidates[wl_key], candidates[d_key], candidates[sig_key], (wl_key, d_key, sig_key)


# ─── FITS loader ──────────────────────────────────────────────────────
def load_fits(path: Path):
    try:
        from astropy.io import fits
    except ImportError:
        sys.exit("jwst_to_csv.py: astropy is required to read FITS input — `pip install astropy`")
    with fits.open(path) as hdul:
        # Find a BINTABLE HDU.
        for hdu in hdul:
            if hdu.is_image or hdu.data is None:
                continue
            cols = hdu.columns.names
            wl_key  = first_match(cols, FITS_WAVELENGTH_KEYS)
            d_key   = first_match(cols, FITS_DEPTH_KEYS)
            sig_key = first_match(cols, FITS_SIGMA_KEYS)
            if wl_key and d_key and sig_key:
                return (
                    hdu.data[wl_key], hdu.data[d_key], hdu.data[sig_key],
                    (wl_key, d_key, sig_key),
                )
        sys.exit(
            "jwst_to_csv.py: no BINTABLE HDU with wavelength/depth/sigma columns in FITS.\n"
            f"  HDUs scanned: {[h.name for h in hdul]}"
        )


# ─── CSV / DAT loader ─────────────────────────────────────────────────
def load_csv_or_dat(path: Path, sep=None):
    import numpy as np
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in (line.split(",") if sep is None else line.split(sep)) if p.strip()]
            if sep is None and len(parts) < 3:
                # Maybe whitespace-delimited.
                parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"jwst_to_csv.py: no numeric rows in {path}")
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2], ("col0", "col1", "col2")


# ─── Zenodo fetch ─────────────────────────────────────────────────────
def fetch_zenodo(record_id: str, member: Optional[str], cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / f"zenodo_{record_id}.zip"
    if not zip_path.exists():
        url = f"https://zenodo.org/api/records/{record_id}"
        with urllib.request.urlopen(url, timeout=30) as r:
            import json
            meta = json.loads(r.read())
        files = meta.get("files", [])
        if not files:
            sys.exit(f"jwst_to_csv.py: zenodo record {record_id} has no files")
        # Take the first file — most exoplanet-data deposits ship one ZIP.
        link = files[0]["links"]["self"]
        print(f"[fetch] downloading {link} -> {zip_path}", file=sys.stderr)
        urllib.request.urlretrieve(link, zip_path)
    if member:
        with zipfile.ZipFile(zip_path) as zf:
            extracted = cache_dir / "extracted"
            extracted.mkdir(exist_ok=True)
            zf.extract(member, extracted)
            return extracted / member
    return zip_path


def main():
    ap = argparse.ArgumentParser(
        description="Convert JWST exoplanet spectrum (HDF5/FITS/CSV) to Argus CSV."
    )
    ap.add_argument("input", nargs="?", help="Input file (HDF5/FITS/CSV/DAT)")
    ap.add_argument("-o", "--output", default="-", help="Output CSV path (default: stdout)")
    ap.add_argument("--source", default="", help="Citation tag to embed in output header")
    ap.add_argument("--zenodo", help="Fetch from Zenodo by record ID (e.g. 7388032)")
    ap.add_argument("--member", help="If --zenodo input is a ZIP, extract this member")
    ap.add_argument("--cache-dir", default=".jwst_cache",
                    help="Where to cache downloaded archives (default: .jwst_cache)")
    args = ap.parse_args()

    # Resolve input path.
    if args.zenodo:
        in_path = fetch_zenodo(args.zenodo, args.member, Path(args.cache_dir))
    elif args.input:
        in_path = Path(args.input)
        if not in_path.exists():
            sys.exit(f"jwst_to_csv.py: input not found: {in_path}")
    else:
        ap.error("must give an INPUT file or --zenodo")

    fmt = detect_format(in_path)
    print(f"[detect] {in_path.name}  format={fmt}", file=sys.stderr)

    if fmt == "hdf5":
        wl, depth, sig, picked = load_hdf5(in_path)
    elif fmt == "fits":
        wl, depth, sig, picked = load_fits(in_path)
    elif fmt in ("csv", "dat"):
        wl, depth, sig, picked = load_csv_or_dat(in_path)
    else:
        sys.exit(f"jwst_to_csv.py: unknown format {fmt}")

    print(f"[map] wavelength={picked[0]}  depth={picked[1]}  sigma={picked[2]}", file=sys.stderr)
    print(f"[bins] {len(wl)} rows", file=sys.stderr)

    # Unit normalisation.
    wl_factor, wl_unit = detect_wavelength_unit(wl)
    d_factor,  d_unit  = detect_depth_unit(depth)
    print(f"[unit] wavelength {wl_unit} (factor {wl_factor}); "
          f"depth {d_unit} (factor {d_factor})", file=sys.stderr)

    wl_um = wl * wl_factor
    depth_frac = depth * d_factor
    sig_frac   = sig * d_factor       # sigma in same units as depth

    # Sort ascending in wavelength.
    import numpy as np
    order = np.argsort(wl_um)
    wl_um = wl_um[order]
    depth_frac = depth_frac[order]
    sig_frac = sig_frac[order]

    # Drop NaN / inf rows.
    mask = np.isfinite(wl_um) & np.isfinite(depth_frac) & np.isfinite(sig_frac)
    if (~mask).sum() > 0:
        print(f"[drop] {(~mask).sum()} rows with NaN/inf", file=sys.stderr)
    wl_um = wl_um[mask]; depth_frac = depth_frac[mask]; sig_frac = sig_frac[mask]

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        out.write(f"# JWST transmission spectrum — converted by jwst_to_csv.py\n")
        if args.source:
            out.write(f"# Source: {args.source}\n")
        out.write(f"# Input file: {in_path.name}  (format: {fmt})\n")
        out.write(f"# Column mapping: wl={picked[0]}, depth={picked[1]}, sigma={picked[2]}\n")
        out.write(f"# Unit conversion: wavelength {wl_unit} -> μm, depth {d_unit} -> fraction\n")
        out.write(f"# Bins: {len(wl_um)}\n")
        out.write(f"# Columns: wavelength_um, transit_depth, sigma_depth\n")
        for w, d, s in zip(wl_um, depth_frac, sig_frac):
            out.write(f"{w:.6f},{d:.7f},{s:.7f}\n")
    finally:
        if out is not sys.stdout:
            out.close()
            print(f"[write] {args.output}  ({len(wl_um)} bins)", file=sys.stderr)


if __name__ == "__main__":
    main()
