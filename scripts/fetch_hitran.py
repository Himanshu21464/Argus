#!/usr/bin/env python3
"""
fetch_hitran.py — download HITRAN line lists into the local Argus opacity
cache, using HITRAN's official HAPI (HITRAN Application Programming
Interface). Wraps anonymous queries; no login token required.

Default cache: ~/.argus/opacity/ (override with --cache-dir).

USAGE
  scripts/fetch_hitran.py H2O                          # default 1500-20000 cm^-1 (JWST PRISM)
  scripts/fetch_hitran.py H2O CO2 CO Na                # multi-molecule
  scripts/fetch_hitran.py H2O --range 6800 7600        # custom wavenumber window
  scripts/fetch_hitran.py --list                       # show what's already cached

OUTPUT
  Each species lands in {cache_dir}/{NAME}.par — the standard 160-char
  HITRAN .par layout that argus::Hitran::load_file consumes directly.

CITATION
  Uses HAPI (Kochanov et al. 2016, JQSRT 177, 15-30; DOI 10.1016/j.jqsrt.2016.03.005)
  Data sourced from HITRAN-2020 (Gordon+ 2022, JQSRT 277, 107949).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

# HITRAN molecule IDs (subset relevant for exoplanet atmospheres)
MOLECULE_ID = {
    "H2O":  (1,  1),     # 1H2-16O
    "CO2":  (2,  1),     # 12C-16O2
    "O3":   (3,  1),
    "N2O":  (4,  1),
    "CO":   (5,  1),     # 12C-16O
    "CH4":  (6,  1),     # 12CH4
    "O2":   (7,  1),
    "NO":   (8,  1),
    "SO2":  (9,  1),
    "NH3":  (11, 1),
    "HCN":  (23, 1),
    "C2H2": (26, 1),
    # NOTE: Atomic alkalis (Na, K) are NOT in HITRAN — those come from
    # VALD or NIST atomic spectra databases. WASP-39b fits include them
    # via the bundled Na D doublet in argus/test_data.hpp.
}

# Reasonable per-molecule default ranges for JWST-PRISM-band exoplanet work
# (cm^-1; trim huge ranges where HITRAN has 10^6+ lines we don't need).
DEFAULT_RANGE = {
    "H2O":  (1500, 20000),   # 0.5-6.7 μm — covers PRISM + visible
    "CO2":  (1800,  8000),
    "CO":   (1900,  4400),
    "CH4":  (1000,  6300),
    "NH3":  (700,   5300),
}


def main():
    ap = argparse.ArgumentParser(
        description="Download HITRAN line lists into the Argus opacity cache."
    )
    ap.add_argument("species", nargs="*", default=[],
                    help="Species names (H2O, CO2, CO, CH4, NH3, Na, K, …).")
    ap.add_argument("--range", type=float, nargs=2, metavar=("LO_CM", "HI_CM"),
                    help="Override the default wavenumber range (applies to all species).")
    ap.add_argument("--cache-dir", default=None,
                    help=f"Cache directory (default: ~/.argus/opacity)")
    ap.add_argument("--list", action="store_true", help="List cached files and exit.")
    ap.add_argument("--clean", action="store_true",
                    help="Remove the cache before fetching.")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir or
                     os.environ.get("ARGUS_OPACITY_CACHE",
                                    Path.home() / ".argus" / "opacity"))

    if args.list:
        cache_dir.mkdir(parents=True, exist_ok=True)
        files = sorted(cache_dir.glob("*.par"))
        if not files:
            print(f"(empty cache at {cache_dir})")
            return 0
        print(f"Cached HITRAN .par files in {cache_dir}:")
        for f in files:
            print(f"  {f.name:<20} {f.stat().st_size/1024:.1f} KB")
        return 0

    if args.clean and cache_dir.exists():
        print(f"[clean] removing {cache_dir}", file=sys.stderr)
        shutil.rmtree(cache_dir)

    if not args.species:
        ap.error("must give species (or --list / --clean)")

    cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"[cache] {cache_dir}", file=sys.stderr)

    try:
        import hapi  # type: ignore
    except ImportError:
        sys.exit(
            "fetch_hitran.py: HAPI is required.  Install with:\n"
            "    pip install --user hitran-api\n"
            "(See https://hitran.org/hapi/ for full docs.)"
        )

    # HAPI keeps its own local DB folder; point it at our cache.
    hapi.db_begin(str(cache_dir))

    overall_ok = True
    for spec in args.species:
        if spec not in MOLECULE_ID:
            print(f"[skip] unknown species '{spec}' — known: "
                  f"{', '.join(MOLECULE_ID)}", file=sys.stderr)
            overall_ok = False
            continue
        mol_id, iso_id = MOLECULE_ID[spec]
        lo, hi = (args.range if args.range else DEFAULT_RANGE.get(spec, (500, 25000)))
        print(f"[fetch] {spec}  molecule_id={mol_id}  iso={iso_id}  "
              f"range={lo:.0f}-{hi:.0f} cm^-1", file=sys.stderr)
        try:
            hapi.fetch(spec, mol_id, iso_id, lo, hi)
        except Exception as e:
            print(f"[fail] {spec}: {type(e).__name__}: {e}", file=sys.stderr)
            overall_ok = False
            continue

        # HAPI saves the line data to {cache_dir}/{spec}.data + .header.
        # Argus reads the canonical .par; rename for clarity.
        data_path = cache_dir / f"{spec}.data"
        par_path  = cache_dir / f"{spec}.par"
        if data_path.exists():
            shutil.copy(data_path, par_path)
            n_lines = sum(1 for _ in open(par_path))
            print(f"[ok]   {spec}: {n_lines} lines  → {par_path}",
                  file=sys.stderr)
        else:
            print(f"[fail] {spec}: HAPI produced no .data file",
                  file=sys.stderr)
            overall_ok = False

    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
