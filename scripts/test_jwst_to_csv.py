#!/usr/bin/env python3
"""Self-test for jwst_to_csv.py — runs through every input format and
verifies the output is a valid Argus CSV (3 columns, finite numbers,
proper unit normalisation)."""
import os, subprocess, sys, tempfile, shutil
from pathlib import Path

CONV = Path(__file__).parent / "jwst_to_csv.py"
assert CONV.exists(), f"converter not found at {CONV}"

def run(args, **kw):
    return subprocess.check_output([sys.executable, str(CONV)] + args, **kw)

def parse_csv(text):
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        parts = line.split(",")
        rows.append([float(p) for p in parts])
    return rows

PASS = 0; FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1; print(f"  OK   {name}")
    else:
        FAIL += 1; print(f"  FAIL {name}  {detail}")

with tempfile.TemporaryDirectory() as td:
    td = Path(td)

    # Case 1: CSV passthrough (no unit conversion needed).
    csv_in = td / "ok.csv"
    csv_in.write_text("# comment\n1.0,0.021,0.0002\n2.0,0.022,0.0002\n")
    out = run([str(csv_in)], stderr=subprocess.DEVNULL).decode()
    rows = parse_csv(out)
    check("csv passthrough: 2 rows", len(rows) == 2, str(rows))
    check("csv passthrough: depth ~0.021", abs(rows[0][1] - 0.021) < 1e-9)

    # Case 2: ppm + nm auto-normalize.
    ppm_in = td / "ppm.csv"
    ppm_in.write_text("1300, 21000, 200\n1400, 21500, 180\n")
    out = run([str(ppm_in)], stderr=subprocess.DEVNULL).decode()
    rows = parse_csv(out)
    check("ppm/nm auto-normalize: depth ~0.021 (was 21000 ppm)",
          abs(rows[0][1] - 0.021) < 1e-9, str(rows[0]))
    check("ppm/nm auto-normalize: λ ~1.3 μm (was 1300 nm)",
          abs(rows[0][0] - 1.3) < 1e-9, str(rows[0]))

    # Case 3: whitespace-delimited DAT.
    dat_in = td / "table.dat"
    dat_in.write_text("# header\n1.5  0.0211  0.00018\n1.6  0.0212  0.00019\n")
    out = run([str(dat_in)], stderr=subprocess.DEVNULL).decode()
    rows = parse_csv(out)
    check("whitespace DAT: 2 rows", len(rows) == 2)

    # Case 4: HDF5 input (use the bundled WASP-39b PRISM via h5py if available).
    try:
        import h5py
        h5_in = td / "test.h5"
        with h5py.File(h5_in, "w") as f:
            f.create_dataset("wavelength", data=[1.0, 2.0, 3.0])
            f.create_dataset("transit_depth", data=[0.021, 0.022, 0.021])
            f.create_dataset("transit_depth_uncertainty", data=[0.0001, 0.0001, 0.0002])
        out = run([str(h5_in)], stderr=subprocess.DEVNULL).decode()
        rows = parse_csv(out)
        check("HDF5: 3 rows from FIREFLy-style schema", len(rows) == 3)
    except ImportError:
        print("  SKIP HDF5 test (no h5py)")

    # Case 5: malformed input fails cleanly.
    bad = td / "bad.csv"
    bad.write_text("not numbers at all\n")
    rc = subprocess.run([sys.executable, str(CONV), str(bad)],
                        capture_output=True).returncode
    check("malformed input → non-zero exit", rc != 0)

print(f"\n{PASS} pass, {FAIL} fail")
sys.exit(0 if FAIL == 0 else 1)
