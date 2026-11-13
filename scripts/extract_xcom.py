"""
One-off script to extract per-element photon cross-section tables from
SolarSoftWare's bundled NIST XCOM data into CSV files that ship with
roentgen.

Source: $SSW/packages/xray/dbase/xcom.fits, as read by
$SSW/packages/xray/idl/xcom__define.pro. Unlike the McMaster 1969/1970
coefficients (see extract_xray_coef.py), this is current NIST XCOM data,
valid from 1 keV to 1e8 keV (100 GeV), and additionally separates out
coherent (Rayleigh) scattering and pair production.

xcom.fits bundles 126 materials (elements and a few named compounds/crystals,
e.g. BGO, ZnSe); which entries are elements, and at what index, is given by
xcom__define.pro's ``mat_vs_z`` method, whose ``qord`` lookup table is
reproduced here verbatim (index by z-1, -1 means "no table", matching
IDL's own "No table for z = 91" message).

Usage:
    python scripts/extract_xcom.py /usr/local/ssw/packages/xray/dbase/xcom.fits
"""

import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

DEFAULT_SOURCE = "/usr/local/ssw/packages/xray/dbase/xcom.fits"
OUTPUT_DIR = Path(__file__).parent.parent / "roentgen" / "data" / "xcom"

# xcom__define.pro's mat_vs_z: qord[z - 1] -> index into the 126-material
# list (0-based). -1 means no table for that z (z=91).
QORD = [
    54, 50, 64, 10, 13, 19, 76, 78, 38, 73, 101, 67, 2, 99, 82, 103, 28, 5, 86, 17,
    96, 112, 117, 29, 68, 60, 30, 74, 31, 122, 42, 46, 6, 97, 15, 61, 93, 102, 121, 123,
    75, 70, 105, 94, 92, 79, 100, 16, 55, 111, 4, 106, 58, 118, 27, 9, 62, 26, 87, 72,
    88, 95, 36, 41, 107, 34, 53, 35, 110, 120, 66, 49, 104, 115, 91, 77, 59, 83, 47, 69,
    108, 63, 12, 85, 7, 90, 39, 89, 0, 109, -1, 116,
]  # fmt: skip

HEADER_COMMENT = """\
# Data Source: NIST XCOM (https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html)
# Extracted from SSW's $SSW/packages/xray/dbase/xcom.fits (as read by
# xcom__define.pro), one extension per element (see xcom__define.pro's
# mat_vs_z / QORD in scripts/extract_xcom.py for the z -> extension mapping).
#
# energy_kev: photon energy, keV
# coherent, incoherent, photoelectric, pair_nuclear, pair_electron,
# total_with_coherent, total_without_coherent: cm2/g
#
# Valid from 1 keV to 1e8 keV. K/L/M-edges appear as duplicate energy rows
# (low-energy-side then high-energy-side value), same convention as
# roentgen's other NIST-sourced per-element data files.
"""

COLUMNS = [
    "energy_kev",
    "coherent",
    "incoherent",
    "photoelectric",
    "pair_nuclear",
    "pair_electron",
    "total_with_coherent",
    "total_without_coherent",
]


def main(source):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with fits.open(source) as hdul:
        for z, index in enumerate(QORD, start=1):
            if index == -1:
                continue
            data = hdul[1 + index].data
            name = data["XSUBS"][0].tobytes().split(b"\x00")[0].decode("ascii")
            xdat = data["XDAT"][0]  # shape (8, n_energies): see COLUMNS order

            table = Table()
            table["energy_kev"] = xdat[0] * 1000.0
            for i, col in enumerate(COLUMNS[1:], start=1):
                table[col] = xdat[i]

            out_path = OUTPUT_DIR / f"z{z:02d}.csv"
            with open(out_path, "w") as f:
                f.write(HEADER_COMMENT)
                f.write(f"# Element: {name} (z={z})\n")
            with open(out_path, "a") as f:
                table.write(f, format="ascii.csv")
            print(f"z={z:2d} {name:12s} -> {out_path.name} ({len(table)} points)")


if __name__ == "__main__":
    source = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SOURCE
    main(source)
    print(f"Wrote element CSVs to {OUTPUT_DIR}")
