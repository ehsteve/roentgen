"""
One-off script to extract the McMaster et al. (1969/1970) photoelectric,
coherent (Rayleigh), and incoherent (Compton) x-ray cross-section coefficients
used by SSW's xray package (resp_calc.pro and friends) into a CSV
file that can ship with roentgen.

Source: $SSW/packages/xray/dbase/xray_coef.fits, extension 1, as read by
$SSW/packages/xray/idl/coefdata.pro. This is the same data xsec.pro/xcross.pro
use for the "classic" (non-XCOM) cross-section calculation.

Reference:
McMaster W H, Del Grande N K, Mallett J H and Hubbell J H 1969, 1970
Compilation of x-ray cross sections, Lawrence Livermore National Laboratory
Report UCRL-50174.

Usage:
    python scripts/extract_xray_coef.py /usr/local/ssw/packages/xray/dbase/xray_coef.fits
"""

import io
import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

DEFAULT_SOURCE = "/usr/local/ssw/packages/xray/dbase/xray_coef.fits"
OUTPUT = Path(__file__).parent.parent / "roentgen" / "data" / "xray_coef.csv"

HEADER_COMMENT = """\
# Data Source: McMaster W H, Del Grande N K, Mallett J H and Hubbell J H
# 1969, 1970 Compilation of x-ray cross sections, Lawrence Livermore
# National Laboratory Report UCRL-50174.
# Extracted from SSW's $SSW/packages/xray/dbase/xray_coef.fits
# (as read by coefdata.pro), extension 1.
#
# at_den: density in g/cm^3
# conv: cross_section(cm2/g) * conv = cross_section(barns/atom)
# edge0..edge4: photoelectric energy sub-range boundaries, keV
# pe_i_j: photoelectric log/log cubic polynomial coefficients, sub-range i (0-3), power j (0-3)
# scc_0..scc_3: coherent (Rayleigh) scattering log/log cubic polynomial coefficients
# sci_0..sci_3: incoherent (Compton) scattering log/log cubic polynomial coefficients
#
# Valid energy range: 1-1000 keV. Used by xcross(energy, coef) as
# exp(sum_j coef[j] * ln(energy)**j) to give a cross-section in barns/atom.
"""


def main(source):
    with fits.open(source) as hdul:
        data = hdul[1].data

        table = Table()
        table["z"] = data["Z"].astype(int)
        table["at_den"] = data["AT_DEN"].astype(float)
        table["conv"] = data["CONV"].astype(float)
        for i in range(5):
            table[f"edge{i}"] = data["EDGE"][:, i].astype(float)
        # PE is stored (87, 4, 4) with dim order (sub_range, power) per FITS TDIM
        pe_array = data["PE"].astype(float)
        for i in range(4):
            for j in range(4):
                table[f"pe_{i}_{j}"] = pe_array[:, i, j]
        for j in range(4):
            table[f"scc_{j}"] = data["SCC"][:, j].astype(float)
        for j in range(4):
            table[f"sci_{j}"] = data["SCI"][:, j].astype(float)

        table.sort("z")

    buf = io.StringIO()
    table.write(buf, format="ascii.csv")

    with open(OUTPUT, "w") as f:
        f.write(HEADER_COMMENT)
        f.write(buf.getvalue())


if __name__ == "__main__":
    source = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SOURCE
    main(source)
    print(f"Wrote {OUTPUT}")
