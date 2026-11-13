"""
One-off script to extract the incoherent scattering function S(q,Z) into a
CSV file that ships with roentgen.

Source: J.H. Hubbell, Wm.J. Veigele, E.A. Briggs, R.T. Brown, D.T. Cromer
and R.J. Howerton, J. Phys. Chem. Ref. Data 4, 471 (1975); erratum 6, 615
(1977), as distributed with the ``xraylib`` package (``SF_Compt``). xraylib is
only needed to run this script, not to use the resulting data.

q = sin(theta/2) / wavelength, in 1/Angstrom. xraylib has S(q,Z) for z=1-98
and refuses q below ~1e-3 (there S goes as q**2, which `roentgen.response.compton`
extrapolates).

Usage:
    python scripts/extract_scattering_function.py
"""

from pathlib import Path

import numpy as np
import xraylib

OUTPUT = Path(__file__).parent.parent / "roentgen" / "data" / "incoherent_scattering_function.csv"

Q = np.geomspace(1.1e-3, 1e3, 300)
ZS = range(1, 99)

HEADER_COMMENT = """\
# Data Source: Hubbell et al., J. Phys. Chem. Ref. Data 4, 471 (1975); erratum 6, 615 (1977),
# sampled from xraylib's SF_Compt (https://github.com/tschoonj/xraylib) by
# scripts/extract_scattering_function.py
#
# q: sin(theta/2)/wavelength, 1/Angstrom
# zNN: incoherent scattering function S(q, Z=NN), dimensionless (-> Z as q -> infinity)
"""

if __name__ == "__main__":
    columns = [[xraylib.SF_Compt(z, q) for q in Q] for z in ZS]
    with open(OUTPUT, "w") as f:
        f.write(HEADER_COMMENT)
        f.write("q," + ",".join(f"z{z:02d}" for z in ZS) + "\n")
        for i, q in enumerate(Q):
            f.write(f"{q:.6g}," + ",".join(f"{col[i]:.6g}" for col in columns) + "\n")
