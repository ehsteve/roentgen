# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""

from pathlib import Path

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable, Table

try:
    from _version import version as __version__
except ImportError:
    __version__ = "0.0.0"  # Fallback for development mode

# roentgen specific configuration
# load some data files on import
_package_directory = Path(__file__).parent
_data_directory = _package_directory / "data"

elements_file = _data_directory / "elements.csv"
elements = QTable(ascii.read(elements_file, format="csv"))

elements["density"].unit = u.g / (u.cm**3)
elements["i"].unit = u.eV
elements["ionization energy"].unit = u.eV
elements["atomic mass"] = elements["z"] / elements["zovera"] * u.u
elements.add_index("z")
elements.add_index("symbol")

compounds_file = _data_directory / "compounds_mixtures.csv"
compounds = QTable(ascii.read(compounds_file, format="csv", fast_reader=False))
compounds["density"].unit = u.g / (u.cm**3)
compounds.add_index("symbol")

notation_translation = Table(
    ascii.read(
        _data_directory / "siegbahn_to_iupac.csv",
        format="csv",
        fast_reader=False,
    )
)
