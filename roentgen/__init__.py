# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""
import os
from importlib.metadata import version, PackageNotFoundError

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable, Table

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass

__all__ = []

# roentgen specific configuration
# load some data files on import

_package_directory = os.path.dirname(os.path.abspath(__file__))
_data_directory = os.path.abspath(os.path.join(_package_directory, "data"))

elements_file = os.path.join(_data_directory, "elements.csv")
elements = QTable(ascii.read(elements_file, format="csv"))

elements["density"].unit = u.g / (u.cm**3)
elements["i"].unit = u.eV
elements["ionization energy"].unit = u.eV
elements["atomic mass"] = elements["z"] / elements["zovera"] * u.u
elements.add_index("z")

compounds_file = os.path.join(_data_directory, "compounds_mixtures.csv")
compounds = QTable(ascii.read(compounds_file, format="csv", fast_reader=False))
compounds["density"].unit = u.g / (u.cm**3)
compounds.add_index("symbol")

notation_translation = Table(
    ascii.read(
        os.path.join(_data_directory, "siegbahn_to_iupac.csv"),
        format="csv",
        fast_reader=False,
    )
)
