# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""
import os

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable, Table

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"
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

emission_lines = QTable(
    ascii.read(
        os.path.join(_data_directory, "emission_lines.csv"),
        format="csv",
        fast_reader=False,
    )
)
# not sure why i need to fix this otherwise it is \ufenergy
emission_lines.rename_column(emission_lines.colnames[0], "energy")
emission_lines[emission_lines.colnames[0]].unit = u.eV
emission_lines.add_index(emission_lines.colnames[0])
emission_lines.add_index(emission_lines.colnames[1])
emission_lines.meta = {
    "source": "Center for X-ray Optics and Advanced Light Source, X-Ray Data Booklet Table 1-2",
    "publication date": "2009 October",
    "url": "https://xdb.lbl.gov/Section1/Table_1-3.pdf",
}
