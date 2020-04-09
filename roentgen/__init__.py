# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""
import os
from warnings import warn

from pkg_resources import get_distribution, DistributionNotFound

#from astropy.tests.helper import TestRunner
#from astropy.config.configuration import (
#     update_default_config,
#     ConfigurationDefaultMissingError,
#     ConfigurationDefaultMissingWarning)

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:  # package is not installed
    __version__ = "unknown"

# add these here so we only need to cleanup the namespace at the end
config_dir = None

if not os.environ.get('ASTROPY_SKIP_CONFIG_UPDATE', False):
    config_dir = os.path.dirname(__file__)
    config_template = os.path.join(config_dir, __package__ + ".cfg")
    if os.path.isfile(config_template):
        try:
            update_default_config(
                __package__, config_dir, version=__version__)
        except TypeError as orig_error:
            try:
                update_default_config(__package__, config_dir)
            except ConfigurationDefaultMissingError as e:
                wmsg = (e.args[0] +
                        " Cannot install default profile. If you are "
                        "importing from source, this is expected.")
                warn(ConfigurationDefaultMissingWarning(wmsg))
                del e
            except Exception:
                raise orig_error

#test = TestRunner.make_test_runner_in(os.path.dirname(__file__))

# roentgen specific configuration
# load some data files on import

from astropy.io import ascii
from astropy.table import QTable, Table
import astropy.units as u

_package_directory = os.path.dirname(os.path.abspath(__file__))
_data_directory = os.path.abspath(os.path.join(_package_directory, 'data'))

elements_file = os.path.join(_data_directory, 'elements.csv')
elements = QTable(ascii.read(elements_file, format='csv'))

elements['density'].unit = u.g / (u.cm ** 3)
elements['i'].unit = u.eV
elements['ionization energy'].unit = u.eV
elements['atomic mass'] = elements['z'] / elements['zovera'] * u.u
elements.add_index('z')

compounds_file = os.path.join(_data_directory, 'compounds_mixtures.csv')
compounds = QTable(ascii.read(compounds_file, format='csv', fast_reader=False))
compounds['density'].unit = u.g / (u.cm ** 3)
# compounds.add_index('symbol')

notation_translation = Table(ascii.read(os.path.join(_data_directory, 'siegbahn_to_iupac.csv'),
                                        format='csv', fast_reader=False))

emission_energies_file = os.path.join(_data_directory, 'emission_energies.csv')
emission_energies = QTable(ascii.read(emission_energies_file, fill_values=('', '-1')))
for colname in emission_energies.colnames[1:]:
    emission_energies[colname].unit = u.eV
emission_energies.meta = {"source": "Center for X-ray Optics and Advanced Light Source, X-Ray Data Booklet Table 1-2",
                          "publication date": "2009 October",
                          "url": "https://xdb.lbl.gov/Section1/Sec_1-2.html"}
emission_energies.add_index('z')

# TODO: add masking of quantity columns once this is made possible in astropy

emission_intensities = QTable(ascii.read(os.path.join(_data_directory, 'emission_intensities.csv'), format='csv', fast_reader=False))
# not sure why i need to fix this otherwise it is \ufenergy
emission_intensities.rename_column(emission_intensities.colnames[0], 'energy')
emission_intensities[emission_intensities.colnames[0]].unit = u.eV
emission_intensities.add_index(emission_intensities.colnames[0])
emission_intensities.add_index(emission_intensities.colnames[1])
emission_intensities.meta = {"source": "Center for X-ray Optics and Advanced Light Source, X-Ray Data Booklet Table 1-2",
                             "publication date": "2009 October",
                             "url": "https://xdb.lbl.gov/Section1/Table_1-3.pdf"}