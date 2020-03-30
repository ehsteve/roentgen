# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""
import os
from warnings import warn

from pkg_resources import get_distribution, DistributionNotFound

from astropy.tests.helper import TestRunner
from astropy.config.configuration import (
     update_default_config,
     ConfigurationDefaultMissingError,
     ConfigurationDefaultMissingWarning)

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

test = TestRunner.make_test_runner_in(os.path.dirname(__file__))

# roentgen specific configuration
# load some data files on import

from astropy.io import ascii
from astropy.table import QTable
import astropy.units as u

_package_directory = os.path.dirname(os.path.abspath(__file__))
_data_directory = os.path.abspath(os.path.join(_package_directory, 'data'))

elements = ascii.read(os.path.join(_data_directory, 'elements.csv'),
                      format='csv', fast_reader=False)

elemental_densities_file = os.path.join(_data_directory,
                                        'elements_densities.csv')
elemental_densities = QTable(ascii.read(elemental_densities_file))
elemental_densities['density'].unit = u.g / (u.cm ** 3)
# elemental_densities.add_index('symbol')

compounds_file = os.path.join(_data_directory, 'compounds_mixtures.csv')
compounds = QTable(ascii.read(compounds_file,
                              format='csv', fast_reader=False))
compounds['density'].unit = u.g / (u.cm ** 3)
# compounds.add_index('symbol')

emission_energies_file = os.path.join(_data_directory,
                                      'emission_energies.csv')
emission_energies = QTable(ascii.read(emission_energies_file))
for colname in emission_energies.colnames[2:]:
    emission_energies[colname].unit = u.eV
