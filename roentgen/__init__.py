# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A Python package for the quantitative analysis of the interaction of energetic
photons with matter (x-rays and gamma-rays).
"""
import os
from astropy.io import ascii
from astropy.table import QTable
import astropy.units as u
# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
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
    #compounds.add_index('symbol')

    emission_energies_file = os.path.join(_data_directory,
                                          'emission_energies.csv')
    emission_energies = QTable(ascii.read())
    for colname in emission_energies.colnames[2:]:
        emission_energies[colname].unit = u.eV

