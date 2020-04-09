import numpy as np

from astropy.table import QTable, Table
import astropy.units as u

import roentgen
from roentgen.absorption import get_atomic_number

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory

__all__ = [
    "get_lines"
]


def get_lines(energy_low, energy_high, element=None):
    """
    Retrieve all emission lines in an energy range.

    Parameters
    ----------
    energy_low : `astropy.units.Quantity`
        The low end of the energy range

    energy_high : `astropy.units.Quantity`
        The high end of the energy range

    element : str, optional
        Select only lines from a specific element

    Returns
    -------
    line_list : `astropy.table.QTable`
    """
    result = Table()

    em = roentgen.emission_intensities
    energies = em[em.colnames[0]]
    bool_array = (energies < energy_high) * (energies > energy_low)
    if np.sum(bool_array) > 0:
        ind = np.where(bool_array)[0]
        result = em.iloc[ind.min():ind.max()]

    if len(result) > 1 and element is not None:
        # check to see if any lines from selected element exist in energy range
        if np.sum(result['z'] == get_atomic_number(element)) > 0:
            result = result.loc['z', get_atomic_number(element)]

    return result
