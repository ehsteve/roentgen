import numpy as np

import astropy.units as u
from astropy.table import Table

import roentgen
from roentgen.util import get_atomic_number

__all__ = ["get_lines"]


@u.quantity_input(energy_low=u.keV, energy_high=u.keV, equivalencies=u.spectral())
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
    result = Table()  # this is the default result

    em = roentgen.emission_lines
    energies = em[em.colnames[0]]
    bool_array = (energies < energy_high) * (energies > energy_low)
    if np.any(bool_array):
        result = em.loc[energy_low.to("eV").value : energy_high.to("eV").value]

    if len(result) > 1 and element is not None:
        # check to see if any lines from selected element exist in energy range
        if np.any(result["z"] == get_atomic_number(element)):
            result = result.loc["z", get_atomic_number(element)]

    return result
