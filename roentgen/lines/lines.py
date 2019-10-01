from astropy.table import QTable
import astropy.units as u

import roentgen
from roentgen.absorption import get_element_symbol, is_an_element

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory


def get_lines_for_element(element_str):
    """
    Retrieve all emission lines for an element.

    Parameters
    ----------
    element_str : str
        The name of the element (e.g. Zn, silicon)

    Returns
    -------
    line_list : `astropy.table.QTable`
    """
    if is_an_element(element_str):
        symbol_str = get_element_symbol(element_str)
        index = list(roentgen.emission_energies['symbol']).index(symbol_str)
        row = roentgen.emission_energies[index]
        cols = list(roentgen.emission_energies[1].columns[2:])
        line_table = QTable()
        energies = []
        for this_col in cols:
            energies.append(row[this_col])
        energies = u.Quantity(energies)
        line_table['line'] = cols
        line_table['energy'] = u.Quantity(energies)
        mask = energies.value == 0
        line_table['energy'].mask = mask
        return line_table
    else:
        ValueError(f"Element not recognized.")


def get_lines(energy_low, energy_high):
    """
    Retrieve all emission lines in a energy range.

    Parameters
    ----------
    energy_low : `astropy.units.Quantity`
        The low end of the energy range

    energy_high : `astropy.units.Quantity`
        The high end of the energy range

    Returns
    -------
    line_list : `astropy.table.QTable`
    """
    line_cols = list(roentgen.emission_energies[0].columns[2:])

    line_energy = []
    element_name = []
    line_type = []

    for z in roentgen.emission_energies:
        this_element = z['symbol']
        energies = []
        for this_col in line_cols:
            if (z[this_col] > energy_low) and (z[this_col] < energy_high):
                line_energy.append(z[this_col])
                element_name.append(this_element)
                line_type.append(this_col)

    if len(line_energy) > 0:
        result = QTable()
        result['energy'] = u.Quantity(line_energy)
        result['element'] = element_name
        result['type'] = line_type
    else:
        result = None

    return result
