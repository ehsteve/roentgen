import os

import numpy as np

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable

import roentgen
from roentgen.util import get_atomic_number, get_element_symbol

__all__ = ["get_lines", "get_edges", "emission_lines"]

emission_lines = QTable(
    ascii.read(
        os.path.join(roentgen._data_directory, "emission_lines.csv"),
        format="csv",
        fast_reader=False,
    )
)
# not sure why i need to fix this otherwise it is \ufenergy
emission_lines.rename_column(emission_lines.colnames[0], "energy")
emission_lines[emission_lines.colnames[0]].unit = u.eV
emission_lines.add_index(emission_lines.colnames[0])
emission_lines.add_index(emission_lines.colnames[1])
emission_lines.add_column(
    [get_element_symbol(int(z)) for z in emission_lines["z"]], name="symbol", index=2
)
emission_lines.meta = {
    "source": "Center for X-ray Optics and Advanced Light Source, X-Ray Data Booklet Table 1-3",
    "publication date": "2009 October",
    "url": "https://xdb.lbl.gov/Section1/Table_1-3.pdf",
}

binding_energies = QTable(
    ascii.read(
        os.path.join(roentgen._data_directory, "electron_binding_energies.csv"),
        format="csv",
        fast_reader=False,
    )
)

for this_col in binding_energies.colnames[2:]:
    binding_energies[this_col].unit = u.eV
binding_energies.add_index(binding_energies.colnames[0])
binding_energies.add_index(binding_energies.colnames[1])
binding_energies.meta = {
    "source": "Center for X-ray Optics and Advanced Light Source, X-Ray Data Booklet Table 1-1",
    "publication date": "2009 October",
    "url": "https://xdb.lbl.gov/Section1/Table_1-1.pdf",
}


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
    result = QTable()  # this is the default result

    energies = emission_lines[emission_lines.colnames[0]]
    bool_array = (energies < energy_high) * (energies > energy_low)
    if np.any(bool_array):
        result = emission_lines.loc[energy_low.to("eV").value : energy_high.to("eV").value]

    if len(result) > 1 and element is not None:
        # check to see if any lines from selected element exist in energy range
        if np.any(result["z"] == get_atomic_number(element)):
            result = result.loc["z", get_atomic_number(element)]

    return result


def get_edges(element):
    """
    Retrieve the absorption edges for a given element.
    Edges occur at the electron binding energies of the element.

    Parameters
    ----------
    element : str

    Returns
    -------
    edge_list : `astropy.table.QTable`
    """
    z = get_atomic_number(element)

    energies = []
    columns = []
    for this_colname, this_element in zip(binding_energies.colnames, binding_energies.loc[z]):
        if isinstance(this_element, u.Quantity) and (this_element.value > 0):
            columns.append(this_colname.split(" ")[0])
            energies.append(this_element)
    result = QTable(
        [energies, columns], names=("energy", "edge name"), meta={"element": f"{element} z={z}"}
    )

    return result
