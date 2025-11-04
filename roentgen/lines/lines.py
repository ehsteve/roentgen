"""A module to access atomic line emission, binding energies, and transmission edges."""

import numpy as np

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable

import roentgen
from roentgen.util import get_atomic_number, get_element_symbol

__all__ = ["get_lines", "get_edges", "emission_lines", "binding_energies"]

emission_lines = QTable(
    ascii.read(
        roentgen._data_directory / "emission_lines.csv",
        format="csv",
        fast_reader=False,
    )
)
# not sure why i need to fix this otherwise it is \ufenergy
# remove unit from column title to make it shorter
emission_lines.rename_column(emission_lines.colnames[0], "energy_ev")
emission_lines["energy_ev"].unit = u.eV
emission_lines.add_column(
    np.round(emission_lines["energy_ev"].to("keV"), 5), name="energy", index=0
)
emission_lines["width [eV]"].unit = u.eV
emission_lines.remove_column("energy_ev")
emission_lines.rename_column("width [eV]", "width")

emission_lines.add_index("energy")
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
        roentgen._data_directory / "electron_binding_energies.csv",
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
def get_lines(energy_low, energy_high, element=None, min_intensity: int = 0):
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

    min_intensity : int, optional
        Select only lines above or equal to a given intensity

    Returns
    -------
    line_list : `astropy.table.QTable`
    """
    bool_array = (emission_lines["energy"] < energy_high) * (
        emission_lines["energy"] > energy_low
    )
    if element is not None:
        bool_array *= emission_lines["z"] == get_atomic_number(element)

    if min_intensity > 0:
        bool_array *= emission_lines["intensity"] >= min_intensity

    return emission_lines[bool_array]


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
    for this_colname, this_element in zip(
        binding_energies.colnames, binding_energies.loc[z]
    ):
        if isinstance(this_element, u.Quantity) and (this_element.value > 0):
            columns.append(this_colname.split(" ")[0])
            energies.append(this_element)
    result = QTable(
        [energies, columns],
        names=("energy", "edge name"),
        meta={"element": f"{element} z={z}"},
    )

    return result
