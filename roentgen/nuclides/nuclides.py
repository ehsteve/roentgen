"""Support for the emission of radiation from radionuclides"""

from pathlib import Path

import numpy as np

from astropy.table import QTable, Table, vstack
import astropy.units as u

import roentgen
from roentgen import elements

__all__ = ["Nuclide"]

_lara_directory = Path(roentgen._data_directory) / "lara"
_lara_files = list(_lara_directory.glob("*.txt"))

nuclides_list = Table(data={"filename": [this_file.name for this_file in _lara_files]})
nuclides_list["symbol"] = [this_row["filename"].split("-")[0] for this_row in nuclides_list]
nuclides_list["mass_number"] = [
    int(this_row["filename"].split("-")[1].split("_")[0]) for this_row in nuclides_list
]
nuclides_list["name"] = [
    elements.loc["symbol", this_row["symbol"]]["name"] for this_row in nuclides_list
]
nuclides_list.add_index("symbol")
nuclides_list.add_index("mass_number")


class Nuclide(object):
    """An object to support radionuclides that emit x-ray radiation.

    Parameters
    ----------
    material_str : str
        A string representation of the material which includes an element symbol
        (e.g. Si), an element name (e.g. Silicon).
        For all supported radionuclides see :download:`lara <../../roentgen/data/lara/>`.
    mass_number : int
        The mass number of the radionuclide (e.g. 55 for Fe-55)

    Attributes
    ----------
    symbols : `list`
        A list of material symbol
    material_names : `list`
        A list of material names
    name : `str`
        A name for the material

    Methods
    -------
    get_lines(energy_range)
        X-ray lines in the energy range.

    Examples
    --------
    >>> from roentgen.nuclides.material import Nuclide
    >>> import astropy.units as u
    >>> fe55 = Nuclide('fe', 55)
    >>> fe55.get_lines(1 * u.kev, 10 * u.keV)
    """

    def __init__(self, element: str, mass_number: int):
        filename = get_lara_file(element, mass_number)
        self._line_tables = read_lara_tables(filename)
        self.lines = vstack(self._line_tables)
        self.meta = read_lara_header(filename)
        self.lines.sort("energy")
        self.name = self.meta["Nuclide"]
        self.element = self.meta["Element"]
        self.half_life = self.meta["Half-life (s)"]

    def __str__(self):
        return f"{self._text_summary()}{self.lines.__repr__()}"

    def __repr__(self):
        return f"{object.__repr__(self)}\n{self}"

    def _text_summary(self):
        num_lines = len(self.lines)
        result = f"Nuclide {self.name}, half_life={self.half_life} - ({num_lines:,} lines)\n"
        if len(self._line_tables) == 1:
            decay_chain = f"{self.name}->{self._line_tables[0]['origin'][-1].lstrip()}"
        else:
            decay_chain = "->".join([this_table["parent"][0] for this_table in self._line_tables])
        result += f"Decay chain: {decay_chain}\n"
        return result

    def get_lines(self, energy_low, energy_high, min_intensity: float = 0.0):
        """Returns a list of all emission lines in the energy range."""
        bool_array = (self.lines["energy"] < energy_high) * (self.lines["energy"] > energy_low)
        if min_intensity > 0:
            bool_array *= self.lines["intensity"] > min_intensity
        return self.lines[bool_array]


def get_lara_file(element: str, mass_number: int):
    these_nuclides = nuclides_list.loc["symbol", element]
    if isinstance(these_nuclides, Table):  # multiple isotopes exist
        this_nuclide = these_nuclides.loc["mass_number", mass_number]
    else:
        this_nuclide = these_nuclides
    return _lara_directory / str(this_nuclide["filename"])


def read_lara_tables(file_path: str | Path) -> list:
    """Return a table of all emissions from all origins.

    Returns
    -------
    result : list"""
    if isinstance(file_path, str):
        file_path = Path(file_path)
    with open(file_path, "r") as fp:
        lines = [line.rstrip() for line in fp]
    # find the emission tables
    table_line_index = []
    for i, this_line in enumerate(lines):
        if this_line.count("---------") > 1:
            table_line_index.append(i)
    result = []
    for j, this_index in enumerate(table_line_index):
        this_table = QTable()
        energy = []
        intensity = []
        origin = []
        # get the parent nuclide, may be self
        table_separator = lines[this_index]
        if table_separator.count(" "):
            parent = table_separator.split(" ")[1]
        else:
            parent = file_path.name.split("_")[0]
        skip_num = 2  # first table has the header line
        if j > 0:
            skip_num = 1
        for this_line in lines[this_index + skip_num :]:
            tokens = this_line.split(";")
            if len(tokens) > 1:
                energy.append(float(tokens[0]))
                intensity.append(float(tokens[2]))
                origin.append(tokens[5])
            else:
                break
        energy = u.Quantity(energy, "keV")
        this_table["energy"] = energy
        this_table["intensity"] = intensity
        this_table["origin"] = np.array(origin, dtype="<U7")
        this_table["parent"] = np.array([parent] * len(origin), dtype="<U7")
        this_table.add_index("origin")
        result.append(this_table)
    return result


def read_lara_header(file_path: str | Path) -> dict:
    """Return a table of all emissions from all origins.

    Returns
    -------
    result : list"""
    if isinstance(file_path, str):
        file_path = Path(file_path)
    with open(file_path, "r") as fp:
        lines = fp.readlines()[:12]
        result = {}
    for this_line in lines:
        tokens = this_line.rstrip().split(";")
        if len(tokens) > 1:
            key = tokens[0].rstrip()
            value = tokens[1].rstrip().lstrip()
            try:
                value = float(value)
            except ValueError:
                pass
            if key.count("Q+"):
                value = value * u.keV
            elif key.count("Half-life (a)"):
                value = value * u.year
            elif key.count("Half-life (s)"):
                print(value)
                value = value * u.s
            elif key.count("Decay"):
                value = value / u.s
            elif key.count("activity"):
                value = value * u.Bq / u.g
            result.update({key: value})
    return result
