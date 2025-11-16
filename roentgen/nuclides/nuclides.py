"""A module to provide access to x-ray and gamma-ray radiation from radionuclides"""

from pathlib import Path

import numpy as np

from astropy.table import QTable, vstack
import astropy.units as u
from astropy.io import ascii

import roentgen

__all__ = ["Nuclide", "get_nuclide_mass_numbers", "nuclides_list"]

_lara_directory = Path(roentgen._data_directory) / "lara"

nuclides_list = QTable(
    ascii.read(
        Path(roentgen._data_directory) / "nuclides_list.csv",
        format="csv",
        fast_reader=False,
    )
)
nuclides_list["half_life"] = nuclides_list["half_life [year]"]
nuclides_list.remove_column("half_life [year]")
nuclides_list["half_life"].unit = u.yr

nuclides_list.add_index("symbol")
nuclides_list.add_index("mass_number")


class Nuclide(object):
    """An object to support radionuclides that emit x-ray radiation.

    Parameters
    ----------
    material_str : str
        A string representation of the material which includes an element symbol
        (e.g. Si), an element name (e.g. Silicon).
        For all supported radionuclides see :download:`lara <../../roentgen/data/nuclides_list.csv/>`.
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
    >>> from roentgen.nuclides import Nuclide
    >>> import astropy.units as u
    >>> fe55 = Nuclide('fe', 55)
    >>> fe55.get_lines(1 * u.keV, 10 * u.keV)
    <QTable length=3>
    energy intensity origin parent
     keV
    float64  float64   str7   str7
    ------- --------- ------ ------
    5.88772       8.4    Mn   Fe-55
    5.89881     16.48    Mn   Fe-55
    6.513        3.38    Mn   Fe-55
    """

    def __init__(self, element: str, mass_number: int, descriptor: str = ""):
        try:
            filename = get_lara_file(element.capitalize(), mass_number, descriptor)
        except KeyError:
            raise KeyError(
                f"No match for mass_number {mass_number} for {element}. Valid mass numbers are {get_nuclide_mass_numbers(element)}"
            )
        self._line_tables = read_lara_tables(filename)
        if len(self._line_tables) > 1:
            self.lines = vstack(self._line_tables)
            self.lines.sort("energy")
        elif len(self._line_tables) == 1:
            self.lines = self._line_tables[0]
        self.meta = read_lara_header(filename)
        self.name = self.meta["Nuclide"]
        self.element = self.meta["Element"]
        self.half_life = self.meta["Half-life (s)"].to("yr")
        self.mass_number = mass_number
        self.data_sheet = f"http://www.lnhb.fr/nuclides/{self.name}_tables.pdf"
        if len(self.lines) == 0 or len(self._line_tables) == 1:
            self.decay_chain = f"{self.name}->{self.meta['Daughter(s)']}"
        else:
            self.decay_chain = "->".join(
                [this_table["parent"][0] for this_table in self._line_tables]
            )

    def __str__(self):
        return f"{self._text_summary()}{self.lines.__repr__()}"

    def __repr__(self):
        return f"{object.__repr__(self)}\n{self}"

    def _text_summary(self):
        num_lines = len(self.lines)
        result = f"Nuclide: {self.name}, ({self.element}) half life={self.half_life.to('year')} - ({num_lines:,} lines)\n"
        result += f"Decay chain: {self.decay_chain}\n"
        return result

    def get_lines(
        self,
        energy_low: u.Quantity[u.keV],
        energy_high: u.Quantity[u.keV],
        min_intensity: float = 0.0,
    ) -> QTable:
        """Returns a list of all emission lines in the energy range."""
        bool_array = (self.lines["energy"] < energy_high) * (
            self.lines["energy"] > energy_low
        )
        if min_intensity > 0:
            bool_array *= self.lines["intensity"] > min_intensity
        return self.lines[bool_array]


def get_nuclide_mass_numbers(element: str) -> list:
    """Return all available nuclide mass numbers for a given element."""
    bool_array = nuclides_list["symbol"] == element
    if sum(bool_array) > 0:
        mass_numbers = nuclides_list["mass_number"][bool_array].data
        return mass_numbers
    else:
        raise ValueError(f"No nuclide data found for element {element}")


def get_lara_file(element: str, mass_number: int, descriptor: str = "") -> Path:
    """Return path to the specified lara data file.

    Parameters
    ----------
    element : str
        The element or symbol name for the nuclide
    mass_number : int
        The mass number of the nuclide

    Returns
    -------
    file_path : Path
    """
    these_nuclides = nuclides_list.loc["symbol", element]
    if isinstance(these_nuclides, QTable):
        filename = f"{these_nuclides['symbol'][0]}-{mass_number}{descriptor}.lara.txt"
    else:
        filename = f"{these_nuclides['symbol']}-{mass_number}{descriptor}.lara.txt"
    # check if file exists
    if any(nuclides_list["filename"] == filename):
        return _lara_directory / filename
    else:
        raise FileNotFoundError(f"{filename} not found in nuclides_list.")
    # if isinstance(these_nuclides, Table):  # multiple isotopes exist
    #    this_nuclide = these_nuclides.loc["mass_number", mass_number]
    #    if len(descriptor) > 0:
    # else:
    #    this_nuclide = these_nuclides
    # return _lara_directory / str(this_nuclide["filename"])


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
    if len(table_line_index) > 0:
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
                parent = file_path.name.split(".")[0]
            skip_num = 2  # first table has the header line
            if j > 0:
                skip_num = 1
            for this_line in lines[this_index + skip_num :]:
                tokens = this_line.split(";")
                if len(tokens) > 1:
                    # check to see if an energy range is provided
                    if tokens[0].count("-") == 1:
                        energy1, energy2 = tokens[0].split(" - ")
                        # find average energy
                        energy.append(0.5 * (float(energy1) + float(energy2)))
                    else:
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
            this_table.sort("energy")
            result.append(this_table)
    else:
        # generate an empty table
        this_table = QTable()
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
            if key.count("Daughter"):
                value = tokens[2].rstrip().lstrip()
            else:
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
                value = value * u.s
            elif key.count("Decay"):
                value = value / u.s
            elif key.count("activity"):
                value = value * u.Bq / u.g
            result.update({key: value})
    return result
