from astropy.table import Table
import astropy.units as u

from roentgen import elements
from roentgen.nuclides.nuclides import _lara_directory, read_lara_header

_lara_files = list(_lara_directory.glob("*.txt"))

nuclides_list = Table(data={"filename": [this_file.name for this_file in _lara_files]})
nuclides_list["symbol"] = [
    this_row["filename"].split("-")[0] for this_row in nuclides_list
]
mass_number = []
descriptor = []
for this_row in nuclides_list:
    this_mass_number = this_row["filename"].split("-")[1].split(".")[0]
    if not this_mass_number[-1].isdigit():
        mass_number.append(int(this_mass_number[:-1]))
        descriptor.append(this_mass_number[-1])
    else:
        mass_number.append(int(this_mass_number))
        descriptor.append(None)

nuclides_list["mass_number"] = mass_number
nuclides_list["descriptor"] = descriptor
nuclides_list["name"] = [
    elements.loc["symbol", this_row["symbol"]]["name"] for this_row in nuclides_list
]
nuclides_list.sort("mass_number")
print(nuclides_list)
half_life = []
# add half life to table
for this_nuc in nuclides_list:
    header = read_lara_header(_lara_directory / this_nuc["filename"])

    half_life.append(header["Half-life (s)"].to("year"))

nuclides_list["half_life [year]"] = u.Quantity(half_life, u.yr)
nuclides_list = nuclides_list[
    "symbol", "name", "mass_number", "half_life [year]", "descriptor", "filename"
]

nuclides_list.write("nuclides_list.csv", delimiter=",", format="ascii", overwrite=True)
