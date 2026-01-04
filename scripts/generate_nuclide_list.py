from astropy.table import Table
import astropy.units as u

from roentgen import elements
from roentgen.nuclides.nuclides import _lara_directory, read_lara_header

_lara_files = list(_lara_directory.glob("*.txt"))

nuclides_list = Table(data={"filename": [this_file.name for this_file in _lara_files]})
nuclides_list["symbol"] = [this_row["filename"].split("-")[0] for this_row in nuclides_list]
mass_number = []
metastable = []
for this_row in nuclides_list:
    this_mass_number = this_row["filename"].split("-")[1].split(".")[0]
    if not this_mass_number[-1].isdigit():
        mass_number.append(int(this_mass_number[:-1]))
        metastable.append(True)
    else:
        mass_number.append(int(this_mass_number))
        metastable.append(False)

nuclides_list["mass_number"] = mass_number
nuclides_list["metastable"] = metastable
nuclides_list["name"] = [
    elements.loc["symbol", this_row["symbol"]]["name"] for this_row in nuclides_list
]
nuclides_list.sort("mass_number")

half_life = []
decay_lists = []


def decay_list_to_str(decay_list):
    result = ""
    for this_daughter in decay_list:
        result += f"({this_daughter['decay_mode']};{this_daughter['element']};{this_daughter['mass_number']};{this_daughter['branching_ratio']})"
    return result


# add half life to table
for this_nuc in nuclides_list:
    header = read_lara_header(_lara_directory / this_nuc["filename"])
    half_life.append(header["Half-life (s)"].to("year"))
    decay_list = header["Daughter(s)"]
    if len(decay_list) > 0:
        decay_lists.append(decay_list_to_str(decay_list))
    else:
        decay_lists.append("")

nuclides_list["decay_chain"] = decay_lists
nuclides_list["half_life [year]"] = u.Quantity(half_life, u.yr)
nuclides_list = nuclides_list[
    "symbol", "name", "mass_number", "half_life [year]", "metastable", "filename", "decay_chain"
]

print(nuclides_list)

nuclides_list.write("nuclides_list.csv", delimiter=",", format="ascii", overwrite=True)
