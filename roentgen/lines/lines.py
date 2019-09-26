
import roentgen
from roentgen.absorption import get_atomic_number, is_an_element
import pandas as pd
import os

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory


def load_lines_for_element(element):
    if is_an_element(element):
        datafile_path = os.path.join(
            _data_directory, "emission_energies.csv")
        data = pd.read_csv(datafile_path, index_col=1)
    return data.loc[element]




data = load_lines_for_element('Si')
print(data)