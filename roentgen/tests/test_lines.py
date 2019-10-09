import pytest

import numpy as np
import astropy.units as u

from roentgen.lines import get_lines, get_lines_for_element
import roentgen


# remove H and He
all_elements = list(roentgen.elements['symbol'])[2:]

transition_names = roentgen.emission_energies.colnames[2:]
count_lines = 0
for this_trans in transition_names:
    count_lines += np.sum(roentgen.emission_energies[this_trans] != 0)


@pytest.mark.parametrize("element_str", all_elements)
def test_line(element_str):
    """Check that all elements return at least one line"""
    line_table = get_lines_for_element(element_str)
    assert len(line_table) > 0
    assert len(line_table) < 10  # can't have more than the transitions in the list


def test_get_lines():
    """Check that lines are returned in a range"""
    line_table = get_lines(4 * u.keV, 4.5 * u.keV)
    assert len(line_table) > 0


def test_get_all_lines():
    """Check that all lines are returned if the range is large enough."""
    line_table = get_lines(0 * u.keV, 200 * u.keV)
    assert len(line_table) == count_lines
