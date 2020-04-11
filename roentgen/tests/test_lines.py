import pytest

import astropy.units as u

import roentgen
from roentgen.lines import get_lines

# remove H and He
all_elements = list(roentgen.elements['symbol'])[2:]


@pytest.mark.parametrize("element_str", all_elements)
def test_line(element_str):
    """Check that all elements return at least one line"""
    line_table = get_lines(0 * u.keV, 100 * u.keV, element=element_str)
    assert len(line_table) > 0


def test_get_lines():
    """Check that lines are returned in a range"""
    line_table = get_lines(4 * u.keV, 6 * u.keV)
    assert len(line_table) > 0


def test_get_all_lines():
    """Check that all lines are returned if the range is large enough."""
    line_table = get_lines(0 * u.keV, 200 * u.keV)
    assert len(line_table) == len(roentgen.emission_lines)


def test_get_right_number_of_lines():
    """Check a few specific cases for known number of lines"""
    assert len(get_lines(0 * u.eV, 520 * u.eV)) == 15
    assert len(get_lines(74900 * u.eV, 81000 * u.eV)) == 11
    assert len(get_lines(6000 * u.eV, 8000 * u.eV, element='Fe')) == 3
    assert len(get_lines(8000 * u.eV, 8100 * u.eV, element='Cu')) == 2

