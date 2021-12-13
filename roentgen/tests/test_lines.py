import pytest

import astropy.units as u

import roentgen
from roentgen.lines import get_lines

# remove H and He
all_elements = list(roentgen.elements["symbol"])[2:]


@pytest.mark.parametrize("element_str", all_elements)
def test_line(element_str):
    """Check that all elements return at least one line"""
    assert len(get_lines(0 * u.keV, 100 * u.keV, element=element_str)) > 0


def test_get_lines():
    """Check that lines are returned in a range"""
    assert len(get_lines(4 * u.keV, 6 * u.keV)) > 0


def test_get_all_lines():
    """Check that all lines are returned if the range is large enough."""
    assert len(get_lines(0 * u.keV, 200 * u.keV)) == len(roentgen.emission_lines)


@pytest.mark.parametrize(
    "energy_low,energy_high,element,result",
    [
        (0 * u.eV, 520 * u.eV, None, 15),
        (74.9 * u.keV, 81.0 * u.keV, None, 11),
        (6 * u.keV, 8 * u.keV, "Fe", 3),
        (8 * u.keV, 8.1 * u.keV, "Cu", 2),
    ],
)
def test_get_right_number_of_lines(energy_low, energy_high, element, result):
    """Check a few specific cases for known number of lines"""
    assert len(get_lines(energy_low, energy_high, element=element)) == result
