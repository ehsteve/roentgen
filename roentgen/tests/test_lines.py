import pytest

import astropy.units as u
from astropy.table import QTable

import roentgen
import roentgen.lines
from roentgen.lines import get_edges, get_lines

# remove H and He
all_elements = list(roentgen.elements["symbol"])[2:]


@pytest.mark.parametrize("element_str", all_elements)
def test_line(element_str):
    """Check that all elements return at least one line"""
    line_list = get_lines(0 * u.keV, 100 * u.keV, element=element_str)
    assert len(line_list) > 0


def test_get_lines():
    """Check that lines are returned in a range"""
    assert len(get_lines(4 * u.keV, 6 * u.keV)) > 0


def test_get_all_lines():
    """Check that all lines are returned if the range is large enough."""
    assert len(get_lines(0 * u.keV, 200 * u.keV)) == len(roentgen.lines.emission_lines)


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


@pytest.mark.parametrize("element_str", all_elements)
def test_get_edges(element_str):
    """Test that all element return one row at least."""
    edge_list = get_edges(element_str)
    assert isinstance(edge_list, QTable)
    assert len(edge_list) > 0


@pytest.mark.parametrize(
    "element_str,edge_index, energy",
    [
        ("H", 0, 13.60 * u.eV),
        ("Fe", 0, 7112.00 * u.eV),
        ("Cu", 2, 952.30 * u.eV),
        ("Au", 8, 2206.00 * u.eV),
    ],
)
def test_get_edges_values(element_str, edge_index, energy):
    """Check a few specific cases"""
    assert get_edges(element_str)[edge_index]["energy"] == energy
