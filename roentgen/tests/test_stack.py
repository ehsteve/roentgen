import numpy as np
import pytest

import astropy.units as u

import roentgen
from roentgen.absorption.material import Material, Stack

all_materials = list(roentgen.elements["symbol"]) + list(roentgen.compounds["symbol"])
energy_array = u.Quantity(np.arange(1, 100, 1), "keV")


@pytest.fixture(params=all_materials)
def material(request):
    return Material(request.param, 500 * u.micron)


def test_twomaterials_to_stack(material):
    # check that adding two materials returns a compound
    assert isinstance(material + Material("Si", 500 * u.micron), Stack)


def test_threematerials_to_stack(material):
    # check that adding three materials returns a compound
    assert isinstance(
        material + Material("Ge", 500 * u.micron) + Material("cdte", 100 * u.micron), Stack
    )


def test_add_material_to_stack():
    stack1 = Material("Ge", 500 * u.micron) + Material("Ge", 100 * u.micron)
    assert isinstance(stack1 + Material("Ge", 500 * u.micron), Stack)


def test_add_two_stacks():
    material_str = ["Ge", "cdte", "Si", "Al"]
    stack1 = Material(material_str[0], 500 * u.micron) + Material(material_str[1], 100 * u.micron)
    stack2 = Material(material_str[2], 500 * u.micron) + Material(material_str[2], 500 * u.micron)
    total_stack = stack1 + stack2
    all_materials = [this_mat.symbol for this_mat in total_stack.materials]
    assert isinstance(total_stack, Stack)
    # check that the stack has all the same number of materials
    assert len(material_str) == len(total_stack.materials)
    # check that the stack has all the same material names
    assert all_materials.sort() == material_str.sort()


def test_stack_calculations(material):
    comp = Material("Ge", 500 * u.micron) + Material("cdte", 100 * u.micron)
    # test that it returns the same number of elements as energy array
    assert len(comp.absorption(energy_array)) == len(energy_array)
    assert len(comp.transmission(energy_array)) == len(energy_array)


def test_bad_add_to_stack():
    with pytest.raises(TypeError):
        Material("Ge", 500 * u.micron) + Material("Si", 100 * u.micron) + ["foo", "bar"]


def test_mat_add_stack():
    stack1 = Material("Ge", 500 * u.micron) + Material("Si", 100 * u.micron)
    assert isinstance(Material("Ge", 500 * u.micron) + stack1, Stack)
