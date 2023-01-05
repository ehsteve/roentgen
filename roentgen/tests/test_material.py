import numpy as np
import pytest

import astropy.units as u

import roentgen
from roentgen.absorption.material import MassAttenuationCoefficient, Material
from roentgen.util import get_material_density, is_an_element

all_materials = list(roentgen.elements["symbol"]) + list(roentgen.compounds["symbol"])
energy_array = u.Quantity(np.arange(1, 100, 1), "keV")


@pytest.fixture(params=all_materials)
def mass_atten(request):
    return MassAttenuationCoefficient(request.param)


def test_is_an_element_symbol():
    """Test that function is okay with all known elements symbols"""
    for el in roentgen.elements["symbol"]:
        assert is_an_element(el)


def test_is_an_element_caseinsensitive_symbol():
    """Test that searching for an element symbol is not case sensitive"""
    for el in roentgen.elements["symbol"]:
        assert is_an_element(el.upper())
        assert is_an_element(el.lower())
        assert is_an_element(el.capitalize())


def test_is_an_element_name():
    """Test that function is okay with all known elements names"""
    for el in roentgen.elements["name"]:
        assert is_an_element(el)


def test_is_an_element_caseinsensitive_name():
    """Test that searching for an element name is not case sensitive"""
    for el in roentgen.elements["name"]:
        assert is_an_element(el.upper())
        assert is_an_element(el.lower())
        assert is_an_element(el.capitalize())


def test_mass_atten(mass_atten):
    # check that all materials can load a MassAttenuationCoefficient object
    assert isinstance(mass_atten, MassAttenuationCoefficient)


def test_returns_quantity(mass_atten):
    assert isinstance(mass_atten.func(1 * u.keV), u.Quantity)


def test_number_of_energies(mass_atten):
    energy = u.Quantity(np.arange(1, 1000), "keV")
    atten = mass_atten.func(energy)
    assert len(energy) == len(atten)


@pytest.fixture(params=all_materials)
def material(request):
    return Material(request.param, 500 * u.micron)


def test_material(material):
    assert isinstance(material, Material)


def test_material_give_density():
    assert isinstance(
        Material(all_materials[40], 500 * u.micron, density=1 * u.kg / u.m**3), Material
    )


@pytest.mark.parametrize(
    "material_wrong",
    [(1), (["Si", "He"]), (4.5)],
)
def test_material_bad_input(material_wrong):
    with pytest.raises(TypeError):
        Material(material_wrong, 500 * u.micron)


def test_bad_add_to_material():
    with pytest.raises(TypeError):
        Material("Ge", 500 * u.micron) + ["foo", "bar"]


@pytest.fixture(params=all_materials)
def thick_material(request):
    return Material(request.param, 500 * u.Mm)


def test_opaque(thick_material):
    # check that extremely large amounts of material mean no transmission
    assert thick_material.transmission(1 * u.keV) < 1e-6


@pytest.fixture(params=all_materials)
def thin_material(request):
    return Material(request.param, 1e-6 * u.pm)


def test_transparent(thin_material):
    # check that extremely large amounts of material mean no transmission
    assert thin_material.transmission(1 * u.keV) > 0.90


@pytest.mark.parametrize(
    "material,energy,thickness",
    [
        ("Si", 10 * u.keV, 1 * u.mm),
        ("Au", 1 * u.keV, 1 * u.micron),
        ("Air (dry)", 1 * u.keV, 1 * u.m),
        ("blood", 5 * u.keV, 1 * u.mm),
    ],
)
def test_linear_attenuation_coefficient(material, energy, thickness):
    mat = Material(material, thickness)
    assert np.isclose(
        mat.linear_attenuation_coefficient(energy),
        mat.mass_attenuation_coefficient(energy) * get_material_density(material),
        rtol=1e-4,
    )


@pytest.mark.parametrize(
    "material",
    [
        ("Si"),
        ("Au"),
        ("Air (dry)"),
        ("blood"),
        ({"Fe": 0.98, "C": 0.02}),  # steel
        ({"Cu": 0.88, "Sn": 0.12}),  # bronze
        ({"water": 0.97, "Na": 0.015, "Cl": 0.015}),  # salt water
    ],
)
def test_material_number_of_energies(material):
    mat = Material(material, 1 * u.m)
    energy = u.Quantity(np.arange(1, 1000), "keV")
    assert len(mat.transmission(energy)) == len(energy)
    assert len(mat.absorption(energy)) == len(energy)
    assert len(mat.mass_attenuation_coefficient(energy)) == len(energy)
    assert len(mat.linear_attenuation_coefficient(energy)) == len(energy)


@pytest.mark.parametrize(
    "material",
    [
        ("Si"),
        ("Au"),
        ("Air (dry)"),
        ("blood"),
        ({"Fe": 0.98, "C": 0.02}),  # steel
        ({"Cu": 0.88, "Sn": 0.12}),  # bronze
        ({"water": 0.97, "Na": 0.015, "Cl": 0.015}),  # salt water
    ],
)
def test_material_scalar_energy(material):
    """If a scalar energy is given then the results should NOT be arrays"""
    mat = Material(material, 1 * u.m)
    energy = 1 * u.keV
    assert isinstance(mat.transmission(energy), float)
    assert isinstance(mat.absorption(energy), float)
    assert mat.mass_attenuation_coefficient(energy).isscalar
    assert mat.linear_attenuation_coefficient(energy).isscalar


@pytest.mark.parametrize(
    "material",
    [
        ("Si"),
        ("Au"),
        ("Air (dry)"),
        ("blood"),
        ({"Fe": 0.98, "C": 0.02}),  # steel
        ({"Cu": 0.88, "Sn": 0.12}),  # bronze
        ({"water": 0.97, "Na": 0.015, "Cl": 0.015}),  # salt water
    ],
)
def test_mass_atten_calculation_in_material(material):
    mat = Material(material, 1 * u.m)
    energy = u.Quantity(np.arange(1, 1000), "keV")
    result1 = mat.mass_attenuation_coefficient(energy)[10]
    result2 = np.average(
        u.Quantity(
            [this_atten.func(energy[10]) for this_atten in mat.mass_attenuation_coefficients]
        ),
        weights=mat.fractional_masses,
    )
    assert np.isclose(result1, result2)


@pytest.mark.parametrize(
    "material_dict",
    [
        ({"Fe": 0.98, "C": 0.02}),  # steel
        ({"Cu": 0.88, "Sn": 0.12}),  # bronze
        ({"water": 0.97, "Na": 0.015, "Cl": 0.015}),  # salt water
    ],
)
def test_dict_input(material_dict):
    assert isinstance(Material(material_dict, 5 * u.mm), Material)


def test_density_calculation():
    a = {"Fe": 0.98, "C": 0.02}
    this_mat = Material(a, 5 * u.m)
    calc_density = get_material_density("Fe") * 0.98 + get_material_density("C") * 0.02
    assert np.isclose(this_mat.density, calc_density)
    a = {"Fe": 0.5, "C": 0.5}
    this_mat = Material(a, 5 * u.m)
    calc_density = get_material_density("Fe") * 0.5 + get_material_density("C") * 0.5
    assert np.isclose(this_mat.density, calc_density)
    a = {"Fe": 1, "C": 1}
    this_mat = Material(a, 5 * u.m)
    calc_density = get_material_density("Fe") * 0.5 + get_material_density("C") * 0.5
    assert np.isclose(this_mat.density, calc_density)
