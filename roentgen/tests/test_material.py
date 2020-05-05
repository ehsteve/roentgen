import os
import pytest

import numpy as np

import astropy.units as u
from astropy.table import Table
from astropy.io import ascii

import roentgen
from roentgen.absorption.material import (MassAttenuationCoefficient, Material,
                                          Compound)
from roentgen.util import is_an_element

all_materials = list(roentgen.elements['symbol']) + list(roentgen.compounds['symbol'])
energy_array = u.Quantity(np.arange(1, 100, 1), 'keV')


@pytest.fixture(params=all_materials)
def mass_atten(request):
    return MassAttenuationCoefficient(request.param)


def test_is_an_element_symbol():
    """Test that function is okay with all known elements symbols"""
    for el in roentgen.elements['symbol']:
        assert(is_an_element(el))


def test_is_an_element_caseinsensitive_symbol():
    """Test that searching for an element symbol is not case sensitive"""
    for el in roentgen.elements['symbol']:
        assert is_an_element(el.upper())
        assert is_an_element(el.lower())
        assert is_an_element(el.capitalize())


def test_is_an_element_name():
    """Test that function is okay with all known elements names"""
    for el in roentgen.elements['name']:
        assert(is_an_element(el))


def test_is_an_element_caseinsensitive_name():
    """Test that searching for an element name is not case sensitive"""
    for el in roentgen.elements['name']:
        assert is_an_element(el.upper())
        assert is_an_element(el.lower())
        assert is_an_element(el.capitalize())


def test_mass_atten(mass_atten):
    # check that all materials can load a MassAttenuationCoefficient object
    assert isinstance(mass_atten, MassAttenuationCoefficient)


def test_returns_quantity(mass_atten):
    assert isinstance(mass_atten.func(1 * u.keV), u.Quantity)


def test_number_of_energies(mass_atten):
    energy = u.Quantity(np.arange(1, 1000), 'keV')
    atten = mass_atten.func(energy)
    assert len(energy) == len(atten)


@pytest.fixture(params=all_materials)
def material(request):
    return Material(request.param, 500 * u.micron)


def test_material(material):
    assert isinstance(material, Material)


def test_twomaterials_to_compound(material):
    # check that adding two materials returns a compound
    assert isinstance(material + Material('Si', 500 * u.micron), Compound)


def test_threematerials_to_compound(material):
    # check that adding three materials returns a compound
    assert isinstance(material + Material('Ge', 500 * u.micron) + Material('cdte', 100 * u.micron), Compound)


def test_compound_calculations(material):
    comp = Material('Ge', 500 * u.micron) + Material('cdte', 100 * u.micron)
    # test that it returns the same number of elements as energy array
    assert len(comp.absorption(energy_array)) == len(energy_array)
    assert len(comp.transmission(energy_array)) == len(energy_array)


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
