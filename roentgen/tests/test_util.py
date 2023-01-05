import numpy as np
import pytest

import astropy.units as u
from astropy.units import cds, imperial

import roentgen
import roentgen.util as util

# this enables the atm unit
cds.enable()

all_materials = list(roentgen.elements["symbol"]) + list(roentgen.compounds["symbol"])
not_real_materials = ["adamantium", "ice-nine", "kryptonite", "redstone", "unobtainium"]


@pytest.mark.parametrize("element", roentgen.elements["symbol"])
def test_is_an_element_symbol(element):
    assert util.is_an_element(element)


@pytest.mark.parametrize("element", roentgen.elements["name"])
def test_is_an_element_name(element):
    assert util.is_an_element(element)


@pytest.mark.parametrize("element", not_real_materials)
def test_is_not_an_element(element):
    assert util.is_an_element(element) is False


@pytest.mark.parametrize("name,symbol", zip(roentgen.elements["name"], roentgen.elements["symbol"]))
def test_get_element_symbol(name, symbol):
    assert util.get_element_symbol(name) == symbol


@pytest.mark.parametrize("symbol", roentgen.elements["symbol"])
def test_get_element_symbol_same_return(symbol):
    assert util.get_element_symbol(symbol) == symbol


@pytest.mark.parametrize("element", not_real_materials)
def test_get_element_symbol_error(element):
    with pytest.raises(ValueError):
        util.get_element_symbol(element)


@pytest.mark.parametrize("symbol,z", zip(roentgen.elements["symbol"], roentgen.elements["z"]))
def test_get_atomic_number_symbol(symbol, z):
    assert util.get_atomic_number(symbol) == z


@pytest.mark.parametrize("name,z", zip(roentgen.elements["name"], roentgen.elements["z"]))
def test_get_atomic_number_name(name, z):
    assert util.get_atomic_number(name) == z


@pytest.mark.parametrize("element", not_real_materials)
def test_get_atomic_number_error(element):
    with pytest.raises(ValueError):
        util.get_atomic_number(element)


@pytest.mark.parametrize("symbol", roentgen.compounds["symbol"])
def test_is_in_known_compounds_symbol(symbol):
    assert util.is_in_known_compounds(symbol)


@pytest.mark.parametrize("name", roentgen.compounds["name"])
def test_is_in_known_compounds_name(name):
    assert util.is_in_known_compounds(name)


@pytest.mark.parametrize("element", not_real_materials)
def test_is_in_known_compounds_not(element):
    assert util.is_in_known_compounds(element) is False


@pytest.mark.parametrize("symbol", roentgen.compounds["symbol"])
def test_get_compound_index_symbol(symbol):
    assert util.get_compound_index(symbol) == roentgen.compounds.loc_indices[symbol]


@pytest.mark.parametrize("name", roentgen.compounds["name"])
def test_get_compound_index_name(name):
    assert util.get_compound_index(name) == roentgen.compounds["name"].tolist().index(name)


@pytest.mark.parametrize("element", not_real_materials)
def test_get_compound_index_error(element):
    with pytest.raises(ValueError):
        util.get_compound_index(element)


@pytest.mark.parametrize("element", not_real_materials)
def test_get_symbol_error(element):
    with pytest.raises(ValueError):
        util.get_material_symbol(element)


@pytest.mark.parametrize("element", not_real_materials)
def test_get_material_name(element):
    with pytest.raises(ValueError):
        util.get_material_name(element)


@pytest.mark.parametrize("element", not_real_materials)
def test_get_density_error(element):
    with pytest.raises(ValueError):
        util.get_material_density(element)


@pytest.mark.parametrize(
    "symbol,density", zip(roentgen.elements["symbol"], roentgen.elements["density"])
)
def test_get_density_element(symbol, density):
    assert util.get_material_density(symbol) == density


@pytest.mark.parametrize(
    "compound,density", zip(roentgen.compounds["symbol"], roentgen.compounds["density"])
)
def test_get_density_compound(compound, density):
    assert util.get_material_density(compound) == density


@pytest.mark.parametrize(
    "material", list(roentgen.compounds["symbol"]) + list(roentgen.compounds["symbol"])
)
def test_get_density_quantity(material):
    assert isinstance(util.get_material_density(material), u.Quantity)


@pytest.mark.parametrize(
    "pressure,temperature,result",
    [
        (1 * cds.atm, 15 * u.Celsius, 1.225 * u.kg / u.m**3),
        (101325 * u.pascal, 15 * u.Celsius, 1.225 * u.kg / u.m**3),
        (1 * cds.atm, 288.15 * u.K, 1.225 * u.kg / u.m**3),
        (1 * cds.atm, 59 * imperial.deg_F, 1.225 * u.kg / u.m**3),
        (5 * cds.atm, 20 * u.Celsius, 6.020 * u.kg / u.m**3),
    ],
)
def test_density_ideal_gas(pressure, temperature, result):
    assert np.isclose(util.density_ideal_gas(pressure, temperature), result, rtol=1e-4)


@pytest.mark.parametrize(
    "pressure,temperature",
    [
        (5 * u.m, 15 * u.Celsius),
        (5 * u.s, 15 * u.Celsius),
        (1 * cds.atm, 5 * u.m),
        (1 * cds.atm, 5 * u.gauss),
    ],
)
def test_density_ideal_gas_bad_units(pressure, temperature):
    with pytest.raises(u.UnitsError):
        util.density_ideal_gas(pressure, temperature)
