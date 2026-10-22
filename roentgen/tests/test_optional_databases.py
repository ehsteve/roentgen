import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.cross_section_database import get_database
from roentgen.absorption.material import Material

pytest.importorskip("xraylib")
pytest.importorskip("xraydb")

_DATABASES = ["xraylib", "xraydb"]
_ENERGY = u.Quantity([5, 10, 30, 100], "keV")


@pytest.mark.parametrize("database", _DATABASES)
def test_components_sum_to_total(database):
    fe = Material("Fe", 1 * u.mm, database=database)
    total = fe.mass_attenuation_coefficient(_ENERGY)
    parts = (
        fe.photoelectric_mass_attenuation_coefficient(_ENERGY)
        + fe.compton_mass_attenuation_coefficient(_ENERGY)
        + get_database(database).rayleigh("Fe", _ENERGY)
    )
    assert_allclose(total.to_value("cm2/g"), parts.to_value("cm2/g"), rtol=1e-6)


def test_xraylib_and_xraydb_agree_with_each_other():
    # both wrap the Elam tables
    a, b = (get_database(d).total("Fe", _ENERGY) for d in _DATABASES)
    assert_allclose(a.to_value("cm2/g"), b.to_value("cm2/g"), rtol=1e-3)


@pytest.mark.parametrize("database", _DATABASES)
def test_scalar_energy_returns_scalar(database):
    result = Material("Si", 1 * u.mm, database=database).mass_attenuation_coefficient(30 * u.keV)
    assert result.isscalar
    assert 0 < result.to_value("cm2/g") < 1e4


@pytest.mark.parametrize("database", _DATABASES)
@pytest.mark.parametrize("energy", [0.01 * u.keV, 1e4 * u.keV])
def test_out_of_range_energy_raises(database, energy):
    with pytest.raises(ValueError):
        get_database(database).total("Fe", energy)
