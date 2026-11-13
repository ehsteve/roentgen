import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.material import Material
from roentgen.absorption.xcom import (
    MISSING_ELEMENTS,
    coherent_mass_attenuation,
    incoherent_mass_attenuation,
    photoelectric_mass_attenuation,
    total_mass_attenuation,
)

# real resp_calc.pro/xcom_obj->xsec() output, generated via hissw + local
# sswidl (same validation approach as as the McMaster/response_matrix tessts):
#   o = obj_new('xcom')
#   pe[i] = o->xsec(energy, z[i], 'PE')
_ELEMENTS = ["H", "C", "Si", "Fe", "Ge", "Cd", "Te", "Au"]
_PE_5KEV_IDL = [
    2.9317839e-02,
    1.8656542e01,
    2.4379736e02,
    1.3737897e02,
    2.4408688e02,
    7.6394366e02,
    8.9671320e02,
    6.5796136e02,
]
_PE_30KEV_IDL = [
    6.1658873e-05,
    5.7057630e-02,
    1.1606640e00,
    7.7628818e00,
    1.3353038e01,
    3.6884659e01,
    7.0745950e00,
    2.6113611e01,
]


def test_photoelectric_matches_idl_5kev():
    result = [photoelectric_mass_attenuation(el, 5 * u.keV).to_value("cm2/g") for el in _ELEMENTS]
    assert_allclose(result, _PE_5KEV_IDL, rtol=1e-5)


def test_photoelectric_matches_idl_30kev():
    result = [photoelectric_mass_attenuation(el, 30 * u.keV).to_value("cm2/g") for el in _ELEMENTS]
    assert_allclose(result, _PE_30KEV_IDL, rtol=1e-5)


def test_missing_elements():
    # SSW's own bundled xcom.fits has no table for z=91 either (mat_vs_z
    # raises "No table for z = 91"); z=93,94 are simply outside the range
    # mat_vs_z's lookup table (QORD) covers at all (it only maps z=1-92).
    assert MISSING_ELEMENTS == [91, 93, 94]


def test_missing_element_raises():
    with pytest.raises(ValueError):
        photoelectric_mass_attenuation(91, 30 * u.keV)


def test_energy_out_of_range_raises():
    with pytest.raises(ValueError):
        photoelectric_mass_attenuation("Si", 1e-3 * u.keV)
    with pytest.raises(ValueError):
        photoelectric_mass_attenuation("Si", 1e10 * u.keV)


@pytest.mark.parametrize("z", [z for z in range(1, 93) if z not in MISSING_ELEMENTS])
def test_all_elements_finite_and_positive(z):
    energy = u.Quantity([1, 10, 100, 1000, 1e5], "keV")
    for func in (
        photoelectric_mass_attenuation,
        coherent_mass_attenuation,
        incoherent_mass_attenuation,
        total_mass_attenuation,
    ):
        result = func(z, energy)
        assert np.all(np.isfinite(result))
        assert np.all(result.value > 0)


def test_total_is_sum_of_components():
    energy = u.Quantity([5, 30, 100], "keV")
    total = total_mass_attenuation("Ge", energy)
    components = (
        photoelectric_mass_attenuation("Ge", energy)
        + coherent_mass_attenuation("Ge", energy)
        + incoherent_mass_attenuation("Ge", energy)
    )
    assert_allclose(total.value, components.value, rtol=1e-6)


def test_material_database_xcom():
    detector = Material("Ge", 500 * u.um, database="xcom")
    pe = detector.photoelectric_mass_attenuation_coefficient(30 * u.keV)
    assert_allclose(pe.to_value("cm2/g"), 13.353038, rtol=1e-5)

    transmission = detector.transmission(u.Quantity([10, 30, 60], "keV"))
    assert np.all((transmission >= 0) & (transmission <= 1))


def test_material_default_database_rejects_component_access():
    detector = Material("Ge", 500 * u.um)  # default database="nist"
    with pytest.raises(NotImplementedError):
        detector.photoelectric_mass_attenuation_coefficient(30 * u.keV)
