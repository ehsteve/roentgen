import numpy as np
import pytest

import astropy.units as u

import roentgen
from roentgen.absorption import Material, Response

all_materials = list(roentgen.elements["symbol"]) + list(roentgen.compounds["symbol"])
energy_array = u.Quantity(np.arange(1, 100, 1), "keV")
detector = Material("Si", 500 * u.micron)
thin_material = Material("air", thickness=1e-30 * u.um)
thick_material = Material("gold", thickness=1 * u.Mm)


@pytest.mark.parametrize("wrong_type", [(None, [1, 2, 3], 1, "Si")])
def test_response_detector_is_none_error(wrong_type):
    # check that the response class raises an error is not Stack or Material
    with pytest.raises(TypeError):
        Response(thin_material, detector=wrong_type)


def test_thin_detector():
    # check if a super thin material means a very small response
    resp = Response(optical_path=thin_material, detector=thin_material)
    assert resp.response(u.Quantity([1], "keV")) < 0.01


def test_thick_optical_path():
    # check if a super thick optical path means a very small response
    resp = Response(optical_path=thick_material, detector=thin_material)
    assert resp.response(u.Quantity([1], "keV")) < 0.01


def test_optical_path_stack_input():
    thin_material_stack = Material("air", thickness=1e-30 * u.um)
    for i in range(10):
        thin_material_stack += thin_material
    assert isinstance(Response(thin_material_stack, detector=detector), Response)


def test_bad_optical_path():
    with pytest.raises(TypeError):
        Response(optical_path="Si", detector=thin_material)


def test_raise_outside_of_data_range():
    """Test that ValueError is raised is trying to get values outside of data range 1 keV to 20 MeV."""
    resp = Response(optical_path=thin_material, detector=thin_material)

    energy = u.Quantity(np.arange(0.1, 10, 0.1), "keV")
    # below 1 keV
    with pytest.raises(ValueError):
        resp.response(energy)

    # above 20 MeV
    energy = u.Quantity(np.arange(10, 23, 0.1), "MeV")
    with pytest.raises(ValueError):
        resp.response(energy)


def test_repr_str():
    resp = Response(optical_path=Material("air", thickness=1e-30 * u.um), detector=thin_material)
    assert isinstance(resp.__repr__(), str)
    assert isinstance(resp.__str__(), str)
