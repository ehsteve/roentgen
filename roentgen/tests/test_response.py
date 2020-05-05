import numpy as np
import pytest

import astropy.units as u

import roentgen
from roentgen.absorption import Material, Response

all_materials = list(roentgen.elements['symbol']) + list(roentgen.compounds['symbol'])
energy_array = u.Quantity(np.arange(1, 100, 1), 'keV')

thin_material = Material("air", thickness=1e-30 * u.um)


def test_response_detector_is_none():
    # check that the response class accepts None for detector
    assert isinstance(Response(thin_material, detector=None), Response)
    # detector is none should assume perfectly absorbing material
    # make sure that response for a thin material and detector is None
    resp = Response(optical_path=thin_material, detector=None)
    assert resp.response(u.Quantity([1], 'keV')) > 0.99


def test_thin_detector():
    # check if a super thin material means no response
    resp = Response(optical_path=thin_material, detector=thin_material)
    assert resp.response(u.Quantity([1], 'keV')) < 0.01


def test_optical_path_input_list():
    thin_material_list = [Material("air", thickness=1e-30 * u.um)]
    assert isinstance(Response(thin_material_list, detector=None), Response)
    for i in range(10):
        thin_material_list += [thin_material]
    assert isinstance(Response(thin_material_list, detector=None), Response)


def test_detector_input_error():
    with pytest.raises(TypeError) as e_info:
        Response(thin_material, detector=0)
