import numpy as np
import pytest

import astropy.units as u

import roentgen
from roentgen.absorption import Substance

def test_substance_from_string():
    si1 = Substance('si')
    assert si1.name == 'Silicon'
    assert np.allclose(si1.mass_attenuation_coefficient([4] * u.keV), 452.9 * (u.cm ** 2) / u.g)
    assert np.allclose(si1.density, 2.33 * u.g / u.cm ** 3)

    si2 = Substance('si', density=4.66 * u.g / u.cm ** 3)
    assert si2.density == 4.66 * u.g / u.cm ** 3


def test_subtance_from_dict():
    with pytest.raises(ValueError):
        kapton = Substance({'H': 1.0})

    e = np.linspace(10, 1000, 1001) * u.keV

    h_dict = Substance({'H': 1.0}, name='h_dict')
    h = Substance('H')
    assert h_dict.name == 'h_dict'
    assert np.allclose(h_dict.mass_attenuation_coefficient(e), h.mass_attenuation_coefficient(e))

    mylar_calc = Substance({'H': 0.041960, 'C': 0.625016, 'O': 0.333024}, name='mylar_calc',
                           density=1.379 * u.g / u.cm**3)
    mylar_lookup = Substance('mylar')
    assert mylar_calc.name == 'mylar_calc'
    assert np.allclose(mylar_calc.mass_attenuation_coefficient(e),
                       mylar_lookup.mass_attenuation_coefficient(e))
