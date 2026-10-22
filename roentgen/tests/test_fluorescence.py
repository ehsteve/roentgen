import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.material import Material
from roentgen.response import fluorescence
from roentgen.response.fluorescence import fluorescence_escape, k_fluorescence_yield

_ENERGY = np.array([20.0, 40.0, 60.0, 100.0])


@pytest.mark.parametrize(
    "z, expected",
    # the K yields hard coded in fluorescence.pro
    [(48, 0.84), (52, 0.875), (53, 0.882), (55, 0.90), (32, 0.54)],
)
def test_k_fluorescence_yield_matches_fluorescence_pro(z, expected):
    assert_allclose(k_fluorescence_yield(z), expected, rtol=1e-2)


@pytest.mark.parametrize(
    "material, density",
    [
        ("Ge", None),
        ({"Na": 22.99, "I": 126.9}, 3.67),
        ({"Cs": 132.9, "I": 126.9}, 4.51),
        ({"Cd": 112.4, "Te": 127.6}, 5.85),
    ],
)
@pytest.mark.parametrize("database", ["mcmaster", "xcom"])
def test_generic_matches_ported_constants(material, density, database):
    # same physics as the hand ported blocks (which always read McMaster whatever the database
    # is), so the total escape probability should agree up to the database difference
    kwargs = {} if density is None else {"density": density * u.g / u.cm**3}
    detector = Material(material, 1 * u.cm, database=database, **kwargs)
    ported, _ = fluorescence_escape(_ENERGY * u.keV, detector)
    generic, _ = fluorescence._generic(_ENERGY, detector, 0.0)
    assert_allclose(generic.sum(axis=1), ported.sum(axis=1), rtol=0.1, atol=1e-4)


def test_generic_zero_below_k_edge():
    detector = Material("Au", 1 * u.mm, database="mcmaster")  # K edge 80.7 keV
    probability, line_energy = fluorescence_escape(np.array([50.0, 100.0]) * u.keV, detector)
    assert np.all(probability[0] == 0)
    assert np.all(probability[1] > 0)
    assert np.all(line_energy < 80.7)


def test_generic_other_database():
    pytest.importorskip("xraylib")
    detector = Material("Au", 1 * u.mm, database="xraylib")
    probability, _ = fluorescence_escape(np.array([100.0]) * u.keV, detector)
    assert 0 < probability.sum() < 1
