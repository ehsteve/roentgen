import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.material import Material, Response
from roentgen.response.compton import incoherent_scattering_function, klein

ENERGY = 30.0
MU = np.cos(np.array([0.3, 1.0, 2.5]))


def test_limits():
    assert_allclose(incoherent_scattering_function(26, 1e4), 26, rtol=1e-6)
    # S goes as q**2 below the table
    s1, s2 = incoherent_scattering_function(26, np.array([1e-4, 2e-4]))
    assert_allclose(s2 / s1, 4)


def test_unsupported_z():
    with pytest.raises(ValueError, match="z=99"):
        incoherent_scattering_function(99, 1.0)


def test_klein_default_unchanged_and_hydrogen_nearly_free():
    free = klein(ENERGY, MU)
    assert_allclose(klein(ENERGY, MU, {1: 1.0}), free, rtol=5e-2)
    assert np.all(klein(ENERGY, MU, {26: 1.0}) < free)


def test_matches_xraylib():
    xraylib = pytest.importorskip("xraylib")
    angle = np.arccos(MU)
    for z in (6, 26, 79):
        expected = [xraylib.SF_Compt(z, ENERGY / 12.3984 * np.sin(th / 2)) / z for th in angle]
        assert_allclose(klein(ENERGY, MU, {z: 1.0}) / klein(ENERGY, MU), expected, rtol=2e-3)


def test_response_matrix_scattering_function_changes_continuum_shape_only():
    detector = Material("Ge", 0.5 * u.cm, database="mcmaster")
    response = Response(Material("N", 1e-30 * u.um, database="mcmaster"), detector)
    edges = np.geomspace(10, 100, 60) * u.keV

    def matrix(**kwargs):
        return response.response_matrix(edges, **kwargs).eloss_matrix

    plain = matrix()
    corrected = matrix(scattering_function=True)
    assert not np.allclose(plain, corrected)
    # same total interaction probability per input energy, only redistributed (away from the top
    # of the output grid, where the continuum is truncated)
    assert_allclose(plain.sum(axis=0)[:50], corrected.sum(axis=0)[:50], rtol=1e-2)
