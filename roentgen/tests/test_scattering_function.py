import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.material import Material, Response
from roentgen.response.compton import compton_scatter_matrix, incoherent_scattering_function, klein

ENERGY = 30.0
MU = np.cos(np.array([0.3, 1.0, 2.5]))


def test_limits():
    assert_allclose(incoherent_scattering_function(26, 1e4 / u.angstrom), 26, rtol=1e-6)
    # S goes as q**2 below the table
    s1, s2 = incoherent_scattering_function(26, np.array([1e-4, 2e-4]) / u.angstrom)
    assert_allclose(s2 / s1, 4)


def test_unsupported_z():
    with pytest.raises(ValueError, match="z=99"):
        incoherent_scattering_function(99, 1.0 / u.angstrom)


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


def _klein_nishina_total(energy):
    """Analytic total Klein-Nishina cross-section, cm2 per electron."""
    a = energy / 510.99895
    r2 = 7.940787e-26
    return (
        2
        * np.pi
        * r2
        * (
            (1 + a) / a**2 * (2 * (1 + a) / (1 + 2 * a) - np.log(1 + 2 * a) / a)
            + np.log(1 + 2 * a) / (2 * a)
            - (1 + 3 * a) / (1 + 2 * a) ** 2
        )
    )


@pytest.mark.parametrize("energy", [30.0, 150.0, 500.0])
def test_exact_klein_integrates_to_total_cross_section(energy):
    mu = np.linspace(-1, 1, 20001)
    exact = 2 * np.pi * np.trapezoid(klein(energy, mu, exact=True), mu)
    assert_allclose(exact, _klein_nishina_total(energy), rtol=1e-4)
    # klein.pro's version is low by a factor E'/E
    assert 2 * np.pi * np.trapezoid(klein(energy, mu, exact=False), mu) < 0.97 * exact


def test_response_matrix_ssw_compatible_changes_continuum():
    detector = Material("Ge", 0.5 * u.cm, database="mcmaster")
    response = Response(Material("N", 1e-30 * u.um, database="mcmaster"), detector)
    edges = np.geomspace(10, 100, 60) * u.keV
    corrected = response.response_matrix(edges).eloss_matrix
    ssw = response.response_matrix(edges, ssw_compatible=True).eloss_matrix
    assert not np.allclose(corrected, ssw)


def test_compton_last_bin_width_not_oversized():
    # 300 keV photons have a maximum electron energy of 162 keV, well above the last channel edge
    channel_edges = np.linspace(1, 50, 99)
    photon_edges = np.array([299.0, 301.0])
    corrected = compton_scatter_matrix(channel_edges, photon_edges, np.array([1.0]))
    ssw = compton_scatter_matrix(channel_edges, photon_edges, np.array([1.0]), exact=False)
    # the part of the continuum inside the channels can't exceed the Compton probability
    assert corrected.sum() < 1.0
    assert ssw.sum() > 1.0
