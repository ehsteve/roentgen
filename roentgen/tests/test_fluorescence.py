import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u

from roentgen.absorption.cross_section_database import get_database
from roentgen.absorption.material import Material
from roentgen.response import fluorescence
from roentgen.lines.lines import k_fluorescence_yield
from roentgen.response.fluorescence import fluorescence_escape

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
    name = fluorescence.SUPPORTED_DETECTORS[
        frozenset(s.capitalize() for s in detector.list_symbols)
    ]
    ported, _ = fluorescence._DISPATCH[name](_ENERGY, 0.0)
    generic, _ = fluorescence._generic(_ENERGY, detector, 0.0)
    assert_allclose(generic.sum(axis=1), ported.sum(axis=1), rtol=0.1, atol=1e-4)


def test_generic_zero_below_k_edge():
    detector = Material("Au", 1 * u.mm, database="mcmaster")  # K edge 80.7 keV
    probability, line_energy = fluorescence_escape(np.array([50.0, 100.0]) * u.keV, detector)
    assert np.all(probability[0] == 0)
    assert np.all(probability[1] > 0)
    assert np.all(line_energy < 80.7 * u.keV)


def test_generic_other_database():
    pytest.importorskip("xraylib")
    detector = Material("Au", 1 * u.mm, database="xraylib")
    probability, _ = fluorescence_escape(np.array([100.0]) * u.keV, detector)
    assert 0 < probability.sum() < 1


def test_ported_constants_only_for_mcmaster():
    energy = np.array([30.0]) * u.keV
    _, mcmaster_lines = fluorescence_escape(energy, Material("Ge", 1 * u.cm, database="mcmaster"))
    _, xcom_lines = fluorescence_escape(energy, Material("Ge", 1 * u.cm, database="xcom"))
    assert len(mcmaster_lines) == 2  # fluorescence.pro's Ka/Kb
    assert len(xcom_lines) > 2  # generic model


@pytest.mark.parametrize("database", ["mcmaster", "xcom", "xraylib", "xraydb"])
@pytest.mark.parametrize("element", ["Cd", "I", "Au"])
def test_k_shell_data_agrees_between_databases(database, element):
    pytest.importorskip("xraylib")
    pytest.importorskip("xraydb")
    k_edge, k_yield, energies, intensity = get_database(database).k_shell(element)
    ref_edge, ref_yield, _, _ = get_database("xraylib").k_shell(element)
    assert_allclose(k_edge.to_value(u.keV), ref_edge.to_value(u.keV), rtol=1e-3)
    assert_allclose(k_yield, ref_yield, rtol=0.05)
    assert energies.size == intensity.size > 0


@pytest.mark.parametrize("database", ["xraylib", "xraydb"])
def test_detector_database_supplies_k_shell_data(database, monkeypatch):
    pytest.importorskip(database)
    detector = Material("Au", 1 * u.mm, database=database)
    calls = []
    original = detector._database.k_shell
    monkeypatch.setattr(
        detector._database, "k_shell", lambda element: calls.append(element) or original(element)
    )
    probability, _ = fluorescence_escape(np.array([100.0]) * u.keV, detector)
    assert calls == ["Au"]
    assert 0 < probability.sum() < 1


@pytest.mark.parametrize(
    "z, fitted",
    # "fitted values" column of Table III.IV, Bambynek et al., Rev. Mod. Phys. 44, 716 (1972)
    [
        (14, 0.047),
        (26, 0.347),
        (32, 0.540),
        (48, 0.840),
        (53, 0.882),
        (79, 0.964),
        (82, 0.968),
        (92, 0.976),
    ],
)
def test_k_fluorescence_yield_reproduces_bambynek_table(z, fitted):
    assert_allclose(k_fluorescence_yield(z), fitted, atol=6e-4)
