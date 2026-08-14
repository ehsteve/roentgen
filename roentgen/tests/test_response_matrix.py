import json
from pathlib import Path

import numpy as np
import pytest

import astropy.units as u

from roentgen.absorption.material import Material, Response
from roentgen.response import compton
from roentgen.util import get_element_symbol, get_material_density

_SSW_TEST_DIR = Path(__file__).parent / "data" / "resp_calc_golden"
_SSW_CASES = sorted(p.stem for p in _SSW_TEST_DIR.glob("*.json"))

# density (g/cm3) for detector materials made fo more than one element,
# matching the constants in in fluorescence.pro/det_xsec.pro
_COMPOUND_DETECTORS = {
    "NaI": ({"Na": 22.99, "I": 126.90}, 3.67),
    "CsI": ({"Cs": 132.9, "I": 126.9}, 4.51),
    "CdTe": ({"Cd": 112.4, "Te": 127.6}, 5.85),
}


def _make_detector(name, thickness_cm):
    if name in _COMPOUND_DETECTORS:
        composition, density = _COMPOUND_DETECTORS[name]
        return Material(
            composition, thickness_cm * u.cm, density=density * u.g / u.cm**3, database="mcmaster"
        )
    return Material(name, thickness_cm * u.cm, database="mcmaster")


def _make_optical_path(z_list, gmcm_list):
    materials = []
    for z, gmcm in zip(z_list, gmcm_list):
        symbol = get_element_symbol(int(z))
        density = get_material_density(symbol)
        thickness = (gmcm * u.g / u.cm**2) / density
        materials.append(Material(symbol, thickness, density=density, database="mcmaster"))
    optical_path = materials[0]
    for m in materials[1:]:
        optical_path = optical_path + m
    return optical_path


def _load_case(name):
    with open(_SSW_TEST_DIR / f"{name}.json") as f:
        params = json.load(f)
    arrays = np.load(_SSW_TEST_DIR / f"{name}.npz")
    return params, arrays


@pytest.mark.parametrize("case_name", _SSW_CASES)
def test_response_matrix_matches_idl(case_name, monkeypatch):
    """response_matrix() output matches real resp_calc.pro output (generated
    via hissw + local sswidl), for a representative set of detectors/windows/
    angles. See the ``.json`` files under roentgen/tests/data/resp_calc_golden/
    for parameters.

    roentgen.response.compton uses precise astropy.constants values for the
    electron rest energy and classical electron radius, but klein.pro et al.
    hardcode rounded 3 signricant figures (511.0, 2.82e-13) so monkeypatched
    back here so this stays a true comparison against SSW."""
    monkeypatch.setattr(compton, "_ELECTRON_REST_ENERGY", 511.0)
    monkeypatch.setattr(compton, "_CLASSICAL_ELECTRON_RADIUS", 2.82e-13)

    params, golden = _load_case(case_name)

    detector = _make_detector(params["detector"], params["d"])
    optical_path = _make_optical_path(params["z"], params["gmcm"])
    response = Response(optical_path, detector)

    kwargs = dict(
        theta=params["theta"] * u.deg,
        area=params["area"] * u.cm**2,
        nosigma=params["nosigma"],
    )
    if not params["nosigma"]:
        p0, p1 = params["func_par"]
        kwargs["fwhm"] = lambda e, p0=p0, p1=p1: p0 + p1 * np.sqrt(e)

    rm = response.response_matrix(
        params["nflux"], params["elo"] * u.keV, params["ehi"] * u.keV, **kwargs
    )

    np.testing.assert_allclose(rm.eloss_matrix, golden["eloss_matrix"], atol=1e-5, rtol=1e-3)
    if not params["nosigma"]:
        np.testing.assert_allclose(
            rm.pulse_height_matrix, golden["pulse_height_matrix"], atol=1e-5, rtol=1e-3
        )
        np.testing.assert_allclose(rm.smatrix, golden["smatrix"], atol=1e-5, rtol=1e-3)


def test_response_matrix_requires_mcmaster_detector():
    detector = Material("Si", 500 * u.micron)  # default database="nist"
    optical_path = Material("air", 1e-30 * u.um)
    response = Response(optical_path, detector)
    with pytest.raises(NotImplementedError):
        response.response_matrix(10, 5 * u.keV, 50 * u.keV, nosigma=True)


def test_response_matrix_requires_fwhm_unless_nosigma():
    detector = Material("Si", 500 * u.micron, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)
    with pytest.raises(ValueError):
        response.response_matrix(10, 5 * u.keV, 50 * u.keV)


def test_response_matrix_unsupported_fluorescence_material():
    # a McMaster-capable element that isn't in the fixed fluorescence list
    detector = Material("Au", 100 * u.micron, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)
    with pytest.raises(ValueError):
        response.response_matrix(10, 5 * u.keV, 50 * u.keV, nosigma=True)


@pytest.mark.parametrize("detector_name", ["Si", "Ge", "NaI", "CsI", "CdTe"])
def test_response_matrix_sanity(detector_name):
    """Basic shape/bounds sanity checks for all suported detector
    materials including CdTe, whose fluorescence normalisation
    deliberately diverges from resp_calc.pro see
    roentgen/response/fluorescence.py."""
    detector = _make_detector(detector_name, 0.1)
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)

    nflux = 8
    rm = response.response_matrix(nflux, 10 * u.keV, 80 * u.keV, nosigma=True)

    assert rm.eloss_matrix.shape == (nflux, nflux)
    assert np.all(np.isfinite(rm.eloss_matrix))
    assert np.all(rm.eloss_matrix >= 0)
    assert np.all(rm.efficiency >= 0) and np.all(rm.efficiency <= 1)


def test_response_matrix_with_broadening_normalises_reasonably():
    detector = Material("Si", 0.05 * u.cm, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)

    fwhm = lambda e: 0.15 + 0.0 * e  # noqa: E731

    rm = response.response_matrix(12, 5 * u.keV, 60 * u.keV, fwhm=fwhm)

    assert rm.pulse_height_matrix.shape == (12, 12)
    assert rm.smatrix.shape == (12, 12)
    assert np.all(np.isfinite(rm.pulse_height_matrix))
    assert np.all(np.isfinite(rm.smatrix))
    # broadening should conserve (approximately) the total counts per input bin
    np.testing.assert_allclose(
        rm.pulse_height_matrix.sum(axis=0), rm.eloss_matrix.sum(axis=0), atol=1e-3, rtol=1e-2
    )
