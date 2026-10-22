import json
from pathlib import Path

import numpy as np
import pytest

import astropy.units as u

from roentgen.absorption.material import Material, Response
from roentgen.response import compton
from roentgen.util import get_element_symbol, get_material_density

_SSW_TEST_DIR = Path(__file__).parent / "data" / "resp_calc_sswidl"
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


def _edges(elo, ehi, n):
    return np.geomspace(elo, ehi, n + 1) * u.keV


def _load_case(name):
    with open(_SSW_TEST_DIR / f"{name}.json") as f:
        params = json.load(f)
    arrays = np.load(_SSW_TEST_DIR / f"{name}.npz")
    return params, arrays


@pytest.mark.parametrize("case_name", _SSW_CASES)
def test_response_matrix_matches_idl(case_name, monkeypatch):
    """response_matrix() output matches real resp_calc.pro output (generated
    via hissw + local sswidl), for a representative set of detectors/windows/
    angles. See the ``.json`` files under roentgen/tests/data/resp_calc_sswidl/
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
    )
    if not params["nosigma"]:
        p0, p1 = params["func_par"]
        kwargs["fwhm"] = lambda e, p0=p0, p1=p1: p0 + p1 * np.sqrt(e)

    rm = response.response_matrix(_edges(params["elo"], params["ehi"], params["nflux"]), **kwargs)

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
        response.response_matrix(_edges(5, 50, 10))


def test_response_matrix_generic_fluorescence_material():
    # not in the fixed fluorescence list, so uses the generic model
    detector = Material("Au", 100 * u.micron, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)
    rm = response.response_matrix(_edges(5, 150, 10))
    assert np.all(np.isfinite(rm.eloss_matrix))
    assert np.all(rm.eloss_matrix >= 0)


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
    rm = response.response_matrix(_edges(10, 80, nflux))

    assert rm.eloss_matrix.shape == (nflux, nflux)
    assert np.all(np.isfinite(rm.eloss_matrix))
    assert np.all(rm.eloss_matrix >= 0)
    assert np.all(rm.efficiency >= 0) and np.all(rm.efficiency <= 1)


def test_response_matrix_with_broadening_normalises_reasonably():
    detector = Material("Si", 0.05 * u.cm, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)

    fwhm = lambda e: 0.15 + 0.0 * e  # noqa: E731

    rm = response.response_matrix(_edges(5, 60, 12), fwhm=fwhm)

    assert rm.pulse_height_matrix.shape == (12, 12)
    assert rm.smatrix.shape == (12, 12)
    assert np.all(np.isfinite(rm.pulse_height_matrix))
    assert np.all(np.isfinite(rm.smatrix))
    # broadening should conserve (approximately) the total counts per input bin
    np.testing.assert_allclose(
        rm.pulse_height_matrix.sum(axis=0), rm.eloss_matrix.sum(axis=0), atol=1e-3, rtol=1e-2
    )


def test_response_matrix_edges_out_differs_from_edges_in():
    detector = Material("Ge", 0.5 * u.cm, database="mcmaster")
    optical_path = Material("N", 1e-30 * u.um, database="mcmaster")
    response = Response(optical_path, detector)

    edges_in = _edges(10, 100, 20)
    edges_out = np.linspace(10, 100, 91) * u.keV
    rm = response.response_matrix(edges_in, edges_out)

    assert rm.eloss_matrix.shape == (90, 20)
    # the photopeak of an input bin lands in the output bin holding its centre energy
    centre_in = np.sqrt(edges_in[10].value * edges_in[11].value)
    centre_out = 0.5 * (edges_out[:-1] + edges_out[1:]).value
    assert abs(centre_out[np.argmax(rm.eloss_matrix[:, 10])] - centre_in) <= 1.0


@pytest.mark.parametrize("bad", [np.geomspace(100, 10, 11), np.array([10, 20, 20, 40, 80])])
def test_response_matrix_rejects_non_increasing_edges(bad):
    detector = Material("Si", 500 * u.micron, database="mcmaster")
    response = Response(Material("N", 1e-30 * u.um, database="mcmaster"), detector)
    with pytest.raises(ValueError, match="strictly increasing"):
        response.response_matrix(bad * u.keV)
    with pytest.raises(ValueError, match="strictly increasing"):
        response.response_matrix(_edges(10, 80, 8), bad * u.keV)
