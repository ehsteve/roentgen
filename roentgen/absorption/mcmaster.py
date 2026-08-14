"""
Photoelectric, coherent (Rayleigh), and incoherent (Compton) x-ray mass attenuation cross-sections.

This is a port of SSW's ``xsec.pro``/``xcross.pro`` (classic, non-XCOM path), backed by the McMaster et al. 1969/1970
coefficient table (see :download:`xray_coef.csv <../../roentgen/data/xray_coef.csv>`, extracted from
``$SSW/packages/xray/dbase/xray_coef.fits``). Unlike roentgen's default NIST-based data
(see `MassAttenuationCoefficient`), this table gives access to the individual interaction components, which is what a
detector response matrix (photopeak vs. Compton-scatter continuum) needs. For the same kind of component split from
current NIST data instead of this 1969/1970 table, see `roentgen.absorption.xcom`.

.. warning::
    Valid only from 1 to 1000 keV, and only for elements with McMaster
    coefficients (87 of z=1-94; missing z are listed in `MISSING_ELEMENTS`).
"""

import numpy as np

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable

import roentgen
from roentgen.util import get_atomic_number

__all__ = [
    "MISSING_ELEMENTS",
    "VALID_ENERGY_RANGE",
    "xcross",
    "photoelectric_mass_attenuation",
    "rayleigh_mass_attenuation",
    "compton_mass_attenuation",
    "total_mass_attenuation",
]

VALID_ENERGY_RANGE = u.Quantity([1.0, 1000.0], "keV")

_coef_table = QTable(
    ascii.read(
        roentgen._data_directory / "xray_coef.csv",
        format="csv",
        comment="#",
    )
)
_coef_table.add_index("z")

MISSING_ELEMENTS = sorted(set(range(1, 95)) - set(int(z) for z in _coef_table["z"]))


def _resolve_z(element):
    return get_atomic_number(element) if isinstance(element, str) else int(element)


def _check_element(z):
    if z not in _coef_table["z"]:
        raise ValueError(
            f"No McMaster cross-section coefficients for z={z}. "
            f"Missing elements (z): {MISSING_ELEMENTS}."
        )


@u.quantity_input(energy=u.keV)
def _energy_kev_array(energy):
    if np.any(energy < VALID_ENERGY_RANGE[0]) or np.any(energy > VALID_ENERGY_RANGE[1]):
        raise ValueError(
            f"Energy out of range. McMaster cross-sections are valid from "
            f"{VALID_ENERGY_RANGE[0]} to {VALID_ENERGY_RANGE[1]}."
        )
    scalar_input = energy.isscalar
    energy_kev = np.atleast_1d(energy.to_value(u.keV)).astype(float)
    return energy_kev, scalar_input


def xcross(energy_kev, coef):
    """
    Evaluate the McMaster log/log cubic polynomial cross-section fit.

    Parameters
    ----------
    energy_kev : `numpy.ndarray`
        Energy in keV (plain array, no units).
    coef : array-like
        4 polynomial coefficients, lowest power (constant term) first.

    Returns
    -------
    cross_section : `numpy.ndarray`
        Cross-section in barns/atom.
    """
    log_energy = np.log(energy_kev)
    # np.polyval expects coefficients highest power first
    return np.exp(np.polyval(np.asarray(coef, dtype=float)[::-1], log_energy))


def _finish(result, scalar_input):
    result = result[0] if scalar_input else result
    return u.Quantity(result, "cm2/g")


def photoelectric_mass_attenuation(element, energy):
    """
    The photeolectric mass attenuation coefficient for an element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 to 1000 keV.

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    z = _resolve_z(element)
    _check_element(z)
    energy_kev, scalar_input = _energy_kev_array(energy)
    row = _coef_table.loc[z]
    pe_range_edges = np.array([row[f"edge{i}"] for i in range(5)])
    result = np.zeros_like(energy_kev)
    for i in range(4):
        in_range = (energy_kev >= pe_range_edges[i]) & (energy_kev < pe_range_edges[i + 1])
        if np.any(in_range):
            pe_coef = [row[f"pe_{i}_{j}"] for j in range(4)]
            result[in_range] = xcross(energy_kev[in_range], pe_coef)
    result = result / row["conv"]
    return _finish(result, scalar_input)


def rayleigh_mass_attenuation(element, energy):
    """
    The coherent (Rayleigh) scattering mass attenuation coefficient for an element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 to 1000 keV.

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    z = _resolve_z(element)
    _check_element(z)
    energy_kev, scalar_input = _energy_kev_array(energy)
    row = _coef_table.loc[z]
    rayleigh_coef = [row[f"scc_{j}"] for j in range(4)]
    result = xcross(energy_kev, rayleigh_coef) / row["conv"]
    return _finish(result, scalar_input)


def compton_mass_attenuation(element, energy):
    """
    The incoherent (Compton) scattering mass attenuation coefficient for an element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 to 1000 keV.

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    z = _resolve_z(element)
    _check_element(z)
    energy_kev, scalar_input = _energy_kev_array(energy)
    row = _coef_table.loc[z]
    compton_coef = [row[f"sci_{j}"] for j in range(4)]
    result = xcross(energy_kev, compton_coef) / row["conv"]
    return _finish(result, scalar_input)


def total_mass_attenuation(element, energy):
    """
    The total mass attenuation coefficient (photoelectric + coherent + incoherent) for an element, matching
    ``xsec.pro``'s ``'AB'`` type.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 to 1000 keV.

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    return (
        photoelectric_mass_attenuation(element, energy)
        + rayleigh_mass_attenuation(element, energy)
        + compton_mass_attenuation(element, energy)
    )
