"""
NIST XCOM photon cross-sections.

Current NIST data (unlike the McMaster et al. 1969/1970 coefficients in
`roentgen.absorption.mcmaster`), valid over a much wider energy range
(1 keV to 100 GeV) and additionally separating out coherent
(Rayleigh) scattering, incoherent (Compton) scattering, photoelectric
absorption, and nuclear/electron-field pair production.

Source: SSW's ``$SSW/packages/xray/dbase/xcom.fits`` (see
:download:`scripts/extract_xcom.py <../../scripts/extract_xcom.py>`), the
same data ``xsec.pro``'s ``use_xcom`` path reads via its own ``xcom`` IDL object.
Covers z=1-92 except z=91, for which SSW's own bundled data has no table
either (see `MISSING_ELEMENTS`).

Per-element tables are log/log interpolated (same convention as
`~roentgen.absorption.material.MassAttenuationCoefficient`): out-of-range
energies raise `ValueError` naturally from the interpolator, using each
element's own tabulated range rather than one shared cutoff.
"""

import numpy as np
from scipy import interpolate

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable

import roentgen
from roentgen.util import get_atomic_number

__all__ = [
    "MISSING_ELEMENTS",
    "photoelectric_mass_attenuation",
    "coherent_mass_attenuation",
    "incoherent_mass_attenuation",
    "total_mass_attenuation",
]

_xcom_directory = roentgen._data_directory / "xcom"

MISSING_ELEMENTS = sorted(
    set(range(1, 95)) - {int(p.stem[1:]) for p in _xcom_directory.glob("z*.csv")}
)

_COMPONENT_COLUMNS = ("coherent", "incoherent", "photoelectric", "total_with_coherent")


class _XCOMTable:
    """
    Per-element XCOM cross-section components, log/log interpolated.
    """

    def __init__(self, z):
        path = _xcom_directory / f"z{z:02d}.csv"
        if not path.exists():
            raise ValueError(f"No XCOM data for z={z}. Missing elements (z): {MISSING_ELEMENTS}.")
        table = QTable(ascii.read(path, format="csv", comment="#"))
        energy = np.asarray(table["energy_kev"], dtype=float)

        # duplicate energies mark absorption edges nudge the low-energy
        # side down so the grid is strictly increasing, same convention as
        # MassAttenuationCoefficient._remove_double_vals_from_data
        uniq, count = np.unique(energy, return_counts=True)
        for dup in uniq[count > 1]:
            idx = np.nonzero(energy == dup)[0]
            energy[idx[0]] -= 1e-6

        log_energy = np.log10(energy)
        self._interp = {
            column: interpolate.interp1d(
                log_energy,
                np.log10(np.asarray(table[column], dtype=float)),
                bounds_error=True,
                assume_sorted=True,
            )
            for column in _COMPONENT_COLUMNS
        }

    def evaluate(self, column, energy_kev):
        return 10 ** self._interp[column](np.log10(energy_kev))


_CACHE = {}


def _get_table(z):
    table = _CACHE.get(z)
    if table is None:
        table = _CACHE[z] = _XCOMTable(z)
    return table


def _resolve_z(element):
    return get_atomic_number(element) if isinstance(element, str) else int(element)


@u.quantity_input(energy=u.keV)
def _evaluate(element, energy, column):
    table = _get_table(_resolve_z(element))
    scalar_input = energy.isscalar
    energy_kev = np.atleast_1d(energy.to_value(u.keV)).astype(float)
    result = table.evaluate(column, energy_kev)
    result = result[0] if scalar_input else result
    return u.Quantity(result, "cm2/g")


def photoelectric_mass_attenuation(element, energy):
    """
    The photoelectric mass attenuation coefficient for a element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 keV to 1e8 keV for this element (out-of-range
        raises `ValueError`).

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    return _evaluate(element, energy, "photoelectric")


def coherent_mass_attenuation(element, energy):
    """
    The coherent (Rayleigh) scattering mass attenuation coefficient for an element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 keV to 1e8 keV for this element (out-of-range
        raises `ValueError`).

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    return _evaluate(element, energy, "coherent")


def incoherent_mass_attenuation(element, energy):
    """
    The incoherant (Compton) scattering mass attenuation coefficient for an element.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 keV to 1e8 keV for this element (out-of-range
        raises `ValueError`).

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    return _evaluate(element, energy, "incoherent")


def total_mass_attenuation(element, energy):
    """
    The total mass attenuation coefficient (photoelectric + coherent +
    incoherent), matching ``xsec.pro``'s ``use_xcom``, ``'AB'`` type.

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy, valid from 1 keV to 1e8 keV for this element (out-of-range
        raises `ValueError`).

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient in cm2/g.
    """
    return _evaluate(element, energy, "total_with_coherent")
