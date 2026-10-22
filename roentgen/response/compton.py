"""
Compton scattering physics: the Klein-Nishina differential cross-section and
the single-Compton-scatter contribuion to a detector response matrix.

Port of SSW's ``klein.pro``, ``klein_elctrn.pro``,
``compton_max_electron_energy.pro``, and ``resp_calc_compton_scatter.pro``.
Works on plain `numpy` arrays of energy values in keV, matching the
the reference implementation closely to ease validation.

``klein.pro``'s ``tot``/``enew``/``electron_max`` keyword outputs aren't
ported: ``resp_calc_compton_scatter.pro`` never uses them, so only
the differential cross-section itself is needed here.

By default `klein` is the plain free-electron Klein-Nishina cross-section,
same as ``klein.pro``. Passing ``elements`` multiplies it by the incoherent
scattering function S(q,Z) of J.H. Hubbell, Wm.J. Veigele, E.A. Briggs, R.T. Brown,
D.T. Cromer and R.J. Howerton, J. Phys. Chem. Ref. Data 4, 471 (1975);
erratum 6, 615 (1977), which accounts for bound-electron effects (suppresses
low-momentum-transfer scattering, more significant at low energy / high Z).
The tabulated data (z=1-98, see :download:`scripts/extract_scattering_function.py
<../../scripts/extract_scattering_function.py>`) is independent of the
cross-section database used for the interaction probabilities.
"""

import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.io import ascii

import roentgen
from roentgen.util.edges import edge_products, interp_linear_extrap

__all__ = [
    "incoherent_scattering_function",
    "compton_max_electron_energy",
    "klein",
    "klein_electron",
    "compton_scatter_matrix",
]

_ELECTRON_REST_ENERGY = (const.m_e * const.c**2).to_value(u.keV)  # ~511 keV
_CLASSICAL_ELECTRON_RADIUS = (const.alpha**2 * const.a0).to_value(u.cm)  # ~2.818e-13 cm
_HC_KEV_ANGSTROM = (const.h * const.c).to_value(u.keV * u.angstrom)  # ~12.398

_SCATTERING_FUNCTION_FILE = roentgen._data_directory / "incoherent_scattering_function.csv"
_scattering_function_table = None


def incoherent_scattering_function(z, q):
    """
    The incoherent scattering function S(q,Z) of an atom.

    Log-log interpolated from Hubbell et al. (1975), tends to ``z`` at large
    ``q`` and goes as ``q**2`` below the tabulated range.

    Parameters
    ----------
    z : `int`
        Atomic number, 1-98.
    q : `numpy.ndarray` or float
        Momentum transfer variable sin(theta/2)/wavelength, in 1/Angstrom.

    Returns
    -------
    `numpy.ndarray`
        S(q,Z), dimensionless.
    """
    global _scattering_function_table
    if _scattering_function_table is None:
        _scattering_function_table = ascii.read(
            _SCATTERING_FUNCTION_FILE, format="csv", comment="#"
        )
    column = f"z{int(z):02d}"
    if column not in _scattering_function_table.colnames:
        raise ValueError(f"No incoherent scattering function data for z={z} (supported: 1-98).")
    table_q = np.asarray(_scattering_function_table["q"], dtype=float)
    table_s = np.asarray(_scattering_function_table[column], dtype=float)

    q = np.asarray(q, dtype=float)
    # interpolate in log-log, the clipped q below table_q[0] is then scaled back with S ~ q**2
    clipped = np.clip(q, table_q[0], table_q[-1])
    s = np.exp(np.interp(np.log(clipped), np.log(table_q), np.log(table_s)))
    return np.where(q < table_q[0], table_s[0] * (q / table_q[0]) ** 2, s)


def _scattering_factor(photon_energy, mu, elements):
    """S(q,Z) averaged over ``elements`` per electron, so it multiplies a per-electron
    Klein-Nishina cross-section."""
    q = photon_energy / _HC_KEV_ANGSTROM * np.sqrt((1.0 - mu) / 2.0)
    atoms = sum(weight * incoherent_scattering_function(z, q) for z, weight in elements.items())
    electrons = sum(weight * z for z, weight in elements.items())
    return atoms / electrons


def compton_max_electron_energy(photon_energy):
    """
    Maximum kinetic energy of a Compton-scattered electron.

    From finematics (conservation of energy and momentum) alone.

    Parameters
    ----------
    photon_energy : `numpy.ndarray` or float
        Incident photon energy in keV.

    Returns
    -------
    Maximum scattered-electron kinetic energy, in keV.
    """
    alpha = photon_energy / _ELECTRON_REST_ENERGY
    return photon_energy / (0.5 / alpha + 1.0)


def klein(photon_energy, mu, elements=None):
    """
    Klein-Nishina differential Compton scattering cross-section, in
    cm2 electron-1 sr-1.

    Free-electron formula, unless ``elements`` is given in which case it is
    multiplied by the incoherent scattering function S(q,Z) (see the module
    docstring).

    Parameters
    ----------
    photon_energy : float
        Incident photon energy in keV (scalar).
    mu : `numpy.ndarray`
        Cosine of the angle between the incident and scattered photons.
    elements : `dict`, optional
        ``{atomic number: relative number of atoms}`` of the scatterer
        (e.g. mass fraction / atomic mass), enables the S(q,Z) correction.
    """
    alpha = photon_energy / _ELECTRON_REST_ENERGY
    one_minus_mu = 1.0 - mu
    alpha_one_minus_mu = alpha * one_minus_mu
    one_plus_alpha_term = 1.0 + alpha_one_minus_mu
    one_plus_mu2 = 1.0 + mu**2
    result = (
        _CLASSICAL_ELECTRON_RADIUS**2
        * (1.0 / one_plus_alpha_term) ** 3
        * (one_plus_mu2 / 2.0)
        * (1.0 + alpha_one_minus_mu**2 / one_plus_mu2 / one_plus_alpha_term)
    )
    if elements is not None:
        result = result * _scattering_factor(photon_energy, mu, elements)
    return result


def klein_electron(photon_energy, electron_energy, elements=None):
    """
    Klein-Nishina differential Compton cross-section as a function of a the kinetic energy of the scattered electron.

    Parameters
    ----------
    photon_energy : float
        Incident photon energy in keV (scalar).
    electron_energy : `numpy.ndarray`
        Kinetic energy of the scattered electron, in keV.
    elements : `dict`, optional
        As for `klein`.

    Returns
    -------
    `numpy.ndarray`
        Same shape as ``electron_energy``; zero wherever kinematically invalid.
    """
    alpha = photon_energy / _ELECTRON_REST_ENERGY
    electron_energy = np.asarray(electron_energy, dtype=float)
    differential_cross_section = np.zeros_like(electron_energy)

    mu = np.full_like(electron_energy, 2.0)  # out-of-range cosine sentinal, recieved from klein()
    kinematically_possible = photon_energy > electron_energy
    mu[kinematically_possible] = (
        1.0
        - electron_energy[kinematically_possible]
        / (photon_energy - electron_energy[kinematically_possible])
        / alpha
    )

    valid = np.abs(mu) <= 1.0
    if np.any(valid):
        dsigma_dmu = 2.0 * np.pi * klein(photon_energy, mu[valid], elements)
        differential_cross_section[valid] = (
            dsigma_dmu * (1.0 + alpha * (1.0 - mu[valid])) ** 2 / alpha / photon_energy
        )
    return differential_cross_section


def _value_locate(edges_1d, values):
    """Index i such that edges_1d[i] <= value < edges_1d[i+1] (ascending edges_1d)."""
    return np.searchsorted(edges_1d, values, side="right") - 1


def compton_scatter_matrix(channel_edges, photon_edges, compton_probability, elements=None):
    """The single-Compton-scatterer contribution to a detector response
    matrix: integrates the Klein-Nishina cross-section through the output
    channel edges for each input photon energy.

    Parameters
    ----------
    channel_edges : `numpy.ndarray`
        Output (energy-loss) bin edges, keV.
    photon_edges : `numpy.ndarray`
        Input photon energy bin edges, keV.
    compton_probability : `numpy.ndarray`
        Probability of a Compton interaction in the detector, one per input
        photon energy bin.
    elements : `dict`, optional
        As for `klein`; enables the S(q,Z) correction to the continuum shape.

    Returns
    -------
    `numpy.ndarray`
        ``(n_output_bins, n_photon_bins)`` matrix of the Compton-scatter
        contribution to the response.
    """
    photon_bins = edge_products(photon_edges)
    photon_center = photon_bins.gmean

    channel_bins = edge_products(channel_edges)
    channel_center = channel_bins.gmean
    channel_width = channel_bins.width
    channel_edges_1 = channel_bins.edges_1
    channel_edges_2 = channel_bins.edges_2

    n_photon_bins = photon_center.size
    n_channel_bins = channel_width.size
    compton_matrix = np.zeros((n_channel_bins, n_photon_bins))

    max_electron_energy = compton_max_electron_energy(photon_center)
    limited_max_energy = np.minimum(max_electron_energy, channel_edges_2[1, n_channel_bins - 1])
    max_channel_index = np.minimum(
        _value_locate(channel_edges_1, max_electron_energy), n_channel_bins - 1
    )
    valid_photon_bins = np.flatnonzero(max_electron_energy > channel_edges_1[0])
    if valid_photon_bins.size == 0:
        return compton_matrix

    for i in range(valid_photon_bins[0], n_photon_bins):
        channel_center_i = channel_center.copy()
        channel_width_i = channel_width.copy()
        if max_channel_index[i] >= 0:
            channel_center_i[max_channel_index[i]] = np.sqrt(
                channel_edges_2[0, max_channel_index[i]] * limited_max_energy[i]
            )
            channel_width_i[max_channel_index[i]] = (
                max_electron_energy[i] - channel_edges_2[0, max_channel_index[i]]
            )

        fine_grid_edges = np.geomspace(0.1, max_electron_energy[i] * 1.02, 101)
        fine_bins = edge_products(fine_grid_edges)
        fine_width, fine_center = fine_bins.width, fine_bins.gmean

        differential_cross_section = klein_electron(photon_center[i], fine_center, elements)
        total_cross_section = np.sum(fine_width * differential_cross_section)

        n_valid_channels = max_channel_index[i] + 1
        channel_cross_section = np.zeros(n_channel_bins)
        channel_cross_section[:n_valid_channels] = interp_linear_extrap(
            channel_center_i[:n_valid_channels], fine_center, differential_cross_section
        )

        channel_cross_section = (
            (channel_cross_section / total_cross_section) * channel_width_i * compton_probability[i]
        )
        compton_matrix[:, i] = channel_cross_section

    return compton_matrix
