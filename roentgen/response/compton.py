"""
Compton scattering physics: the Klein-Nishina differential cross-section and
the single-Compton-scatter contribuion to a detector response matrix.

Port of SSW's ``klein.pro``, ``klein_elctrn.pro``,
``compton_max_electron_energy.pro``, and ``resp_calc_compton_scatter.pro``.
The public `compton_scatter_matrix` takes `astropy.units.Quantity` energies;
the underscore-prefixed kernels below work on plain `numpy` arrays of energy
values in keV, matching the reference implementation closely to ease
validation, and are internal.

``klein.pro``'s ``tot``/``enew``/``electron_max`` keyword outputs aren't
ported: ``resp_calc_compton_scatter.pro`` never uses them, so only
the differential cross-section itself is needed here.

.. todo::
    `_klein` uses the plain free-electron Klein-Nishina cross-section, same
    as ``klein.pro``. The physically more accurate version multiplies this
    by the incoherent scattering function S(q,Z), which accounts for
    bound-electron effects (suppresses low-momentum-transfer scattering,
    more significant at low energy / high Z). Source:
    J.H. Hubbell, Wm.J. Veigele, E.A. Briggs, R.T. Brown, D.T. Cromer, and
    R.J. Howerton, "Atomic Form Factors, Incoherent Scattering Functions,
    and Photon Scattering Cross Sections", J. Phys. Chem. Ref. Data 4, 471
    (1975); erratum 6, 615 (1977). Tabulated S(q,Z) data (digitised from
    that paper) is available via the ``xraylib`` package's
    ``SF.dat``.

.. todo::
    `_klein` raises the ``1 / (1 + alpha (1 - mu))`` term to the third power,
    as ``klein.pro`` does. The standard Klein-Nishina number cross-section
    uses the second power which under-weights large-angle scatters and so the
    region near the Compton edge. The second-power form matches ``xraylib.DCS_KN``
     and integrates to ``xraylib.CS_KN``
"""

import numpy as np

import astropy.constants as const
import astropy.units as u

from roentgen.util.edges import edge_products, interp_linear_extrap

__all__ = ["compton_scatter_matrix"]

_ELECTRON_REST_ENERGY = (const.m_e * const.c**2).to_value(u.keV)  # ~511 keV
_CLASSICAL_ELECTRON_RADIUS = (const.alpha**2 * const.a0).to_value(u.cm)  # ~2.818e-13 cm


def _compton_max_electron_energy(photon_energy):
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


def _klein(photon_energy, mu):
    """
    Klein-Nishina differential Compton scattering cross-section, in
    cm2 electron-1 sr-1.

    Free-electron formula only -- does not include the incoherent
    scattering function S(q,Z) bound-electron correction. See the module
    docstring note for the reference and a practical data source.

    Parameters
    ----------
    photon_energy : float
        Incident photon energy in keV (scalar).
    mu : `numpy.ndarray`
        Cosine of the angle between the incident and scattered photons.
    """
    alpha = photon_energy / _ELECTRON_REST_ENERGY
    one_minus_mu = 1.0 - mu
    alpha_one_minus_mu = alpha * one_minus_mu
    one_plus_alpha_term = 1.0 + alpha_one_minus_mu
    one_plus_mu2 = 1.0 + mu**2
    return (
        _CLASSICAL_ELECTRON_RADIUS**2
        * (1.0 / one_plus_alpha_term) ** 3
        * (one_plus_mu2 / 2.0)
        * (1.0 + alpha_one_minus_mu**2 / one_plus_mu2 / one_plus_alpha_term)
    )


def _klein_electron(photon_energy, electron_energy):
    """
    Klein-Nishina differential Compton cross-section as a function of a the kinetic energy of the scattered electron.

    Parameters
    ----------
    photon_energy : float
        Incident photon energy in keV (scalar).
    electron_energy : `numpy.ndarray`
        Kinetic energy of the scattered electron, in keV.

    Returns
    -------
    `numpy.ndarray`
        Same shape as ``electron_energy``; zero wherever kinematically invalid.
    """
    alpha = photon_energy / _ELECTRON_REST_ENERGY
    electron_energy = np.asarray(electron_energy, dtype=float)
    differential_cross_section = np.zeros_like(electron_energy)

    mu = np.full_like(electron_energy, 2.0)  # out-of-range cosine sentinal, recieved from _klein()
    kinematically_possible = photon_energy > electron_energy
    mu[kinematically_possible] = (
        1.0
        - electron_energy[kinematically_possible]
        / (photon_energy - electron_energy[kinematically_possible])
        / alpha
    )

    valid = np.abs(mu) <= 1.0
    if np.any(valid):
        dsigma_dmu = 2.0 * np.pi * _klein(photon_energy, mu[valid])
        differential_cross_section[valid] = (
            dsigma_dmu * (1.0 + alpha * (1.0 - mu[valid])) ** 2 / alpha / photon_energy
        )
    return differential_cross_section


def _value_locate(edges_1d, values):
    """Index i such that edges_1d[i] <= value < edges_1d[i+1] (ascending edges_1d)."""
    return np.searchsorted(edges_1d, values, side="right") - 1


@u.quantity_input
def compton_scatter_matrix(
    channel_edges: u.Quantity[u.keV], photon_edges: u.Quantity[u.keV], compton_probability
):
    """The single-Compton-scatterer contribution to a detector response
    matrix: integrates the Klein-Nishina cross-section through the output
    channel edges for each input photon energy.

    Parameters
    ----------
    channel_edges : `astropy.units.Quantity`
        Output (energy-loss) bin edges.
    photon_edges : `astropy.units.Quantity`
        Input photon energy bin edges.
    compton_probability : `numpy.ndarray`
        Probability of a Compton interaction in the detector, one per input
        photon energy bin.

    Returns
    -------
    `numpy.ndarray`
        Dimensionless ``(n_output_bins, n_photon_bins)`` matrix of the
        Compton-scatter contribution to the response (a probability).
    """
    photon_edges = photon_edges.to_value(u.keV)
    channel_edges = channel_edges.to_value(u.keV)
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

    max_electron_energy = _compton_max_electron_energy(photon_center)
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

        differential_cross_section = _klein_electron(photon_center[i], fine_center)
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
