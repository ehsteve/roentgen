"""
Detector energy-resolution (Gaussian pulse-height) broardening.

Port of SSW's ``pulse_spread.pro``. IDL's ``gaussint(x)`` is the
standard normal CDF, i.e. `scipy.stats.norm.cdf`.

Unlike the IDL version, the FWHM-vs-energy relationship is passed as a
plain Python callable (``fwhm_func(energy_kev) -> fwhm_kev``) rather than an
IDL function-name string plus a seperate parameter vector, which is the the
more idiomatic choice in Python.
"""

import numpy as np
from scipy.stats import norm

from roentgen.util.edges import edge_products

__all__ = ["pulse_spread_matrix"]

# FWHM -> Gaussian sigma. IDL literally divides by 2.36 (not the more precise
# 2 * sqrt(2 * ln(2)) = 2.3548...), kept as-is for numerical equivalence.
# TODO switch to more accurate value onc verified
_FWHM_TO_SIGMA = 2.36


def pulse_spread_matrix(eout_edges, fwhm_func):
    """The Gaussian pulse-spread matrix: the probability that a photon
    depositing its energy in output bin ``i`` is measured in output bin
    ``j``, for every ``(i, j)``.

    Parameters
    ----------
    eout_edges : `numpy.ndarray`
        Output (PHA channel) bin edges, keV, 1D ``(n+1,)`` or 2D ``(2, n)``.
    fwhm_func : callable
        ``fwhm_func(energy_kev) -> fwhm_kev``, the detector's energy
        resolution (FWHM) as a function of energy, both in keV.

    Returns
    -------
    `numpy.ndarray`
        ``(n_output, n_output)`` pulse-spread matrix. Multiply into an
        energy-loss matrix of shape ``(n_output, n_photon)`` with ``@`` to
        apply the broadening.
    """
    eout_bins = edge_products(eout_edges)
    eout_center = eout_bins.mean
    eout_edges_2 = eout_bins.edges_2
    n_output = eout_bins.width.size

    sigma = np.asarray(fwhm_func(eout_center), dtype=float) / _FWHM_TO_SIGMA

    pulse_spread = np.zeros((n_output, n_output))
    for i in range(n_output):
        z_upper = (eout_edges_2[1, :] - eout_center[i]) / sigma[i]
        z_lower = (eout_edges_2[0, :] - eout_center[i]) / sigma[i]
        pulse_spread[:, i] = norm.cdf(z_upper) - norm.cdf(z_lower)
    return pulse_spread
