"""
Utilities for working with contiguous channel/bin boundary arrays.

Port of SSW's ``edge_products.pro``/``get_edges.pro``. Works on
plain `numpy` arrays of energy values in keV (no `astropy.units.Quantity`)
since it's an internal building block for
`~roentgen.response.response_matrix._compute_response_matrix` and
`~roentgen.response.compton.compton_scatter_matrix`, matching the
reference implementation closely to ease validation.

The ``/contiguous`` degenerate-edge-merging option present in the IDL
version isn't ported: no call site in ``resp_calc.pro`` uses it.
"""

from typing import NamedTuple

import numpy as np

__all__ = ["EdgeProducts", "edge_products", "interp_linear_extrap"]


class EdgeProducts(NamedTuple):
    """The commonly used derived quantities from a set of channel boundaries.

    A `NamedTuple` rather than a dataclass: every call site either destructures
    a couple of fields positionally (``mean, _, width, _, _ = edge_products(x)``)
    or reaches for one or two by name (``ep.mean``) -- both work here, without
    forcing every caller through attribute access for what's fundamentally a
    plain bundle of arrays.
    """

    # arithmetic mean of each channel's edges
    mean: np.ndarray
    # geometric mean of each channel's edges
    gmean: np.ndarray
    # absolute width of each channel
    width: np.ndarray
    # 2xn array of [low, high] edges per channel
    edges_2: np.ndarray
    # n+1 array of the edges of n contiguous channels
    edges_1: np.ndarray


def edge_products(edges):
    """From a vector of contiguous channel boundaries, return the commonly
    used derived quantities for plotting, scaling, and rebinning.

    Parameters
    ----------
    edges : array-like
        If 1D, n+1 contiguous channel boundaries for n channels.
        If 2D, ``(2, n)`` array of ``[low(i), high(i)]``, assumed contiguous
        for the purpose of deriving ``edges_1``.

    Returns
    -------
    `EdgeProducts`
    """
    edges = np.asarray(edges, dtype=float)

    if edges.size == 1:
        return EdgeProducts(
            mean=edges, gmean=edges, width=np.zeros_like(edges), edges_2=edges, edges_1=edges
        )

    if edges.ndim == 2 and edges.shape[0] == 2:
        edges_2 = edges
        n_channels = edges.shape[1]
        edges_1 = np.concatenate([edges_2[0, :], edges_2[1, n_channels - 1 : n_channels]])
    else:
        edges_1 = edges
        edges_2 = np.vstack([edges_1[:-1], edges_1[1:]])

    mean = np.sum(edges_2, axis=0) / 2.0
    gmean = np.sqrt(np.clip(edges_2[0] * edges_2[1], 0.0, None))
    width = np.abs(edges_2[1] - edges_2[0])

    return EdgeProducts(mean=mean, gmean=gmean, width=width, edges_2=edges_2, edges_1=edges_1)


def interp_linear_extrap(x, xp, fp):
    """Linear interpolation with linear extrapolation beyond the data range
    (matches IDL's INTERPOL, unlike `numpy.interp` which clips to the
    boundary values).

    Parameters mirror `numpy.interp`: query points ``x``, and a function
    given as ascending sample points ``xp`` with values ``fp`` at those
    points. Points below ``xp[0]``/above ``xp[-1]`` are extrapolated using
    the slope of the first/last segment, rather than clamped to
    ``fp[0]``/``fp[-1]``.

    Parameters
    ----------
    x : array-like
        Query points.
    xp : `numpy.ndarray`
        1D array of ascending sample points.
    fp : `numpy.ndarray`
        Values at ``xp``, same length as ``xp``.

    Returns
    -------
    `numpy.ndarray`
        Interpolated (and, outside ``[xp[0], xp[-1]]``, extrapolated)
        values at ``x``.
    """
    x = np.asarray(x, dtype=float)
    y = np.interp(x, xp, fp)
    below = x < xp[0]
    if np.any(below):
        slope = (fp[1] - fp[0]) / (xp[1] - xp[0])
        y[below] = fp[0] + slope * (x[below] - xp[0])
    above = x > xp[-1]
    if np.any(above):
        slope = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
        y[above] = fp[-1] + slope * (x[above] - xp[-1])
    return y
