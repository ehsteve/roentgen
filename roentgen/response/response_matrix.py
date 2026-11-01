"""
Detector response-matrix calculation a port of SSW's ``resp_calc.pro``.

Builds the full energy-redistribution matrix for a detector: for a grid of
input photon energies, the probability of a photon being measured in each
output (PHA) energy bin, including the photopeak (photoelectric capture),
single Compton-scatter escape, K-shell fluorescence escape peaks, and
(optionally) Gaussian resolution broadening.

Two divergences from ``resp_calc.pro``, both intentional and documented at
the call sites that make them: window/filter attenuation reuses
`~roentgen.absorption.material.Material.transmission`/`Stack.transmission`
(whichever cross-section database the optical path `Material` was built
with) rather than re-deriving it it from the McMaster data the way
``xsec.pro``'s ``'AB'`` type does; and the CdTe fluorescence branch fixes a
latent indexing bug in ``fluorescence.pro`` (see
`roentgen.response.fluorescence`).

Note this divergence isn't neutralised by building the window `Material`
with ``database="xcom"``: although ``xsec.pro`` *can* compute its ``'AB'``
type from XCOM data (via its ``use_xcom`` keyword), ``resp_calc.pro`` never
sets that keyword -- its window attenuation (and everything else) is always
McMaster-based. Only ``database="mcmaster"`` on the window actually matches
``resp_calc.pro``'s real numbers; ``"xcom"`` is a third, independently
different source, not an alternate match.
"""

from dataclasses import dataclass

import numpy as np

import astropy.units as u
from astropy.units import Quantity

from roentgen.response.compton import compton_scatter_matrix
from roentgen.response.fluorescence import fluorescence_escape
from roentgen.response.pulse_spread import pulse_spread_matrix
from roentgen.util.edges import edge_products, interp_linear_extrap

__all__ = ["ResponseMatrix"]


@dataclass
class ResponseMatrix:
    """
    Responce Matrix and and related data
    """

    # input photon energy bin edges, shape (n_photon + 1,)
    ein: u.Quantity
    # output (PHA channel) energy bin edges, shape (n_output + 1,)
    eout: u.Quantity
    # (n_output, n_photon) energy-loss matrix, before resolution broadening, in ct cm2 / ph:
    # ``eloss_matrix @ fluence`` gives counts per output bin for the photon fluence
    # (ph / cm2) in each input bin
    eloss_matrix: u.Quantity
    # (n_output, n_photon) pulse-height matrix, after resolution broadening
    # same units (``None`` if no ``fwhm`` was given)
    pulse_height_matrix: u.Quantity | None
    # pulse_height_matrix per unit output energy, in ct cm2 / (ph keV), i.e. counts/keV from
    # ``smatrix @ fluence`` (``None`` if no ``fwhm`` was given)
    smatrix: u.Quantity | None
    # relative (dimensionless) photopeak efficiency per input photon energy bin
    efficiency: np.ndarray
    # (n_output, n_output) dimensionless Gaussian resolution pulse-spread matrix (``None`` if no ``fwhm`` was given)
    pulse_shape: np.ndarray | None


def _edges_kev(value):
    return np.atleast_1d(value.to_value(u.keV)).astype(float)


@u.quantity_input
def _compute_response_matrix(
    response,
    nflux: int,
    elo: Quantity[u.keV],
    ehi: Quantity[u.keV],
    fwhm=None,
    eout: Quantity[u.keV] | None = None,
    edges_in: Quantity[u.keV] | None = None,
    edges_out: Quantity[u.keV] | None = None,
    theta: Quantity[u.deg] = 0 * u.deg,
    area: Quantity[u.cm**2] = 1 * u.cm**2,
    nosigma: bool = False,
    linearpha: bool = False,
):
    """Compute the full detector energy-redistribution (response) matrix.

    Private direct port of ``resp_calc.pro``, keeping its full set of
    options. The public entry point is `Response.response_matrix`, which
    forwards to this.

    Port of SSW's ``resp_calc.pro`` for a grid of input photon
    energies, the probability of a photon beeing measured in each output
    (PHA) energy bin, including the photopeak, single Compton-scatter
    escape, K-shell fluorescence escape peaks, and (unless ``nosigma``)
    Gaussian resolution broadening.

    ``response.detector`` must be an `~roentgen.absorption.material.Material`
    built with ``database="mcmaster"`` or ``database="xcom"`` (a single
    element, or a dict of elements), the photoelectric/Compton split this
    needs isn't available from the default NIST-backed data.

    Fluorescence escape further requires the detector's elements to match
    one of `roentgen.response.fluorescence.SUPPORTED_DETECTORS`
    (Ge, Si, NaI, CsI, CdTe, or SXS).

    Parameters
    ----------
    response : `~roentgen.absorption.material.Response`
        Supplies the the detector (``response.detector``) and window/filter
        stack (``response.optical_path``).
    nflux : `int`
        Number of input photon energy bins (log-spaced from ``elo`` to
        ``ehi``), unless ``edges_in`` is given.
    elo, ehi : `astropy.units.Quantity`
        Lower/upper energy limits of the input photon grid.
    fwhm : callable, optional
        ``fwhm(energy) -> fwhm``, the detector's energy resolution as a
        function of energy, taking and returning `astropy.units.Quantity`.
        Required unless ``nosigma=True``.
    eout : `astropy.units.Quantity`, optional
        Output (PHA channel) bin edges. Defaults to the input photon edges
        (or a linear grid if ``linearpha``).
    edges_in, edges_out : `astropy.units.Quantity`, optional
        Explicit input/output bin edges, overiding ``nflux``/``eout``.
    theta : `astropy.units.Quantity`, optional
        Angle of incidence from the detector normal. Default 0 deg.
    area : `astropy.units.Quantity`, optional
        Detector geometric area. Default 1 cm2.
    nosigma : `bool`, optional
        If `True`, skip resolution broadening (``pulse_height_matrix`` and
        ``smatrix`` on the result will be `None`).
    linearpha : `bool`, optional
        If `True` and ``eout`` isn't given, use a a linear (rather than log)
        output energy grid.

    Returns
    -------
    `ResponseMatrix`
        ``eloss_matrix`` and ``pulse_height_matrix`` are in ct cm2 / ph and
        ``smatrix`` in ct cm2 / (ph keV), as in ``resp_calc.pro``. They include
        the detector ``area`` but not the input bin widths, so
        ``pulse_height_matrix @ (flux * np.diff(ein) * exposure)`` is the counts
        per output bin for a photon flux density ``flux`` in ph / (keV s cm2)
        at the detector and an ``exposure`` time (and ``smatrix @ ...`` the
        counts per output keV). ``efficiency`` and ``pulse_shape`` are
        dimensionless.

    Raises
    ------
    NotImplementedError
        If ``response.detector`` wasn't built with ``database="mcmaster"``.
    ValueError
        If ``response.detector``'s elements aren't one of the fixed set
        `roentgen.response.fluorescence.SUPPORTED_DETECTORS`, or if
        ``fwhm`` isn't given and ``nosigma`` isn't `True`.
    """
    detector = response.detector
    optical_path = response.optical_path
    cos_theta = np.cos(theta.to_value(u.rad))

    elo_kev = elo.to_value(u.keV)
    ehi_kev = ehi.to_value(u.keV)

    if edges_in is not None:
        ein_edges = _edges_kev(edges_in)
    else:
        ein_edges = np.exp(np.arange(nflux + 1) * np.log(ehi_kev / elo_kev) / nflux) * elo_kev

    ein_bins = edge_products(ein_edges)
    ein_center = ein_bins.mean
    nflux = ein_bins.width.size

    if edges_out is not None:
        eout_edges = _edges_kev(edges_out)
    elif eout is not None:
        eout_edges = _edges_kev(eout)
    elif linearpha:
        eout_edges = np.arange(nflux + 1) * (ehi_kev - elo_kev) / nflux + elo_kev
    else:
        eout_edges = ein_edges

    eout_bins = edge_products(eout_edges)
    eout_center, eout_width = eout_bins.mean, eout_bins.width
    n_output_bins = eout_width.size
    output_bin_index = np.arange(n_output_bins, dtype=float)

    eloss_matrix = np.zeros((n_output_bins, nflux))

    # Photopeak diagonal
    peak_position = np.maximum(interp_linear_extrap(ein_center, eout_center, output_bin_index), 0.0)
    peak_flat_index = peak_position + np.arange(nflux) * n_output_bins
    peak_index_lo = np.floor(peak_flat_index).astype(int)
    peak_index_hi = np.minimum(peak_index_lo + 1, nflux * n_output_bins - 1)
    peak_weight_lo = 1.0 - (peak_flat_index - peak_index_lo)
    peak_weight_hi = 1.0 - peak_weight_lo

    ein_center_quantity = u.Quantity(ein_center, u.keV)
    photoelectric_coefficient = detector.photoelectric_linear_attenuation_coefficient(
        ein_center_quantity
    ).to_value(1 / u.cm)
    compton_coefficient = detector.compton_linear_attenuation_coefficient(
        ein_center_quantity
    ).to_value(1 / u.cm)
    total_coefficient = photoelectric_coefficient + compton_coefficient
    detector_thickness_cm = detector.thickness.to_value(u.cm)
    interaction_probability = 1.0 - np.exp(
        -np.minimum(total_coefficient * detector_thickness_cm / cos_theta, 40.0)
    )
    photoelectric_probability = (
        interaction_probability * photoelectric_coefficient / total_coefficient
    )
    compton_probability = interaction_probability * compton_coefficient / total_coefficient

    # window/filter attenuation -- reuses Material/Stack.transmission() directly
    # (see module docstring for why this differs from xsec.pro's 'AB' type)
    window_transmission = optical_path.transmission(ein_center_quantity) ** (1.0 / cos_theta)

    fluorescence_probability, fluorescence_energy = fluorescence_escape(
        ein_center_quantity, detector, theta=theta
    )
    fluorescence_energy = fluorescence_energy.to_value(u.keV)
    n_fluorescence_lines = fluorescence_energy.size
    total_fluorescence_probability = fluorescence_probability.sum(axis=1)

    efficiency = (
        (1.0 - total_fluorescence_probability) * window_transmission * photoelectric_probability
    )

    peak_row_lo, peak_col_lo = peak_index_lo % n_output_bins, peak_index_lo // n_output_bins
    peak_row_hi, peak_col_hi = peak_index_hi % n_output_bins, peak_index_hi // n_output_bins
    eloss_matrix[peak_row_lo, peak_col_lo] = efficiency * peak_weight_lo
    eloss_matrix[peak_row_hi, peak_col_hi] = efficiency * peak_weight_hi

    # Fluorescene escape peaks
    for line in range(n_fluorescence_lines):
        escape_energy = ein_center - fluorescence_energy[line]
        escape_position = np.maximum(
            interp_linear_extrap(escape_energy, eout_center, output_bin_index), 0.0
        )
        escape_flat_index = escape_position + np.arange(nflux) * n_output_bins
        valid = np.flatnonzero(escape_flat_index <= nflux * n_output_bins - 1)
        if valid.size == 0:
            continue
        escape_index_lo = np.floor(escape_flat_index[valid]).astype(int)
        escape_index_hi = np.minimum(escape_index_lo + 1, nflux * n_output_bins - 1)
        escape_weight_lo = 1.0 - (escape_flat_index[valid] - escape_index_lo)
        escape_weight_hi = 1.0 - escape_weight_lo
        escape_contribution = (
            fluorescence_probability[valid, line]
            * efficiency[valid]
            / (1.0 - total_fluorescence_probability[valid])
        )
        escape_row_lo, escape_col_lo = (
            escape_index_lo % n_output_bins,
            escape_index_lo // n_output_bins,
        )
        escape_row_hi, escape_col_hi = (
            escape_index_hi % n_output_bins,
            escape_index_hi // n_output_bins,
        )
        eloss_matrix[escape_row_lo, escape_col_lo] += escape_weight_lo * escape_contribution
        eloss_matrix[escape_row_hi, escape_col_hi] += escape_weight_hi * escape_contribution

    # Compton-scatter continuum
    eloss_matrix = eloss_matrix + compton_scatter_matrix(
        u.Quantity(eout_edges, u.keV), u.Quantity(ein_edges, u.keV), compton_probability
    )

    # Resolution brodening
    # projected detector area, converting photons per cm2 to counts
    drm_scale = area.to(u.cm**2) * cos_theta * u.ct / u.photon
    pulse_height_matrix = smatrix = pulse_shape = None
    if not nosigma:
        if fwhm is None:
            raise ValueError("fwhm must be given unless nosigma=True")
        pulse_shape = pulse_spread_matrix(u.Quantity(eout_edges, u.keV), fwhm)
        pulse_height_matrix = (pulse_shape @ eloss_matrix) * drm_scale
        smatrix = pulse_height_matrix / u.Quantity(eout_width, u.keV)[:, None]

    eloss_matrix = eloss_matrix * drm_scale

    return ResponseMatrix(
        ein=u.Quantity(ein_edges, u.keV),
        eout=u.Quantity(eout_edges, u.keV),
        eloss_matrix=eloss_matrix,
        pulse_height_matrix=pulse_height_matrix,
        smatrix=smatrix,
        efficiency=efficiency,
        pulse_shape=pulse_shape,
    )
