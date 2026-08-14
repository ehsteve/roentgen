"""
K-shell x-ray fluorescence escape from a detector.

Port of SSW's ``fluorescence.pro``. Unlike photoelectric/Compton
cross-sections (`roentgen.absorption.mcmaster`), K-shell fluorescence
yields and branching ratios aren't derivable from the McMaster coefficient
table they're empirical (Bambynek-type) values, and the reference
implementation only tabulates them for six specific detector materials. This
module ports those those constants and combination formulas as literal
values, same limitation as as upstream: only Ge, Si, NaI, CsI, CdTe, and SXS
(Xe/CO2) detectors are supported.

Each detector-specific block below computes, for one or more K-shell escape
lines, a "ratio" of the total photoelectric interaction probability at the
the incident energy to that at the escape lines own energy. Only *relative*
element weights matter for this ratio (any overall multiplicative scale
factor applied identically to the numerator and denominator cancels out),
which is why several combination formulas below don't bother normalising to
a mass fraction the way ``det_xsec.pro`` documents them — they're kept
exactly as SSW computes them, including one apparent copy/paste quirk in the
CsI case (see below), for numerical equivalence with ``resp_calc.pro``.

.. todo::
    Given ``roengen``'s additional line support `roentgen.lines`
    (emission_lines.csv has energy + relative intensity per transition, e.g.
    Kα1,2 vs Kβ) and K-edge energies already in electron_binding_energies.csv
    The one missing piece is the absolute K-shell fluorescence yield w_K(Z)
    (probability of X-ray emission vs. Auger electron) but could use
    semi-empirical formula (Bambynek/Krause-type fit) that covers Z≈5–100
    without needing a new bundled table the fluorescence code could be expanded
    to support more elements over a wider energy range.
"""

import numpy as np

import astropy.units as u

from roentgen.absorption.mcmaster import photoelectric_mass_attenuation

__all__ = ["SUPPORTED_DETECTORS", "fluorescence_escape"]

SUPPORTED_DETECTORS = {
    frozenset({"Ge"}): "Ge",
    frozenset({"Si"}): "Si",
    frozenset({"Na", "I"}): "NaI",
    frozenset({"Cs", "I"}): "CsI",
    frozenset({"Cd", "Te"}): "CdTe",
    frozenset({"Xe", "C", "O"}): "SXS",
}


def _photoelectric(element, energy_kev):
    """Photoelectric mass attenuation coefficient, plain `numpy.ndarray` in/out."""
    return photoelectric_mass_attenuation(element, u.Quantity(energy_kev, "keV")).to_value("cm2/g")


def _finalise(
    energy_kev,
    incident_weighted_pe,
    escape_weighted_pe,
    yield_weight,
    k_edge_energy,
    escape_line_energy,
    theta_deg,
):
    """Shared tail of fluorescence.pro: build escape_probability[num, n_lines]
    from the per-line total-interaction ratio and yield weight, zeroing lines
    below their K edge."""
    cos_theta = np.cos(np.radians(theta_deg))
    n_escape_lines = len(escape_line_energy)
    escape_probability = np.zeros((energy_kev.size, n_escape_lines))
    for line in range(n_escape_lines):
        ratio = (incident_weighted_pe / cos_theta) / escape_weighted_pe[line]
        escape_probability[:, line] = (
            yield_weight[:, line] / 2.0 * (1.0 - (1.0 / ratio) * np.log(1.0 + ratio))
        )
        escape_probability[energy_kev < k_edge_energy[line], line] = 0.0
    return escape_probability, np.asarray(escape_line_energy, dtype=float)


def _ge(energy_kev, theta_deg):
    incident_weighted_pe = _photoelectric("Ge", energy_kev)
    escape_line_energy = [9.87, 10.98]
    k_edge_energy = [11.1, 11.1]
    escape_weighted_pe = _photoelectric("Ge", escape_line_energy)
    yield_weight = np.tile([0.471, 0.0691], (energy_kev.size, 1))
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


def _si(energy_kev, theta_deg):
    incident_weighted_pe = _photoelectric("Si", energy_kev)
    escape_line_energy = [1.839 - 0.118]
    k_edge_energy = [1.839]
    escape_weighted_pe = _photoelectric("Si", escape_line_energy)
    yield_weight = np.tile([0.03], (energy_kev.size, 1))
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


def _nai(energy_kev, theta_deg):
    pe_na = _photoelectric("Na", energy_kev) * 23.0
    pe_i = _photoelectric("I", energy_kev) * 126.9
    incident_weighted_pe = pe_na + pe_i
    iodine_fraction = pe_i / incident_weighted_pe

    escape_line_energy = [29.0]
    k_edge_energy = [33.2]
    escape_weighted_pe = [
        _photoelectric("Na", escape_line_energy)[0] * 23.0
        + _photoelectric("I", escape_line_energy)[0] * 126.9
    ]

    yield_weight = (0.882 * iodine_fraction)[:, None]
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


def _csi(energy_kev, theta_deg):
    pe_cs = _photoelectric("Cs", energy_kev) * 132.9
    pe_i = _photoelectric("I", energy_kev) * 126.9
    incident_weighted_pe = pe_cs + pe_i
    cesium_fraction = pe_cs / incident_weighted_pe
    iodine_fraction = pe_i / incident_weighted_pe

    escape_line_energy = [30.2, 29.0]
    k_edge_energy = [36.0, 33.2]
    # NOTE: fluorescence.pro's CsI escape_weighted_pe uses the Cs weight
    # (132.9) for *both* terms (`xsec(ke,55.,...)*132.9 + xsec(ke,53.,...)*132.9`),
    # not 126.9 for the iodine term. Kept as-is for numerical equivalence.
    escape_weighted_pe = [
        132.9 * (_photoelectric("Cs", [e])[0] + _photoelectric("I", [e])[0])
        for e in escape_line_energy
    ]

    yield_weight = np.stack([0.90 * cesium_fraction, 0.882 * iodine_fraction], axis=1)
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


def _cdte(energy_kev, theta_deg):
    # NOTE: fluorescence.pro's CdTe branch broadcasts mu0/mu1 to (num, 8)
    # arrays and then indexes them in the shared kp loop with a bare scalar
    # `mu0[i]`.
    #
    # So for num > 1 input
    # energies (i.e. any real resp_calc.pro call, which always has more than
    # one photon energy bin) most of the 8 escape lines get normalised
    # against the wrong K-line's mu0. This is a latent bug in the reference
    # implementation. Thi s port implements the clearly-intended
    # per-line indexing instead , so CdTe fluorescence here will disagree with real
    # resp_calc.pro output whenever more than one photon energy is used.

    # local branching fractions use one weight pair (186.6, 211.9)...
    pe_cd_local = _photoelectric("Cd", energy_kev) * 186.6
    pe_te_local = _photoelectric("Te", energy_kev) * 211.9
    local_total = pe_cd_local + pe_te_local
    cadmium_fraction = pe_cd_local / local_total
    tellurium_fraction = pe_te_local / local_total

    cd_line_weights = np.array([53.2, 100.0, 27.3, 5.3])
    cd_line_weights = 0.84 * cd_line_weights / cd_line_weights.sum()
    te_line_weights = np.array([53.7, 100.0, 28.2, 5.8])
    te_line_weights = 0.875 * te_line_weights / te_line_weights.sum()

    escape_line_energy = [22.9841, 23.1736, 26.09, 26.684, 27.2, 27.47, 31.0, 31.7]
    k_edge_energy = [26.711] * 4 + [31.813] * 4

    # ...while the total-interaction ratio (incident/escape weighted pe) uses
    # det_xsec.pro's own CdTe combination, a *different* weight pair
    # (112.4, 127.6). Both are ported verbatim from fluorescence.pro/det_xsec.pro.
    def cdte_mixture_pe(x):
        return _photoelectric("Cd", x) * 112.4 + _photoelectric("Te", x) * 127.6

    incident_weighted_pe = cdte_mixture_pe(energy_kev)
    escape_weighted_pe = cdte_mixture_pe(escape_line_energy)

    yield_weight = np.tile(np.concatenate([cd_line_weights, te_line_weights]), (energy_kev.size, 1))
    yield_weight[:, 0:4] *= cadmium_fraction[:, None]
    yield_weight[:, 4:8] *= tellurium_fraction[:, None]
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


def _sxs(energy_kev, theta_deg):
    pe_xe = _photoelectric("Xe", energy_kev) * 0.9 * 131.3
    pe_c = _photoelectric("C", energy_kev) * 0.1 * 12.0
    pe_o = _photoelectric("O", energy_kev) * 0.1 * 32.0
    incident_weighted_pe = pe_xe + pe_c + pe_o
    xenon_fraction = pe_xe / incident_weighted_pe

    escape_line_energy = [29.6]
    k_edge_energy = [34.58]
    escape_weighted_pe = [
        _photoelectric("Xe", escape_line_energy)[0] * 0.9 * 131.3
        + _photoelectric("C", escape_line_energy)[0] * 0.1 * 12.0
        + _photoelectric("O", escape_line_energy)[0] * 0.1 * 32.0
    ]

    yield_weight = (0.90 * xenon_fraction)[:, None]
    return _finalise(
        energy_kev,
        incident_weighted_pe,
        escape_weighted_pe,
        yield_weight,
        k_edge_energy,
        escape_line_energy,
        theta_deg,
    )


_DISPATCH = {
    "Ge": _ge,
    "Si": _si,
    "NaI": _nai,
    "CsI": _csi,
    "CdTe": _cdte,
    "SXS": _sxs,
}


def fluorescence_escape(energy, detector, theta=0 * u.deg):
    """
    K-shell fluorescence escape probabilities and energies for a detector.

    Parameters
    ----------
    energy : `astropy.units.Quantity`
        Incident photon energy.
    detector : `~roentgen.absorption.material.Material`
        The detector material. Must be a single element or a dict-of-elements
        `Material` matching one of `SUPPORTED_DETECTORS` (Ge, Si, NaI
        i.e. {Na, I}, CsI i.e. {Cs, I}, CdTe i.e. {Cd, Te}, or SXS i.e.
        {Xe, C, O}) — the same fixed set the reference implementation
        supports, since fluorescence yield data isn't tabulated elsewhere.
    theta : `astropy.units.Quantity`, optional
        Angle of incidence from the detector normal. Default 0 deg.

    Returns
    -------
    escape_probability : `numpy.ndarray`
        ``(n_energy, n_lines)`` probability of fluorescence escape at each
        line, per incident photon.
    escape_line_energy : `numpy.ndarray`
        Escape line energies, in keV.
    """
    symbols = frozenset(s.capitalize() for s in detector.list_symbols)
    name = SUPPORTED_DETECTORS.get(symbols)
    if name is None:
        supported = ", ".join(sorted(SUPPORTED_DETECTORS.values()))
        raise ValueError(
            f"No fluorescence escape data for detector material "
            f"{sorted(symbols)}. Supported: {supported}."
        )
    theta_deg = theta.to_value(u.deg)
    energy_kev = np.atleast_1d(energy.to_value(u.keV)).astype(float)
    return _DISPATCH[name](energy_kev, theta_deg)
