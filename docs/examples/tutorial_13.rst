Detector Response Matrices
===========================
The examples so far use ``Material.transmission()``/``Material.absorption()``,
which give a single number per energy: the fraction of photons absorbed.
A real detector doesn't just absorb-or-not - a photon's energy can be
lost. Some of the energy can escape as a Compton-scattered photon, or
as K-shell fluorescence, and the electronics add their own resolution
blurring. ``Response.response_matrix()`` models all of this: for a grid of
input photon energies it returns the full probability of measuring each
output (PHA channel) energy, not just whether the photon was absorbed.

The detector's photoelectric/Compton split has to come from the McMaster
cross-section data (``database="mcmaster"``, see the
:doc:`user guide </guide/index>` for the ``Material`` class), and
fluorescence escape is only tabulated for a handful of detector materials
(Ge, Si, NaI, CsI, CdTe, SXS) -- see `roentgen.response.fluorescence`.

.. plot::
    :context:
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()

    from roentgen.absorption import Material, Response

    detector = Material({'Cd': 112.4, 'Te': 127.6}, 2 * u.mm, density=5.85 * u.g / u.cm**3, database='mcmaster')
    window = Material('Be', 1 * u.cm, database='mcmaster')
    resp = Response(window, detector)

    # a simple, illustrative energy-resolution model
    fwhm = lambda energy: (0.5 + 0.02 * np.sqrt(energy.to_value(u.keV))) * u.keV

    edges = np.geomspace(4, 150, 513) * u.keV
    area = 3 * u.cm**2
    response_matrix = resp.response_matrix(edges, fwhm=fwhm, area=area)

    fig, ax = plt.subplots()
    smatrix = response_matrix.smatrix.to_value(u.ct * u.cm**2 / (u.photon * u.keV))
    vmax = smatrix.max()
    mesh = ax.pcolormesh(
        response_matrix.ein.to_value(u.keV),
        response_matrix.eout.to_value(u.keV),
        smatrix,
        norm=LogNorm(vmin=vmax * 1e-5, vmax=vmax),
        shading='auto',
    )
    ax.set_xlabel('Photon Energy [keV]')
    ax.set_ylabel('Measured Energy [keV]')
    ax.set_title('CdTe Detector Response Matrix')
    fig.colorbar(mesh, label='Counts / keV per photon cm$^{-2}$')
    plt.show()

The bright diagonal ridge is the photopeak: most photons are fully absorbed
and measured at their true energy. The fainter parallel ridges just below it
are Cd and Te K-shell escape peaks - CdTe has several K-alpha/K-beta lines
close in energy, so they show up as a cluster of ridges rather than one.
Below all of that is the Compton-scatter continuum, bounded above by the
Compton edge -- photons that Compton-scatter in the detector and the
scattered photon escapes, so only part of the energy is measured.

Picking out a single column of this matrix (the response to one input photon
energy) looks like the pulse-height spectrum a real detector would measure
for a monoenergetic source, e.g. a 59.5 keV calibration line (Am-241, a
common choice for CdTe).

.. plot::
    :context: close-figs
    :include-source:

    line_energy = 59.5 * u.keV
    line_index = np.searchsorted(response_matrix.ein.to_value(u.keV), line_energy.to_value(u.keV)) - 1
    eout_center = 0.5 * (response_matrix.eout[:-1] + response_matrix.eout[1:])
    line_cut = response_matrix.smatrix[:, line_index].to_value(u.ct * u.cm**2 / (u.photon * u.keV))

    fig, ax = plt.subplots()
    ax.semilogy(eout_center.to_value(u.keV), line_cut)
    ax.set_xlabel('Measured Energy [keV]')
    ax.set_ylabel('Counts / keV per photon cm$^{-2}$')
    ax.set_title(f'CdTe detector response to a {line_energy} photon')
    ax.set_ylim(1e-4, 3)
    plt.show()

Besides the (resolution-broadened) photopeak at 59.5 keV and the Compton
continuum below it, there's a cluster of escape peaks around 25-38 keV --
59.5 keV minus Cd's and Te's K-shell fluorescence energies - from photons
whose K-shell fluorescence photon escaped the detector instead of being
reabsorbed.

Applying the response to a source
---------------------------------

The response matrix is per unit photon fluence: it includes the detector area
but not the width of the input energy bins, as in SSW's ``resp_calc.pro``. A
photon flux density at the detector multiplied by the input bin widths and
the exposure time is therefore the fluence the matrix expects, and the
result is counts per output bin. The flux density follows from the source
model and its distance (for a point source it falls as :math:`1/r^2`):

.. plot::
    :context: close-figs
    :include-source:

    ein_center = np.sqrt(response_matrix.ein[:-1] * response_matrix.ein[1:])
    distance = 2 * u.au
    flux = 40 * (ein_center / (10 * u.keV)) ** -3 * u.photon / (u.keV * u.s * u.cm**2) * (1 * u.au / distance) ** 2
    exposure = 60 * u.s

    counts = response_matrix.pulse_height_matrix @ (flux * np.diff(response_matrix.ein) * exposure)

    # photons arriving on the detector in each bin, to compare with the counts (the input and
    # output bins are the same here)
    incident = flux * np.diff(response_matrix.ein) * exposure * area

    fig, ax = plt.subplots()
    ax.stairs(incident.to_value(u.photon), response_matrix.ein.to_value(u.keV), label='incident photons')
    ax.stairs(counts.to_value(u.ct), response_matrix.eout.to_value(u.keV), label='detected counts')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Measured Energy [keV]')
    ax.set_ylabel('Photons, counts per bin')
    ax.set_title('CdTe detector counts for a 60 s exposure')
    ax.legend()
    plt.show()

Matching SSW with ``ssw_compatible``
------------------------------------

``response_matrix`` corrects a issue in ``resp_calc.pro`` by default.
Pass ``ssw_compatible=True`` to reproduce SSW instead, for
example to check results against an existing SSW analysis. The two differ in
four ways (see `roentgen.absorption.Response.response_matrix`):

* ``klein.pro`` has an extra factor of :math:`E'/E` in the Compton
  cross-section, which puts too little weight on large-angle scatters;
* the last Compton output bin is oversized when the maximum electron energy is
  above the last output edge;
* the window transmission is applied to the photopeak and escape peaks but not
  to the Compton continuum;
* the photopeak and escape peaks of photons outside the output range spill
  into the neighbouring column (or the first output bin), rather than being
  dropped.

The last two are easiest to see with a hard spectrum, a window that blocks most
low energy photons, and an output range that ends below the highest input
energy. Here the input grid goes up to 150 keV but the output only to about 100 keV.
Every output edge must coincide with an input edge, so we select the output
edges from the input grid rather than generating an independent grid.
The power law with photon index 1.5 is chosen for display purposes, it is much
harder than most real sources, so is compared with one with an index of 4, more
typical of a solar flare:

.. plot::
    :context: close-figs
    :include-source:

    edges_out = edges[edges <= 100 * u.keV]
    rm_fixed = resp.response_matrix(edges, edges_out, fwhm=fwhm, area=area)
    rm_ssw = resp.response_matrix(edges, edges_out, fwhm=fwhm, area=area, ssw_compatible=True)
    eout_center = np.sqrt(edges_out[:-1] * edges_out[1:])

    fig, axes = plt.subplots(2, 2, sharex=True, height_ratios=[3, 1], figsize=(10, 5), layout='constrained')
    for (ax, ax_ratio), index, name in zip(axes.T, [1.5, 4], ['hard (display)', 'solar-like']):
        power_law = lambda energy: 40 * (energy / (10 * u.keV)) ** -index * u.photon / (u.keV * u.s * u.cm**2)
        fluence = power_law(ein_center) * np.diff(edges) * exposure
        counts_fixed = rm_fixed.pulse_height_matrix @ fluence
        counts_ssw = rm_ssw.pulse_height_matrix @ fluence

        # photons arriving on the detector in each output bin (the input grid extends higher)
        incident = power_law(eout_center) * np.diff(edges_out) * exposure * area

        ax.stairs(incident.to_value(u.photon), edges_out.to_value(u.keV), color='gray', ls=':', label='incident photons')
        ax.stairs(counts_fixed.to_value(u.ct), edges_out.to_value(u.keV), label='default')
        ax.stairs(counts_ssw.to_value(u.ct), edges_out.to_value(u.keV), ls='--', label='ssw_compatible=True')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('Photons, counts per bin')
        ax.set_title(f'{name}, index {index}, 60 s exposure')
        ax.legend()

        ratio = counts_ssw.to_value(u.ct) / counts_fixed.to_value(u.ct)
        ax_ratio.stairs(ratio, edges_out.to_value(u.keV), color='k')
        ax_ratio.axhline(1, color='gray', lw=0.5)
        ax_ratio.set_ylabel('SSW / default')
        ax_ratio.set_xlabel('Measured Energy [keV]')
    plt.show()

Over most of the range the two agree to well under a percent, since the
photopeak dominates and is the same in both. They differ at the two ends:

* In the lowest channels SSW counts more. At 4.5 keV it is about ten times
  the default for the hard spectrum, and about 35% higher for the solar-like
  one. Most of this is the photopeaks of the photons above the output
  range, which SSW places at the bottom of the output range instead of dropping
  them (see the single photon response below), so it scales with the number of
  photons above 100 keV and is far smaller for a steep spectrum. The effect is
  large because the Be window has removed almost all the real counts at these
  energies. The unattenuated Compton continuum adds a smaller excess on top,
  from low energy photons that the window stops before they reach the detector.
* In the last channel SSW counts about 35% more for the hard spectrum and 13%
  more for the solar-like spectrum. The size and sign of the difference near
  the upper boundary depend on the energy grids, since SSW redistributes
  out-of-range peaks instead of dropping them.

The spilled photopeaks are clearer in the response to a single photon above the
output range. The default has no photopeak for a 120 keV photon, only its
Compton continuum (ending at the Compton edge near 38 keV) and the K-shell
escape peaks near 88-97 keV that land in range. SSW adds the
photopeak near the bottom of the output range, around 5 keV, so its column
carries about ten times as many counts as the default's:

.. plot::
    :context: close-figs
    :include-source:

    high_index = np.searchsorted(edges.to_value(u.keV), 120) - 1
    unit = u.ct * u.cm**2 / u.photon

    fig, ax = plt.subplots()
    ax.semilogy(eout_center.to_value(u.keV), rm_fixed.eloss_matrix[:, high_index].to_value(unit), label='default')
    ax.semilogy(eout_center.to_value(u.keV), rm_ssw.eloss_matrix[:, high_index].to_value(unit), '--',
                label='ssw_compatible=True')
    ax.set_xlabel('Measured Energy [keV]')
    ax.set_ylabel('Counts per photon cm$^{-2}$')
    ax.set_title('Response to a 120 keV photon, above the output range')
    ax.legend()
    plt.show()
