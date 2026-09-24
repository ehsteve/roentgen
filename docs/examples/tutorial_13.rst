Detector Response Matrices
===========================
The examples so far use ``Material.transmission()``/``Material.absorption()``,
which give a single number per energy: the fraction of photons absorbed.
A real detector doesn't just absorb-or-not - a photon's energy can be
mismeasured. Some of the energy can escape as a Compton-scattered photon, or
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

    # a simple, illustrative energy-resolution model: FWHM in keV
    fwhm = lambda energy: 0.5 + 0.02 * np.sqrt(energy)

    response_matrix = resp.response_matrix(
        512, 4 * u.keV, 150 * u.keV, fwhm=fwhm, area=3 * u.cm**2
    )

    fig, ax = plt.subplots()
    vmax = response_matrix.smatrix.max()
    mesh = ax.pcolormesh(
        response_matrix.ein.to_value(u.keV),
        response_matrix.eout.to_value(u.keV),
        response_matrix.smatrix,
        norm=LogNorm(vmin=vmax * 1e-5, vmax=vmax),
        shading='auto',
    )
    ax.set_xlabel('Photon Energy [keV]')
    ax.set_ylabel('Measured Energy [keV]')
    ax.set_title('CdTe Detector Response Matrix')
    fig.colorbar(mesh, label='Counts / keV / photon')
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
common choice for CdTe). No need to recompute ``response_matrix`` -- with
``:context:`` set, this picks up right where the previous plot left off:

.. plot::
    :context: close-figs
    :include-source:

    line_energy = 59.5 * u.keV
    line_index = np.searchsorted(response_matrix.ein.to_value(u.keV), line_energy.to_value(u.keV)) - 1
    eout_center = 0.5 * (response_matrix.eout[:-1] + response_matrix.eout[1:])

    fig, ax = plt.subplots()
    ax.semilogy(eout_center, response_matrix.smatrix[:, line_index])
    ax.set_xlabel('Measured Energy')
    ax.set_ylabel('Counts / keV / photon')
    ax.set_title(f'CdTe detector response to a {line_energy} photon')
    ax.set_ylim(1e-4, 3)
    plt.show()

Besides the (resolution-broadened) photopeak at 59.5 keV and the Compton
continuum below it, there's a cluster of escape peaks around 25-38 keV --
59.5 keV minus Cd's and Te's K-shell fluorescence energies - from photons
whose K-shell fluorescence photon escaped the detector instead of being
reabsorbed.
