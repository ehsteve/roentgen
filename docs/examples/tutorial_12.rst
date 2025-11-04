Ba-133 Source Spectrum
======================
In this example, we show the spectrum of a Ba-133 radiation source as measured
by a detector with a finite energy resolution.


.. plot::
    :include-source:

    import numpy as np

    import astropy.units as u
    from astropy.modeling.models import Gaussian1D
    from astropy import visualization

    from matplotlib import pyplot as plt
    from roentgen.nuclides import Nuclide

    ba133 = Nuclide('Ba', 133)
    lines = ba133.get_lines(1 * u.keV, 100 * u.keV, min_intensity=1)
    energy_resolution = 0.5 * u.keV

    g0 = Gaussian1D(amplitude=lines[0]['intensity'], mean=lines[0]['energy'], stddev=energy_resolution)
    for this_line in lines[1:]:
        g0 += Gaussian1D(amplitude=this_line['intensity'], mean=this_line['energy'], stddev=energy_resolution)

    energy_ax = np.arange(1, 85, 0.5) * u.keV

    with visualization.quantity_support():
        fig, ax = plt.subplots()
        ax.plot(energy_ax, g0(energy_ax))
        ax.set_xlabel(f'Energy')
        ax.vlines(lines['energy'], ymin=[0]*len(lines), ymax=lines['intensity'])
        plt.show()


If this emission were detected by a CdTe detector then escape lines would also be present.
These emission lines are at energies which are downshifted by the energy absorbed by the material.
In this example, we will only select the strongest absorption line.

.. plot::
    :include-source:

    import numpy as np
    from matplotlib import pyplot as plt

    from astropy.modeling.models import Gaussian1D
    from astropy import visualization
    import astropy.units as u

    from roentgen.nuclides import Nuclide
    from roentgen.lines import get_lines

    ba133 = Nuclide('Ba', 133)
    lines = ba133.get_lines(1 * u.keV, 100 * u.keV, min_intensity=1)
    energy_resolution = 0.5 * u.keV

    cd_line = get_lines(10 * u.keV, 100 * u.keV, 'Cd', min_intensity=100)[0]
    te_line = get_lines(10 * u.keV, 100 * u.keV, 'Te', min_intensity=100)[0]

    g0 = Gaussian1D(amplitude=lines[0]['intensity'], mean=lines[0]['energy'], stddev=energy_resolution)
    g0 += Gaussian1D(amplitude=lines[0]['intensity'] * 0.2, mean=lines[0]['energy'] - cd_line['energy'], stddev=energy_resolution)
    g0 += Gaussian1D(amplitude=lines[0]['intensity'] * 0.2, mean=lines[0]['energy'] - te_line['energy'], stddev=energy_resolution)

    for this_line in lines[1:]:
        g0 += Gaussian1D(amplitude=this_line['intensity'], mean=this_line['energy'], stddev=energy_resolution)
        g0 += Gaussian1D(amplitude=this_line['intensity'] * 0.2, mean=this_line['energy'] - cd_line['energy'], stddev=energy_resolution)
        g0 += Gaussian1D(amplitude=this_line['intensity'] * 0.2, mean=this_line['energy'] - te_line['energy'], stddev=energy_resolution)

    energy_ax = np.arange(1, 85, 0.5) * u.keV

    with visualization.quantity_support():
        fig, ax = plt.subplots()
        ax.plot(energy_ax, g0(energy_ax))
        ax.set_xlabel(f'Energy')
        ax.set_xlim(1 * u.keV)
        ax.vlines(lines['energy'], ymin=[0]*len(lines), ymax=lines['intensity'], color='red')
        ax.vlines(lines['energy'] - te_line['energy'], ymin=[0]*len(lines), ymax=lines['intensity'] * 0.2, color='green')
        ax.vlines(lines['energy'] - cd_line['energy'], ymin=[0]*len(lines), ymax=lines['intensity'] * 0.2, color='green')
        plt.show()
