Ba-133 Source Spectrum
======================
In this example, we show the spectrum of a Ba-133 radiation source as measured
by a detector with a finite energy resolution.


.. plot::
    :include-source:

    import astropy.units as u
    import numpy as np

    from matplotlib import pyplot as plt
    from roentgen.nuclides import Nuclide

    ba133 = Nuclide('Ba', 133)
    lines = ba133.get_lines(1 * u.keV, 100 * u.keV, min_intensity=1)
    energy_resolution = 0.5 * u.keV

    g0 = Gaussian1D(amplitude=lines[0]['intensity'], mean=lines[0]['energy'], stddev=energy_resolution)
    for this_line in lines[1:]:
        g0 += Gaussian1D(amplitude=this_line['intensity'], mean=this_line['energy'], stddev=energy_resolution)

    energy_ax = np.arange(1, 100, 0.5) * u.keV

    fig, ax = plt.subplots()
    ax.plot(energy_ax.value, g0(energy_ax))
    ax.set_xlabel(f'Energy [{energy_ax.unit}]')
    ax.vlines(lines['energy'].value, ymin=[0]*len(lines), ymax=lines['intensity'])
    plt.show()

