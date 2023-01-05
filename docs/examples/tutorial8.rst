Finding the transmission through a stack of materials
=====================================================

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()

    from roentgen.absorption import Material

    optical_path = Material('air', 2 * u.m) + Material('mylar', 5 * u.micron) + Material('Al', 5 * u.micron)

    energy = u.Quantity(np.linspace(1, 100, 300), 'keV')

    plt.plot(energy, optical_path.transmission(energy), label='Transmission')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.ylabel('Response')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()
