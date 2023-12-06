X-ray Transmission through air at STP
======================================

.. plot::
    :include-source:

    import numpy as np
    from matplotlib import pyplot as plt
    from astropy import constants as const
    import astropy.units as u

    from roentgen.absorption import Material
    from roentgen.util import density_ideal_gas

    air_density = density_ideal_gas(1 * const.atm, 20 * u.Celsius)
    air = Material('air', 1 * u.m, density=air_density)

    e = np.linspace(1, 30, 1000)*u.keV

    plt.plot(e, air.transmission(e))
    plt.xlabel(f'Energy [{e.unit}]')
    plt.ylabel('X-ray Transmission through 1 m of air at STP')
    plt.ylim(0, 1)
    plt.show()
