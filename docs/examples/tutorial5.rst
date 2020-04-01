Finding the response of an x-ray detector
=========================================

.. plot::
    :include-source:

    from roentgen.absorption import Material, Response
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt

    optical_path = [Material('air', 2 * u.m), Material('mylar', 5 * u.micron), Material('Al', 5 * u.micron)]
    detector = Material('Si', 500 * u.micron)

    resp = Response(optical_path=optical_path, detector=detector)
    energy = u.Quantity(np.arange(1,30), 'keV')

    plt.plot(energy, resp.response(energy))
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.ylabel('Response')
    plt.show()
