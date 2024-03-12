Defining your own custom material
=================================
Though the compounds list provided in this package is extensive it may not provide the exact material that you are interested in.
It easily possible to define your own material if you know it's chemical composition.
The following example, we investigate the absorption of a Gallium Arsenide (GaAs) detector in the same configuration as that in `G. Lioliou & A.M. Barnett (2016) <https://doi.org/10.1016/j.nima.2016.08.047>`__.
Compare our plot with their `Figure 1 <https://www.sciencedirect.com/science/article/pii/S016890021630866X#f0005>`__.

.. plot::
    :include-source:

    import astropy.units as u
    import numpy as np

    from matplotlib import pyplot as plt
    from roentgen.absorption import Material

    gaas_deadlayer = Material({'Ga': 0.518, 'As': 0.482},
                            density=5.32*u.g/u.cm**3, thickness=500*u.nm)

    gaas_detector = Material({'Ga': 0.518, 'As': 0.482},
                            density=5.32*u.g/u.cm**3, thickness=10000*u.nm)
    e = np.linspace(1, 30, 1000)*u.keV

    plt.plot(e, gaas_detector.absorption(e) * gaas_deadlayer.transmission(e))
    plt.xlabel(f'Energy [{e.unit}]')
    plt.ylabel('Quantum Efficiency')
    plt.ylim(0.001, 1)
    plt.yscale('log')
    plt.show()
