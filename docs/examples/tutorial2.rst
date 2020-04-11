Finding the x-ray absorption of different detector materials
============================================================

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    import astropy.units as u
    from roentgen.absorption import Material

    thickness = 500 * u.micron

    material_list = ['cdte', 'Si', 'Ge']
    energy = u.Quantity(np.arange(1, 100, 0.1), 'keV')

    for material in material_list:
        mat = Material(material, thickness)
        plt.plot(energy, mat.absorption(energy), label=mat.name)

    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.ylabel('Absorption')
    plt.legend(loc='lower left')
    plt.show()
