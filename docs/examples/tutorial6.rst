Finding the x-ray transmission through different materials
==========================================================

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()

    from roentgen.absorption import Material

    thickness = 5 * u.mm

    material_list = ['air', 'Al', 'Pb']
    energy = u.Quantity(np.arange(1, 100, 0.1), 'keV')

    for material in material_list:
        mat = Material(material, thickness)
        plt.plot(energy, mat.transmission(energy), label=mat.name)

    plt.title(f"Thickness = {thickness}")
    plt.ylabel('Transmission')
    plt.legend(loc='lower right')
    plt.show()
