Finding the x-ray transmission through different materials
==========================================================

.. plot::
    :include-source:

    from roentgen.absorption import Material
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support
    quantity_support()

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
