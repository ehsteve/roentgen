Finding X-ray K absorption Edges
================================

.. plot::
    :include-source:

    import numpy as np
    from matplotlib import pyplot as plt
    from astropy import constants as const
    import astropy.units as u

    from roentgen.absorption import Material
    from roentgen.lines import get_edges

    materials = ['Mo', 'Ag']
    e = np.linspace(5, 30, 1000)*u.keV

    for this_material in materials:
        mat = Material(this_material, 50 * u.micron)
        edges = get_edges(this_material)

        plt.plot(e, mat.transmission(e), label=f"{this_material} K edge at {edges['energy'][0].to('keV')}")
        # lets only show the K shell
        plt.axvline(edges['energy'][0].to('keV').value, linestyle='dashed')

    plt.xlabel(f'Energy [{e.unit}]')
    plt.ylabel('Transmission')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()
