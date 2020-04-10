Finding all x-ray lines in an energy range
==========================================

.. plot::
    :include-source:

    from roentgen.lines import get_lines
    from roentgen.absorption import get_atomic_number
    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()

    energy_range = [4, 4.2] * u.keV
    lines = get_lines(energy_range[0], energy_range[1])

    for row in lines:
        plt.vlines([row['energy']], [0], row['intensity'], label=f'{row["z"]} {row["transition"]}')

    plt.ylabel('intensity')
    plt.legend(loc='upper left')
    plt.show()
