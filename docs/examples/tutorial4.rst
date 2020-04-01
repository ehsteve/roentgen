Finding all x-ray lines in an energy range
==========================================

.. plot::
    :include-source:

    from roentgen.lines import get_lines
    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()

    energy_range = [4, 4.5] * u.keV
    lines = get_lines(energy_range[0], energy_range[1])

    plt.plot(energy_range, [0, 0])

    plt.xlim(energy_range[0], energy_range[1])
    plt.ylim(0, 1)

    for this_element, line_name, this_e in zip(lines['element'], lines['transition'], lines['energy']):
        plt.axvline(this_e, label=f'{this_element} {line_name} {line_name}')
    plt.xlabel('Energy [keV]')
    plt.show()
