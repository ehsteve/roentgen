Finding all x-ray lines in an energy range
==========================================

.. plot::
    :include-source:

    from roentgen.lines import get_lines
    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support
    quantity_support()

    energy_range = [4, 4.5] * u.keV
    lines = get_lines(energy_range[0], energy_range[1])

    ax = plt.subplot(111)
    ax.set_xlim(energy_range[0], energy_range[1])
    ax.set_ylim(0, 1)
    ax.plot(energy_range, [0, 0])

    for this_element, line_name, this_e in zip(lines['element'], lines['transition'], lines['energy']):
        ax.axvline(this_e, label=f'{this_element} {line_name} {this_line}')
    ax.set_xlabel('Energy [keV]')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
