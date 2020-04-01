Finding the x-ray emission lines for a particular element
=========================================================

.. plot::
    :include-source:

    from roentgen.lines import get_lines_for_element
    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support
    quantity_support()

    this_element = 'Ni'

    plt.plot([0, 10], [0, 0])
    plt.xlim(0, 10)
    plt.ylim(0, 1)

    lines = get_lines_for_element(this_element)
    for this_trans, this_energy in zip(lines['transition'], lines['energy']):
        plt.axvline(this_energy, label=f'{this_element} {this_energy} {this_trans}')
    plt.legend()
    plt.show()
