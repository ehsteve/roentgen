Finding the x-ray emission lines for a particular element
=========================================================

.. plot::
    :include-source:

    from roentgen.lines import get_lines_for_element
    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support
    quantity_support()

    this_element = 'Ni'

    ax = plt.subplot(111)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    ax.plot([0, 10], [0, 0])

    lines = get_lines_for_element(this_element)
    for line_name, this_trans in zip(lines['transition'], lines['energy']):
        ax.axvline(this_line, label=f'{this_element} {line_name} {this_trans}', color=this_color)
    plt.legend()
    plt.show()
