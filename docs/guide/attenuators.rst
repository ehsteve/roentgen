Attenuators
===========

RHESSI
------
The Reuven Ramaty High Energy Solar Spectroscopic Imager (RHESSI) makes use of two aluminum attenuators or shutters to reduce the high fluxes from solar flares at low energies.
These two attenuators are frequently referred to as the thin and thick attenuator.
They can also be inserted at the same time creating a third attenuation scenario.
The attenuators are described in `Smith et al. 2002 <https://link.springer.com/article/10.1023/A:1022400716414>`_.
The shutter parameters can be found at `HESSI SHUTTER PARAMETERS <https://hesperia.gsfc.nasa.gov/rhessidatacenter/instrument/Shutter_parameters.html>`_.
The shutters are circular aluminum disks and have three different thicknesses.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    import roentgen
    from roentgen.absorption import Material
    from roentgen.util import get_material_density

    # the following parameters were found on the shutter parameter page listed above
    shutter_diameters = {'thick': [1.716, 24.3, 58.5] * u.mm,
                         'thin': [6.74, 17.1, 58.5] * u.mm}
    shutter_thickness = {'thick': [0.013, 0.177, 0.414] * u.g / (u.cm)**2 / get_material_density('Al'),
                         'thin': [0.0192, 0.261, 0.112] * u.g / (u.cm)**2 / get_material_density('Al')}

    # calculate the fraction of the shutter with each thickness
    shutter_fractions = {'thick': u.Quantity([(shutter_diameters['thick'][0]) ** 2,
                                      shutter_diameters['thick'][1]**2 - shutter_diameters['thick'][0] ** 2,
                                      shutter_diameters['thick'][2]**2 - shutter_diameters['thick'][1] ** 2]) / shutter_diameters['thick'][2] ** 2,
                         'thin': u.Quantity([(shutter_diameters['thin'][0]) ** 2,
                                      shutter_diameters['thin'][1]**2 - shutter_diameters['thin'][0] ** 2,
                                      shutter_diameters['thin'][2]**2 - shutter_diameters['thin'][1] ** 2]) / shutter_diameters['thin'][2] ** 2}

    def composite_attenuator(energy, thicknesses, fractions, material='Al'):
        """An attenuator with multiple thicknesses and fractions."""
        transmission_fraction = 0.0
        if np.sum(fractions) != 1:
            raise AttributeError('Fractions must add up to 1.')
        for this_thick, this_frac in zip(thicknesses, fractions):
            transmission_fraction += this_frac * Material(material, thickness=this_thick).transmission(energy)
        return transmission_fraction

    energy = np.linspace(3, 30, 100) * u.keV

    plot_type = ['transmission', 'flux']

    fig, axis = plt.subplots(len(plot_type), figsize=(10, 20))
    for this_plot_type, ax in zip(plot_type, axis):

        for this_key in list(shutter_diameters.keys()) + ['thick+thin']:
            if this_key != 'thick+thin':
                this_trans = composite_attenuator(energy, shutter_thickness[this_key], shutter_fractions[this_key])
            else:
                this_trans = composite_attenuator(energy, shutter_thickness['thick'], shutter_fractions['thick']) * composite_attenuator(energy, shutter_thickness['thin'], shutter_fractions['thin'])
            if this_plot_type == 'flux':
                this_trans *= np.exp(-energy.value)
            ax.plot(energy,  this_trans, label=this_key)

        ax.legend()
        ax.set_title(this_plot_type)
        if this_plot_type == 'flux':
            ax.set_yscale('log')
            ax.set_ylim(1e-15, 1e-5)

    plt.show()