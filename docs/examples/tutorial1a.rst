Getting and plotting the Mass Attenuation Coefficient
=====================================================

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    import astropy.units as u
    from astropy.visualization import quantity_support
    quantity_support()
    from roentgen.absorption import MassAttenuationCoefficient

    cdte_atten = MassAttenuationCoefficient('cdte')

    energy = u.Quantity(np.arange(1, 1000), 'keV')
    atten = cdte_atten.func(energy)

    plt.plot(energy, atten)
    plt.plot(cdte_atten.energy, cdte_atten.data, 'o')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.ylabel('Mass attenuation Coefficient [' + str(atten.unit) + ']')
    plt.title(cdte_atten.name)
    plt.show()
