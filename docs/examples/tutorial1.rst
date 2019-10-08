Getting and plotting the Mass Attenuation Coefficient
=====================================================

.. plot::
    :include-source:

    from roentgen.absorption import MassAttenuationCoefficient
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt

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
