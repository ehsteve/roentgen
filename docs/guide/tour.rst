Tour
====
The purpose of this tour is to present an overview of the functionality of this
package.

Absorption
----------
The purpose of this module is to provide straight-forward analysis of x-ray absorption and transmission. The primary
component that mediates the x-ray attenuation through a material is its mass attenuation coefficient. These tabulated
values can be inspected using the :ref:`MassAttenuationCoefficient` object. To create one::

    from roentgen.absorption import MassAttenuationCoefficient
    si_matten = MassAttenuationCoefficient('Si')

Tabulated values for all elements are provided as well as additional specialized materials. Elements can be specificied
by their symbol name or by their full name (e.g. Si, Silicon). The first letter must be capitalized. A list of all of
the elements is provided by::

    roentgen.elements

Specialized materials, referred to as compounds, are also available. A complete list is provided by::

    roentgen.compounds

For compounds all letters are lower case. Here is the mass attenuation coefficient for Silicon.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt
    from roentgen.absorption import MassAttenuationCoefficient

    si_matten = MassAttenuationCoefficient('Si')
    plt.plot(si_matten.energy, si_matten.data)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Energy [' + str(si_matten.energy[0].unit) + ']')
    plt.ylabel('Mass Attenuation Coefficient [' + str(si_matten.data[0].unit) + ']')
    plt.title(si_matten.name)

In order to determine the x-ray attenuation by a material the :ref:`Material` object is provided. This object can easily be
created by providing the thickness of the material through which the x-rays are interacting. For example, a 500 micron
thick layer of Aluminum can be created by::

    al = Material('Al', 500 * u.micron)

An optional density can also be provided. A default density is assumed if none is provided. To inspect the density::

    al.density

The assumed densities are stored in `roentgen.elemental_densities` and in `roentgen.compounds`. Once this object is
created it is now possible to get the absorption and transmission::

    energy = u.Quantity(np.arange(1,30), 'keV')
    al.transmission(energy)
    al.absoprtion(energy)

Here is a plot of the absorption of x-rays through 500 micron of Aluminum, a standard thickness for electronics boxes.

.. plot::
    :include-source:

    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    from roentgen.absorption import Material

    al = Material('Al', 500 * u.micron)
    energy = u.Quantity(np.arange(1,30), 'keV')

    plt.plot(energy, al.absorption(energy), label='Absorption')
    plt.ylabel('Efficiency')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.title(al.name)
    plt.legend(loc='lower left')


From the above plot, one can see that the Al blocks most all x-rays below about 7 keV. The relationship between
transmission and absorption can be seen in the following plot for 500 microns of Silicon, a
standard thickness for an x-ray detector.

.. plot::
    :include-source:

    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    from roentgen.absorption import Material

    si = Material('Si', 500 * u.micron)
    energy = u.Quantity(np.arange(1, 50), 'keV')

    plt.plot(energy, si.absorption(energy), label='Absorption')
    plt.plot(energy, si.transmission(energy), label='Transmission')
    plt.ylabel('Efficiency')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.title(si.name)
    plt.legend(loc='lower left')


One final plot which shows the transmission of x-rays through 10 meters of air.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt
    from roentgen.absorption import Material
    import numpy as np

    thickness = 10 * u.m
    air = Material('air', thickness)
    energy = u.Quantity(np.arange(1,30), 'keV')

    plt.plot(energy, air.transmission(energy), label='Transmission')
    plt.ylabel('Transmission')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.title("{0} {1}".format(str(thickness), air.name))
    # plt.legend(loc='lower left')

This plot shows that air, though not a dense material, can absorb low energy x-rays over long distances.
Materials can be added together to form more complex optical paths. If two materials are added together they form
a new object, a :ref:`Compound`. A simple example might be to consider the transmission through air and then through a
thermal blanket composed of a thin layer of mylar and Aluminum::

    optical_path = Material('air', 2 * u.m) + Material('mylar', 5 * u.micron) + Material('Al', 5 * u.micron)

This new object provides also provides transmission and absorption of the combination of these materials. Here is a
plot of that transmission over energy

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt
    from roentgen.absorption import Material
    import numpy as np

    optical_path = Material('air', 2 * u.m) + Material('mylar', 5 * u.micron) + Material('Al', 5 * u.micron)
    energy = u.Quantity(np.arange(1,30), 'keV')

    plt.plot(energy, optical_path.transmission(energy), label='Transmission')
    plt.ylabel('Efficiency')
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.legend(loc='upper left')


Frequently it is useful to consider the response function of a particular detector which includes absorption through
materials included before the detector. This can be calculated by multiplying the transmission of the materials before
the detector by the absorption of the detector material. The following example uses the same optical path as defined
above and assumes a Silicon detector.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt
    from roentgen.absorption import Material
    import numpy as np

    optical_path = Material('air', 2 * u.m) + Material('mylar', 5 * u.micron) + Material('Al', 5 * u.micron)
    si = Material('Si', 500 * u.micron)
    energy = u.Quantity(np.arange(1,30), 'keV')

    plt.plot(energy, optical_path.transmission(energy) * si.absorption(energy))
    plt.xlabel('Energy [' + str(energy.unit) + ']')
    plt.ylabel('Response')

This plot shows that the peak efficiency is less than 50% and lies around 15 keV.