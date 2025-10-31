X-ray emission from radionuclides
=================================
The purpose of this guide is to present an overview of the `roentgen.nuclides` module which provides access to tables of emission lines from radioactive sources (e.g. Am-241, Ba-133).
radionuclides are unstable atoms that spontaneously and randomly decay into a different nuclide over time, emitting ionizing radiation in the process.
Those new nuclide may undergo a similar process leading to a chain of decay of different elements.
The radiation emitted may taken the form of x-rays or gamma-rays with characteristic energies.
Because these characteristic energies are well-known, these sources are frequently used in laboratories to calibrate x-ray and gamma-ray detectors.

This module provides access to data from the Library for gamma and alpha emissions (`Lara <http://www.lnhb.fr/accueil/donnees-nucleaires/module-lara/>`_) provided by the `Laboratoire National Henri Becquerel <http://www.lnhb.fr>`_.

The list of available nuclides is available through::

    >>> from roentgen.nuclides import nuclides_list
    >>> print(nuclides_list)

To access a particular nuclide first initialize it::

    >>> from roentgen.nuclides import Nuclide
    >>> fe55 = Nuclide('Fe', 55)
    >>> print(fe55)

You can see all of the emission lines::

    >>> fe55.lines

For nuclides with complicated decay chains this may be too much information.
This module therefore provides a convenience function to filter the list.
For example, Am-241 has a long decay chain generating 873 different lines some of them very weak.
We can there filter the list to only the strongest lines within a particular energy range:::

    >>> import astropy.units as u
    >>> am241 = Nuclide('Am', 241)
    >>> strong_lines = am241.get_lines(5 * u.keV, 100 * u.keV, min_intensity=5)
    >>> print(strong_lines)

