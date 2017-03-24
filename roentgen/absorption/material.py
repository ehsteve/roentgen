"""
"""

from __future__ import absolute_import

import numpy as np
import os
import astropy.units as u
from scipy import interpolate
import roentgen

__all__ = ['Material', 'MassAttenuationCoefficient', 'Compound', 'is_an_element', 'get_atomic_number']

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory


class Material(object):
    """An object which provides the properties of a material in x-rays

    Parameters
    ----------
    material_str : str or int
        A representation of the material which includes  (e.g. cdte, be, mylar, si)
    thickness : `astropy.units.Quantity`
        The thickness of the material in the optical path.
    density : `astropy.units.Quantity`
        The density of the material. If None use default values.

    Examples
    --------
    >>> from roentgen.absorption.material import Material
    >>> import astropy.units as u
    >>> detector = Material('cdte', 500 * u.um)
    >>> thermal_blankets = Material('mylar', 0.5 * u.mm)
    """

    @u.quantity_input
    def __init__(self, material_str, thickness: u.m, density=None):
        self.name = material_str
        self.thickness = thickness
        self.mass_attenuation_coefficient = MassAttenuationCoefficient(material_str)
        self.name = self.mass_attenuation_coefficient.name
        if density is None:
            if is_an_element(material_str):
                self.density = roentgen.elemental_densities[get_atomic_number(material_str)-1]['density']
            else:
                # not using loc because table indexing is not yet stable
                # self.density = roentgen.compounds.loc[material_str]['density']
                index = list(roentgen.compounds['symbol']).index(material_str)
                self.density = roentgen.compounds[index]['density']
        else:
            self.density = density

    def __repr__(self):
        """Returns a human-readable representation."""
        txt = '<Material ' + str(self.name) + ' (' + str(self.name) + ') '
        txt += str(self.thickness) + ' ' + str(self.density) + '>'
        return txt

    def __add__(self, other):
        if isinstance(other, Material):
            return Compound([self, other])
        elif isinstance(other, Compound):
            return Compound([self] + other.materials)
        else:
            raise ValueError("Can't add <Material> and " + str(other))

    def transmission(self, energy):
        """Provide the transmission fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        coefficients = self.mass_attenuation_coefficient.func(energy)
        transmission = np.exp(- coefficients * self.density * self.thickness)
        return transmission

    def absorption(self, energy):
        """Provides the absorption fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return 1 - self.transmission(energy)


class Compound(object):
    """An object which provides the x-ray properties of a compound (i.e.
     many materials).

    Parameters
    ----------
    materials : list
        A list of Material objects

    Examples
    --------
    >>> from roentgen.absorption.material import Material, Compound
    >>> import astropy.units as u
    >>> detector = Compound([Material('Pt', 5 * u.um), Material('cdte', 500 * u.um)])
    """

    def __init__(self, materials):
        self.materials = materials

    def __add__(self, other):
        if isinstance(other, Material):
            return Compound(self.materials + [other])
        elif isinstance(other, Compound):
            return Compound(self.materials + other.materials)
        else:
            raise ValueError("Can't add <Compound> and " + str(other))

    def __repr__(self):
        """Returns a human-readable representation."""
        txt = '<Compound '
        for material in self.materials:
            txt += str(material)
        return txt + '>'

    def transmission(self, energy):
        """Provide the transmission fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        transmission = np.ones(len(energy), dtype=np.float)
        for material in self.materials:
            coefficients = material.mass_attenuation_coefficient.func(energy)
            transmission *= np.exp(- coefficients * material.density * material.thickness)
        return transmission

    def absorption(self, energy):
        """Provides the absorption fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return 1 - self.transmission(energy)


class MassAttenuationCoefficient(object):
    """
    The mass attenuation coefficient

    Parameters
    ----------
    material : str
        A string representing a material (e.g. cdte, be, mylar, si)
    """
    def __init__(self, material):
        if is_an_element(material):
            atomic_number = get_atomic_number(material)
            datafile_path = os.path.join(_data_directory, 'elements', 'z' + str(atomic_number).zfill(2) + '.csv')
            symbol = roentgen.elements[atomic_number-1]['symbol']
            name = roentgen.elements[atomic_number-1]['name']
        else:
            datafile_path = os.path.join(_data_directory, 'compounds_mixtures', material.lower().replace(' ', '_') + '.csv')
            index = list(roentgen.compounds['symbol']).index(material)
            symbol = roentgen.compounds[index]['symbol']
            name = roentgen.compounds[index]['name']

        data = np.loadtxt(datafile_path, delimiter=',')
        # find the material in our list
        self.symbol = symbol
        self.name = name
        self.energy = u.Quantity(data[:, 0] * 1000, 'keV')
        self.data = u.Quantity(data[:, 1], 'cm^2/g')

        data_energy_kev = np.log10(self.energy.value)
        data_attenuation_coeff = np.log10(self.data.value)
        self._f = interpolate.interp1d(data_energy_kev, data_attenuation_coeff, bounds_error=False, fill_value=0.0)
        self.func = lambda x: u.Quantity(10 ** self._f(np.log10(x.to('keV').value)), 'cm^2/g')

    def get_filename(material_str):
        if len(material_str) <= 2:
            # likely a symbol of an element
            if material_str in list(elements['symbol']):
                return


def is_an_element(element_str):
    """Returns True is the string represents an element"""
    result = False
    if (len(element_str) <= 2) and (element_str in list(roentgen.elements['symbol'])):
            result = True
    else:
        if element_str in list(roentgen.elements['name']):
            result = True
    return result


def get_atomic_number(element_str):
    """Return the atomic number of the element"""
    # check to see if element_str is symbol
    if is_an_element(element_str):
        if len(element_str) <= 2:
            atomic_number = list(roentgen.elements['symbol']).index(element_str) + 1
        else:
            atomic_number = list(roentgen.elements['name']).index(element_str) + 1
    else:
        atomic_number = None
    return atomic_number
