"""
"""

from __future__ import absolute_import

import numpy as np
import os
import astropy.units as u
from scipy import interpolate
import roentgen
from roentgen.util import get_atomic_number, get_density, get_compound_index, get_element_symbol, is_an_element, is_in_known_compounds

__all__ = [
    "Material",
    "MassAttenuationCoefficient",
    "Compound",
    "Response"
]

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory


class Material(object):
    """An object which provides the properties of a material in x-rays

    Parameters
    ----------
    material_str : str
        A string representation of the material which includes an element symbol
        (e.g. Si), an element name (e.g. Silicon), or the name of a compound
        (e.g. cdte, mylar). Cap-sensitive.
    thickness : `astropy.units.Quantity`
        The thickness of the material in the optical path.
    density : `astropy.units.Quantity`
        The density of the material. If None use default values.
    is_detector : bool
        A property to hold information on whether absorption or transmission
        is important. Used for compounds.

    Examples
    --------
    >>> from roentgen.absorption.material import Material
    >>> import astropy.units as u
    >>> detector = Material('cdte', 500 * u.um)
    >>> thermal_blankets = Material('mylar', 0.5 * u.mm)
    """

    @u.quantity_input
    def __init__(self, material_str, thickness: u.m, density=None,
                 is_detector=False):
        self.name = material_str
        self.thickness = thickness
        self.mass_attenuation_coefficient = MassAttenuationCoefficient(material_str)
        self.name = self.mass_attenuation_coefficient.name
        self.is_detector = is_detector

        if density is None:
            self.density = get_density(material_str)
        else:
            self.density = density

    def __repr__(self):
        """Returns a human-readable representation."""
        txt = "<Material " + str(self.name) + " (" + str(self.name) + ") "
        txt += str(self.thickness) + " " + str(self.density) + ">"
        return txt

    def __add__(self, other):
        if isinstance(other, Material):
            return Compound([self, other])
        elif isinstance(other, Compound):
            return Compound([self] + other.materials)
        else:
            raise ValueError("Can't add <Material> and " + str(other))

    @u.quantity_input(energy=u.keV)
    def transmission(self, energy):
        """Provide the transmission fraction.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        coefficients = self.mass_attenuation_coefficient.func(energy)
        transmission = np.exp(-coefficients * self.density * self.thickness)
        return transmission.value  # remove the dimensionless unit

    @u.quantity_input(energy=u.keV)
    def absorption(self, energy):
        """Provides the absorption fraction.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return 1.0 - self.transmission(energy)


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
        self.name = "".join(
            [m.name + " " + str(m.thickness) + " " for m in materials]
        ).rstrip()

    def __add__(self, other):
        if isinstance(other, Material):
            return Compound(self.materials + [other])
        elif isinstance(other, Compound):
            return Compound(self.materials + other.materials)
        else:
            raise ValueError("Can't add <Compound> and " + str(other))

    def __repr__(self):
        """Returns a human-readable representation."""
        txt = "<Compound "
        for material in self.materials:
            txt += str(material)
        return txt + ">"

    def transmission(self, energy):
        """Provide the transmission fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        transmission = np.ones(len(energy), dtype=np.float)
        for material in self.materials:
            this_transmission = (
                material.transmission(energy)
            )
            transmission *= this_transmission
        return transmission

    def absorption(self, energy):
        """Provides the absorption fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return 1.0 - self.transmission(energy)


class Response(object):
    """
    An object to handle the response of a detector material which includes
    an optical path or filter through which x-rays must first traverse before
    reaching the detector.

    Parameters
    ----------
    optical_path : list
        A list of Material objects which make up the optical path.

    detector : Material or None
        A Material which represents the detector material where the x-rays
        are absorbed. If provided with None, than assume a perfectly absorbing
        detector material.

    Examples
    --------
    >>> from roentgen.absorption.material import Material, Response
    >>> import astropy.units as u
    >>> optical_path = [Material('air', 1 * u.m), Material('Al', 500 * u.mm)]
    >>> resp = Response(optical_path, detector=Material('cdte', 500 * u.um))
    """
    def __init__(self, optical_path, detector):
        # make sure the materials are a list since we iterate over them
        # to calculate the transmission
        if isinstance(optical_path, Material):
            self.optical_path = [optical_path]
        elif isinstance(optical_path, list):
            self.optical_path = optical_path
        else:
            raise TypeError("optical_path must be Material or list of Materials")
        if (type(detector) is Material) or (detector is None):
            self.detector = detector
        else:
            raise TypeError('Detector must be a Material or None')

    def __repr__(self):
        """Returns a human-readable representation."""
        txt = "<Response optical path="
        for material in self.optical_path:
            txt += str(material)
        txt += " detector=" + str(self.detector)
        return txt + ">"

    def response(self, energy):
        """Returns the response as a function of energy"""
        # calculate the transmission
        transmission = np.ones(len(energy), dtype=np.float)
        detector_absorption = np.ones(len(energy), dtype=np.float)
        for material in self.optical_path:
            this_transmission = (
                material.transmission(energy)
            )
            transmission *= this_transmission
        if self.detector is None:
            detector_absorption = np.ones(len(energy), dtype=np.float)
        else:
            detector_absorption = self.detector.absorption(energy)

        return transmission * detector_absorption


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
            datafile_path = os.path.join(
                _data_directory, "elements", "z" + str(atomic_number).zfill(2) + ".csv"
            )
            symbol = roentgen.elements[atomic_number - 1]["symbol"]
            name = roentgen.elements[atomic_number - 1]["name"]
        elif is_in_known_compounds(material):
            compound_index = get_compound_index(material)
            symbol = roentgen.compounds[compound_index]["symbol"]
            name = roentgen.compounds[compound_index]["name"]
            datafile_path = os.path.join(
                _data_directory, "compounds_mixtures", symbol.replace(" ", "_") + ".csv"
            )
        else:
            return NameError("Element or compound not found.")
        data = np.loadtxt(datafile_path, delimiter=",")
        # find the material in our list
        self.symbol = symbol
        self.name = name
        self.energy = u.Quantity(data[:, 0] * 1000, "keV")
        self.data = u.Quantity(data[:, 1], "cm^2/g")

        data_energy_kev = np.log10(self.energy.value)
        data_attenuation_coeff = np.log10(self.data.value)
        self._f = interpolate.interp1d(
            data_energy_kev, data_attenuation_coeff, bounds_error=False,
            fill_value=0.0
        )
        self.func = lambda x: u.Quantity(
            10 ** self._f(np.log10(x.to("keV").value)), "cm^2/g"
        )

    def get_filename(material_str):
        if len(material_str) <= 2:
            # likely a symbol of an element
            if material_str in list(roentgen.elements["symbol"]):
                return


