import numpy as np

import astropy.units as u
from astropy.constants.codata2018 import h, c, e, k_B

__all__ = [
    "thermal_bremsstrahlung",
    "mewe_gaunt",
    "xray_tube_spectrum",
    "cutoff_wavelength"
]


def thermal_bremsstrahlung(wavelength, z, density, temperature):
    """
    This equation assumes that the density of ions is equal to the density of
    electrons.
    """
    frequency = wavelength.to('Hz', equivalencies=u.spectral())
    return z ** 2 * density ** 2 * np.sqrt(1 / (k_B * temperature)) * \
           np.exp(- h * frequency / (k_B * temperature)) * mewe_gaunt(wavelength, temperature)


def mewe_gaunt(wavelength, temperature):
    """Based on Equation 3 in Mewe 1985."""
    return 27.83 * (temperature.to('MK').value + 0.65) ** (-1.33) + 0.15 * wavelength.to('angstrom').value ** (0.34) * temperature.to('MK').value ** (0.422)


def xray_tube_spectrum(wavelength, voltage):
    """Implements the Kramer's law to estimate the continuum part of the
    spectrum of an x-ray tube"""
    f = (wavelength / cutoff_wavelength(voltage) - 1) * 1 / wavelength ** 2
    # remove unphysical negative values
    f[f < 0.0] = 0.0
    return f


def cutoff_wavelength(voltage):
    """
    Implements the Duane-Hunt law which provides the low energy cutoff for the
    emission of an x-ray tube.
    """
    return h * c / (e * voltage)