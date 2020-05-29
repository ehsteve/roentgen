import numpy as np

import astropy.units as u
from astropy.constants.codata2018 import h, c, e, k_B

__all__ = [
    "ThermalBremsstrahlung",
    "GauntFactor",
    "XraytubeEmission",
    "CutoffWavelength"
]

@u.quantity_input(wavelength=u.nm, equivalencies=u.spectral())
def ThermalBremsstrahlung(wavelength, temperature: u.K):
    """
    Provides the radiation spectrum caused by the acceleration of electrons
    in the Coulomb field of another charge (generally an ion).
    This is generally referred to as free-free emission. For the case of
    thermal bremsstrahlung, the assumption is that everything is in thermal
    equilibrium. Assumes an optically thin plasma.

    Parameters
    ----------
    wavelength : `astropy.units.Quantity`
        The wavelength at which to evaluate the emission.

    temperature : `astropy.units.Quantity`
        The temperature of the plasma.

    Returns
    -------
    emission : The unscaled emission, provided without units

    TODO: Consider changing this function to an astropy model.
    """
    frequency = wavelength.to('Hz')

    emission = np.sqrt(1 / (k_B * temperature))
    emission *= np.exp(- h * frequency / (k_B * temperature.to('K')))
    emission *= GauntFactor(wavelength, temperature)

    return emission.to_value()


@u.quantity_input(wavelength=u.nm, equivalencies=u.spectral())
def GauntFactor(wavelength, temperature: u.K):
    """
    Provides the quantum mechanical correction factor for the bremsstrahlung
    emission process in an optically thin plasma.

    Parameters
    ----------
    wavelength : `astropy.units.Quantity`
        The wavelength at which to evaluate the correction factor.

    temperature : `astropy.units.Quantity`
        The temperature of the plasma.

    Returns
    -------
    correction_factor

    References
    ----------
    `Mewe et al. Astron. and Astrophys Suppl. Ser., 62, 197-254 (1985) <https://ui.adsabs.harvard.edu/abs/1985A%26AS...62..197M/abstract>`_
    """
    return 27.83 * (temperature.to('MK').value + 0.65) ** (-1.33) + 0.15 * wavelength.to('angstrom').value ** (0.34) * temperature.to('MK').value ** (0.422)


@u.quantity_input(wavelength=u.nm, equivalencies=u.spectral())
def XraytubeEmission(wavelength, voltage: u.volt):
    """Implements the Kramer's law to estimate the continuum part of the
    spectrum of an x-ray tube.

    Parameters
    ----------
    wavelength : `astropy.units.Quantity`
        The wavelength at which to evaluate the correction factor.

    voltage : `astropy.units.Quantity`
        The voltage applied to the x-ray tube.

    Returns
    -------
    emission : The unscaled emission, provided without units
    """
    f = (wavelength / CutoffWavelength(voltage) - 1) * 1 / wavelength ** 2
    # remove unphysical negative values
    f[f < 0.0] = 0.0
    return f.to_value()


@u.quantity_input
def CutoffWavelength(voltage: u.volt):
    """
    Implements the Duane-Hunt law which provides the low energy (
    or high wavelength) cutoff for the emission of an x-ray tube.

    Parameters
    ----------
    voltage : `astropy.units.Quantity`
        The voltage applied to the x-ray tube.

    Returns
    -------
    cutoff : `astropy.units.Quantity`

    Reference
    ---------
    `Duane-Hunt Law <https://en.wikipedia.org/wiki/Duane–Hunt_law>`_
    """
    return h * c / (e * voltage).decompose()
