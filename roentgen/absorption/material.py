""" """

import numpy as np
from astropy.units import Quantity
from scipy import interpolate

import astropy.units as u

import roentgen
from roentgen.absorption.cross_section_database import get_database
from roentgen.util import (
    get_atomic_number,
    get_compound_index,
    get_material_density,
    get_material_name,
    get_material_symbol,
    is_an_element,
    is_in_known_compounds,
)

__all__ = ["Material", "MassAttenuationCoefficient", "Stack", "Response"]

_package_directory = roentgen._package_directory
_data_directory = roentgen._data_directory


class Material(object):
    """
    An object which enables the calculation of the x-ray transmission and
    absorption of a material (e.g. an element or a compound/mixture).

    A material may be composed of a single atomic element such as Aluminum ('Al'), or composed of number of elements and/or compounds.

    Parameters
    ----------
    material_str : str or dict
        A string representation of the material which includes an element symbol
        (e.g. Si), an element name (e.g. Silicon), or the name of a compound
        (e.g. cdte, mylar). For supported elements see :download:`elements.csv <../../roentgen/data/elements.csv>` and for compounds see :download:`compounds_mixtures.csv <../../roentgen/data/compounds_mixtures.csv>`.
        Can also be a dictionary of element and compounds with fractional masses (ex. {"Cu":0.70, "Zn":0.30})
    thickness : `astropy.units.Quantity`
        The thickness of the material
    density : `astropy.units.Quantity`, optional
        The density of the material.
        If not provided, uses default values which can be found in :download:`elements.csv <../../roentgen/data/elements.csv>` for elements or
        in :download:`compounds_mixtures.csv <../../roentgen/data/compounds_mixtures.csv>` for compounds.
        If many materials are present, calculates the weighted density.
    database : `str`, optional
        The cross-section data source to use: ``"nist"`` (default),
        ``"mcmaster"``, ``"xcom"``, ``"xraylib"`` or ``"xraydb"``. Despite the name, ``"nist"`` isn't
        "the" NIST source -- ``"xcom"`` is *also* NIST data, just from a
        different NIST publication (see
        `~roentgen.absorption.cross_section_database` for the distinction).
        All but ``"nist"`` additionally provide the
        photoelectric/Compton/Rayleigh split (see
        `photoelectric_mass_attenuation_coefficient`,
        `compton_mass_attenuation_coefficient`) but only for materials made
        up of individual elements (not named compounds). ``"mcmaster"`` is
        the McMaster et al. 1969/1970 data (valid 1-1000 keV, matches SSW's
        ``resp_calc.pro``); ``"xcom"`` is current NIST XCOM data (valid
        1 keV-1e8 keV, one more element covered) -- see
        `roentgen.absorption.mcmaster`/`roentgen.absorption.xcom`. ``"xraylib"`` and
        ``"xraydb"`` use the optional packages of those names.

    .. warning::
        Elements beyond z = 92 are not supported by this class.

    Attributes
    ----------
    symbols : `list`
        A list of material symbol
    material_names : `list`
        A list of material names
    name : `str`
        A name for the material
    fractional_masses : `np.ndarray`
        A normalized array of fractional masses
    mass_attenuation_coefficients : `list`
        A list of `MassAttenuationCoefficient`

    Examples
    --------
    >>> from roentgen.absorption.material import Material
    >>> import astropy.units as u
    >>> detector = Material('cdte', 500 * u.um)
    >>> thermal_blankets = Material('mylar', 0.5 * u.mm)
    >>> bronze = Material({"Cu": 0.88, "Sn": 0.12}, 1 * u.mm)
    """

    @u.quantity_input
    def __init__(self, material_input, thickness: u.m, density=None, database="nist"):
        self.thickness = thickness
        self._database = get_database(database)
        if isinstance(density, u.Quantity):
            self.density = density
        if isinstance(material_input, str):
            self.list_names = [get_material_name(material_input)]
            self.list_symbols = [get_material_symbol(material_input)]
            self.mass_attenuation_coefficients = [MassAttenuationCoefficient(material_input)]
            self.symbol = self.mass_attenuation_coefficients[0].symbol
            self.name = self.mass_attenuation_coefficients[0].name
            self.fractional_masses = np.ones(1)
            if density is None:
                self.density = get_material_density(material_input)
        elif isinstance(material_input, dict):
            self.list_names = [
                get_material_name(this_str) for this_str in list(material_input.keys())
            ]
            self.list_symbols = [
                get_material_symbol(this_str) for this_str in list(material_input.keys())
            ]
            # normalize the fractional masses
            fractional_masses = np.array(list(material_input.values()))
            self.fractional_masses = fractional_masses / fractional_masses.sum()
            self.name = "".join(f"{this_name}" for this_name in self.list_names)
            self.symbol = "".join(f"{this_symbol}" for this_symbol in self.list_symbols)
            self.mass_attenuation_coefficients = [
                MassAttenuationCoefficient(e) for e in self.list_names
            ]
            if density is None:
                # calculate the average weighted density
                densities = (
                    [
                        get_material_density(this_material).to("kg/m**3").value
                        for this_material in self.list_names
                    ]
                    * u.kg
                    / u.m**3
                )
                self.density = np.average(densities, weights=self.fractional_masses)
        else:
            raise TypeError("Material input must be a string or a dictionary.")

    def __repr__(self):
        """Returns a developer-relevant representation."""
        # at this point, no reason for this to be different than __str__
        return self.__str__()

    def __str__(self):
        """Returns a human-readable user-focused representation."""
        txt = f"Material('{self.name}' thickness={self.thickness} density={self.density.to('kg/m**3'):2.1f})"
        return txt

    def __add__(self, other):
        if isinstance(other, Material):
            return Stack([self, other])
        elif isinstance(other, Stack):
            return Stack([self] + other.materials)
        else:
            raise TypeError(f"Cannot add {self} and {other}")

    @u.quantity_input(energy=u.keV)
    def mass_attenuation_coefficient(self, energy):
        result = np.sum(
            np.vstack(
                [
                    self._database.total(symbol, energy) * frac_mass
                    for symbol, frac_mass in zip(self.list_symbols, self.fractional_masses)
                ]
            ),
            axis=0,
        )
        if energy.isscalar:
            return result[0]
        else:
            return result

    @u.quantity_input(energy=u.keV)
    def photoelectric_mass_attenuation_coefficient(self, energy):
        """The photoelectric-only mass attenuation coefficient.

        Only available for materials made up of individual elements, and
        only for databases which provide component cross-sections
        (``database="mcmaster"``, valid 1-1000 keV, or ``database="xcom"``,
        valid 1 keV-1e8 keV).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        result = np.sum(
            np.vstack(
                [
                    self._database.photoelectric(symbol, energy) * frac_mass
                    for symbol, frac_mass in zip(self.list_symbols, self.fractional_masses)
                ]
            ),
            axis=0,
        )
        if energy.isscalar:
            return result[0]
        else:
            return result

    @u.quantity_input(energy=u.keV)
    def compton_mass_attenuation_coefficient(self, energy):
        """The incoherent (Compton) scattering-only mass attenuation coefficient.

        Only available for materials made up of individual elements, and
        only for databases which provide component cross-sections
        (``database="mcmaster"``, valid 1-1000 keV, or ``database="xcom"``,
        valid 1 keV-1e8 keV).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV
        """
        result = np.sum(
            np.vstack(
                [
                    self._database.compton(symbol, energy) * frac_mass
                    for symbol, frac_mass in zip(self.list_symbols, self.fractional_masses)
                ]
            ),
            axis=0,
        )
        if energy.isscalar:
            return result[0]
        else:
            return result

    @u.quantity_input(energy=u.keV)
    def transmission(self, energy):
        """Provide the transmission fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        coefficients = self.mass_attenuation_coefficient(energy)
        transmission = np.exp(-coefficients * self.density * self.thickness)
        return transmission.value  # remove the dimensionless unit

    @u.quantity_input(energy=u.keV)
    def absorption(self, energy):
        """Provides the absorption fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        return 1.0 - self.transmission(energy)

    def linear_attenuation_coefficient(self, energy: u.keV):
        """Provides the linear attenuation coefficient as a function of energy.

        linear coeff = mass coeff * density.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        return self.mass_attenuation_coefficient(energy) * self.density

    def photoelectric_linear_attenuation_coefficient(self, energy: u.keV):
        """Provides the photoelectric-only linear attenuation coefficient as a
        function of energy. See `photoelectric_mass_attenuation_coefficient`.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return self.photoelectric_mass_attenuation_coefficient(energy) * self.density

    def compton_linear_attenuation_coefficient(self, energy: u.keV):
        """Provides the incoherent (Compton) scattering-only linear attenuation
        coefficient as a function of energy. See `compton_mass_attenuation_coefficient`.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.
        """
        return self.compton_mass_attenuation_coefficient(energy) * self.density


class Stack(object):
    """
    An object which enables the calculation of the x-ray transmission and
    absorption of a stack of materaials.
    This object is created automatically when `Material` objects are added together.

    Parameters
    ----------
    materials : list
        A list of `Material` objects

    Examples
    --------
    >>> from roentgen.absorption.material import Material, Stack
    >>> import astropy.units as u
    >>> detector = Stack([Material('Pt', 5 * u.um), Material('cdte', 500 * u.um)])
    >>> optical_path = Material('mylar', 50 * u.micron) + Material('Al', 1 * u.mm)
    """

    def __init__(self, materials):
        self.materials = materials

    def __add__(self, other):
        if isinstance(other, Material):
            return Stack(self.materials + [other])
        elif isinstance(other, Stack):
            return Stack(self.materials + other.materials)
        else:
            raise TypeError(f"Cannot add {self} and {other}")

    def __repr__(self) -> str:
        """Returns a developer-relevant representation."""
        # at this point, no reason for this to be different than __str__
        return self.__str__()

    def __str__(self) -> str:
        """Returns a human-readable user-focused representation."""
        txt = "Stack(["
        for this_material in self.materials:
            txt += f"{this_material}, "
        txt = f"{txt[:-2]}])"
        return txt

    def transmission(self, energy):
        """Provides the transmission fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        transmission = np.ones(len(energy), dtype=float)
        for material in self.materials:
            this_transmission = material.transmission(energy)
            transmission *= this_transmission
        return transmission

    def absorption(self, energy):
        """Provides the absorption fraction (0 to 1).

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        return 1.0 - self.transmission(energy)


class Response(object):
    """
    An object to handle the response of a detector material which includes
    an optical path or filter through which x-rays must first traverse before
    reaching the detector.

    Parameters
    ----------
    optical_path : Stack or Material
        A list of Material objects which make up the optical path.

    detector : Material or None
        A Material which represents the detector material where the x-rays
        are absorbed. If provided with None, than assume a perfectly absorbing
        detector material.

    Examples
    --------
    >>> from roentgen.absorption.material import Material, Response, Stack
    >>> import astropy.units as u
    >>> optical_path = Stack([Material('air', 1 * u.m), Material('Al', 500 * u.mm)])
    >>> resp = Response(optical_path, detector=Material('cdte', 500 * u.um))
    """

    def __init__(self, optical_path, detector):
        # make sure the materials are a list since we iterate over them
        # to calculate the transmission
        if isinstance(optical_path, Stack) or isinstance(optical_path, Material):
            self.optical_path = optical_path
        else:
            raise TypeError("optical_path must be a Stack or Material")

        if isinstance(detector, Material):
            self.detector = detector
        else:
            raise TypeError("detector must be a Material")

    def __repr__(self) -> str:
        """Returns a developer-relevant representation."""
        # at this point, no reason for this to be different than __str__
        return self.__str__()

    def __str__(self) -> str:
        """Returns a human-readable user-focused representation."""
        txt = f"Response(optical_path={self.optical_path} detector={self.detector})"
        return txt

    def response(self, energy):
        """Returns the response as a function of energy which corresponds to the
        transmission through the optical path multiplied by the absorption in
        the detector.

        Parameters
        ----------
        energy : `astropy.units.Quantity`
            An array of energies in keV.

        Raises
        ------
        ValueError
            If energy is outside of the interpolation range of 1 keV to 20 MeV.
        """
        # calculate the transmission
        transmission = np.ones(len(energy), dtype=float)
        detector_absorption = np.ones(len(energy), dtype=float)

        transmission = self.optical_path.transmission(energy)
        detector_absorption = self.detector.absorption(energy)

        return transmission * detector_absorption

    @u.quantity_input
    def response_matrix(
        self,
        edges_in: Quantity[u.keV],
        edges_out: Quantity[u.keV] | None = None,
        fwhm=None,
        theta: Quantity[u.deg] = 0 * u.deg,
        area: Quantity[u.cm**2] = 1 * u.cm**2,
        scattering_function: bool = False,
        ssw_compatible: bool = False,
    ):
        """Compute the full detector energy-redistribution (response) matrix.

        Port of SSW's ``resp_calc.pro``: for each input photon energy bin,
        the probability of a photon being measured in each output (PHA)
        energy bin, including the photopeak, single Compton-scatter escape,
        K-shell fluorescence escape peaks, and (if ``fwhm`` is given)
        Gaussian resolution broadening.

        ``self.detector`` must be a `Material` built with
        ``database="mcmaster"`` or ``database="xcom"`` (a single element, or
        a dict of elements) -- the photoelectric/Compton split this needs
        isn't available from the default NIST-backed data. ``"mcmaster"``
        matches ``resp_calc.pro`` exactly (valid 1-1000 keV); ``"xcom"`` uses
        current NIST data over a wider range (1 keV-1e8 keV). Fluorescence
        escape uses the constants from ``resp_calc.pro`` for Ge, Si, NaI, CsI, CdTe and SXS
        (`roentgen.response.fluorescence.SUPPORTED_DETECTORS`) with ``database="mcmaster"``, and
        otherwise a generic model using the database's own K-shell data.

        Parameters
        ----------
        edges_in : `astropy.units.Quantity`
            Input photon energy bin edges, ``n + 1`` values for ``n`` bins.
        edges_out : `astropy.units.Quantity`, optional
            Output (PHA channel) bin edges. Defaults to ``edges_in``.
        fwhm : callable, optional
            ``fwhm(energy_kev) -> fwhm_kev``, the detector's energy
            resolution as a function of energy (both plain arrays in keV).
            If not given, no resolution broadening is applied
            (``pulse_height_matrix`` and ``smatrix`` will be `None`).
        theta : `astropy.units.Quantity`, optional
            Angle of incidence from the detector normal. Default 0 deg.
        area : `astropy.units.Quantity`, optional
            Detector geometric area. Default 1 cm2.
        scattering_function : `bool`, optional
            If `True`, include the incoherent scattering function S(q,Z) in the
            shape of the Compton continuum (see `roentgen.response.compton`),
            otherwise plain Klein-Nishina as in ``resp_calc.pro``. Default `False`.
        ssw_compatible : `bool`, optional
            If `True`, reproduce ``resp_calc.pro`` exactly, including three known defects:
            the Compton continuum uses ``klein.pro``'s cross-section, which has an extra factor
            of E'/E (see `roentgen.response.compton.klein`); the last Compton output bin is
            oversized when the maximum electron energy is above the last channel edge; and the
            window/filter transmission is not applied to the Compton continuum. Default `False`
            corrects all three.

        Returns
        -------
        `~roentgen.response.response_matrix.ResponseMatrix`
        """
        # bin edges come straight from the user and the port silently gives wrong
        # numbers if they are not ascending
        for edges in (edges_in, edges_out):
            if edges is not None and np.any(np.diff(edges.value) <= 0):
                raise ValueError("energy bin edges must be strictly increasing")

        from roentgen.response.response_matrix import _compute_response_matrix

        return _compute_response_matrix(
            self,
            len(edges_in) - 1,
            edges_in[0],
            edges_in[-1],
            fwhm=fwhm,
            edges_in=edges_in,
            edges_out=edges_out,
            theta=theta,
            area=area,
            nosigma=fwhm is None,
            scattering_function=scattering_function,
            ssw_compatible=ssw_compatible,
        )


class MassAttenuationCoefficient(object):
    """
    The mass attenuation coefficient.

    Parameters
    ----------
    material_str : str
        A string representation of the material which includes an element symbol
        (e.g. Si), an element name (e.g. Silicon), or the name of a compound
        (e.g. cdte, mylar).

    Attributes
    ----------
    data : `astropy.units.Quantity` array
        The mass attenuation data values.
    energy : `astropy.units.Quantity`
        The energy values of the mass attenuation values.
    symbol : `str`
        The material symbol
    name : `str`
        The material name
    func : `lambda func`
        A function which returns the interpolated mass attenuation value at
        any given energy. Energies must be given by an `astropy.units.Quantity`.
        The interpolation range is 1 keV to 20 MeV.
        Going outside that range will result in a ValueError.

    Examples
    --------
    >>> from roentgen.absorption.material import MassAttenuationCoefficient
    >>> mass_atten = MassAttenuationCoefficient('air')
    """

    def __init__(self, material):
        """
        Parameters
        ----------
        material : str
            A string representation of the material which includes an element symbol
            (e.g. Si), an element name (e.g. Silicon), or the name of a compound
            (e.g. cdte, mylar).
        """
        if is_an_element(material):
            atomic_number = get_atomic_number(material)
            filename = "z" + str(atomic_number).zfill(2) + ".csv"
            datafile_path = _data_directory / "elements" / filename
            symbol = roentgen.elements[atomic_number - 1]["symbol"]
            name = roentgen.elements[atomic_number - 1]["name"]
        elif is_in_known_compounds(material):
            compound_index = get_compound_index(material)
            symbol = roentgen.compounds[compound_index]["symbol"]
            name = roentgen.compounds[compound_index]["name"]
            filename = symbol.replace(" ", "_") + ".csv"
            datafile_path = _data_directory / "compounds_mixtures" / filename
        else:
            raise ValueError(f"Element or compound {material} not found.")
        data = np.loadtxt(datafile_path, delimiter=",")
        # find the material in our list
        self.symbol = symbol
        self.name = name
        self.energy = u.Quantity(data[:, 0] * 1000, "keV")
        self.data = u.Quantity(data[:, 1], "cm^2/g")

        self._remove_double_vals_from_data()

        data_energy_kev = np.log10(self.energy.value)
        data_attenuation_coeff = np.log10(self.data.value)
        self._f = interpolate.interp1d(
            data_energy_kev,
            data_attenuation_coeff,
            bounds_error=True,
            assume_sorted=True,
        )
        self.func = lambda x: u.Quantity(10 ** self._f(np.log10(x.to("keV").value)), "cm^2/g")

    def __repr__(self) -> str:
        """Returns a developer-relevant representation."""
        # at this point, no reason for this to be different than __str__
        return self.__str__()

    def __str__(self) -> str:
        """Returns a human-readable user-focused representation."""
        txt = f"MassAttenuationCoefficient('{self.name}')"
        return txt

    def _remove_double_vals_from_data(self):
        """Remove double-values energy values. Edges are represented with
        the same energy index and at the bottom and top value of the edge. This
        must be removed to enable correct interpolation."""
        uniq, count = np.unique(self.energy, return_counts=True)
        duplicates = uniq[count > 1]
        for this_dup in duplicates:
            ind = (self.energy == this_dup).nonzero()
            # shift the first instance of the energy, the bottom of the edge
            self.energy[ind[0][0]] -= 1e-3 * u.eV
