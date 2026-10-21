"""
Pluggable cross-section data sources for `~roentgen.absorption.material.Material`.

roentgen ships three data sources. Note that two of them are both NIST data, just from different NIST publications,
 ``"nist"`` is *not* "the" NIST source, it's specifically the older total-only tables:

- ``"nist"`` (the default) — roentgen's existing per-element/compound total
  mass attenuation coefficient tables (see `MassAttenuationCoefficient`),
  from NIST's *X-Ray Mass Attenuation Coefficients* tables (Hubbell &
  Seltzer, physics.nist.gov/PhysRefData/XrayMassCoef/). Only the total
  coefficient is available; there's no split between interaction types in
  this data.
- ``"mcmaster"`` — the McMaster et al. 1969/1970 coefficients used by
  SSW's xray package (see `roentgen.absorption.mcmaster`). Provides
  the total *and* the individual photoelectric, coherent (Rayleigh), and
  incoherent (Compton) components, needed to build a detector response
  matrix (`Response.response_matrix`). Valid 1-1000 keV.
- ``"xcom"`` — current NIST *XCOM* data (physics.nist.gov/PhysRefData/Xcom/,
  a separate, actively-maintained NIST database/tool -- distinct from the
  ``"nist"`` source above despite both being NIST), also bundled with SSW
  (see `roentgen.absorption.xcom`), read via the same component split as
  ``"mcmaster"`` but over a much wider energy range (1 keV to 1e8 keV).
  Covers one more element than ``"mcmaster"`` (z=1-92 except z=91, vs. 87
  elements for ``"mcmaster"``).

Which source a `Material` uses is an explicit, per-`Material` choice (the
``database`` keyword), rather than something a function silently decides on
your behalf. This also gives a slot for adding further sources later
without another API change.
"""

from roentgen.absorption import mcmaster, xcom

__all__ = [
    "CrossSectionDatabase",
    "NISTDatabase",
    "McMasterDatabase",
    "XCOMDatabase",
    "get_database",
]


class CrossSectionDatabase:
    """Interface for a source of per-element mass attenuation coefficients.

    Every method below shares the same signature and units:

    Parameters
    ----------
    element : `str` or `int`
        An element symbol, name, or atomic number.
    energy : `astropy.units.Quantity`
        Energy at which to evaluate the coefficient. Valid range depends on
        the concrete database.

    Returns
    -------
    `astropy.units.Quantity`
        Mass attenuation coefficient, in cm2/g.
    """

    # short name used to select this database, e.g. `Material(..., database="nist")`
    name = None

    def total(self, element, energy):
        """Total mass attenuation coefficient."""
        raise NotImplementedError

    def photoelectric(self, element, energy):
        """Photoelectric-only mass attenuation coefficient, in cm2/g."""
        raise NotImplementedError(
            f"The '{self.name}' database does not provide component "
            "(photoelectric/compton/rayleigh) cross-sections, only the total. "
            "Construct the Material with database='mcmaster' or database='xcom' for component access."
        )

    def compton(self, element, energy):
        """Incoherent (Compton) scattering mass attenuation coefficient, in cm2/g."""
        raise NotImplementedError(
            f"The '{self.name}' database does not provide component "
            "(photoelectric/compton/rayleigh) cross-sections, only the total. "
            "Construct the Material with database='mcmaster' or database='xcom' for component access."
        )

    def rayleigh(self, element, energy):
        """Coherent (Rayleigh) scattering mass attenuation coefficient, in cm2/g."""
        raise NotImplementedError(
            f"The '{self.name}' database does not provide component "
            "(photoelectric/compton/rayleigh) cross-sections, only the total. "
            "Construct the Material with database='mcmaster' or database='xcom' for component access."
        )


class NISTDatabase(CrossSectionDatabase):
    """roentgen's default per-element/compound total mass attenuation data,
    from NIST's *X-Ray Mass Attenuation Coefficients* tables (Hubbell &
    Seltzer) -- not to be confused with `XCOMDatabase`, a *different* NIST
    source (NIST XCOM) that this module also provides."""

    name = "nist"

    def total(self, element, energy):
        # deferred import: material.py imports this module for `get_database`,
        # so importing MassAttenuationCoefficient at module scope would be circular.
        from roentgen.absorption.material import MassAttenuationCoefficient

        # Performace regression this re-reads the CSV and rebuilds the interpolator on every
        # call,  extra overhead compared to a bare MassAttenuationCoefficient list Material used
        # before it went through the database. Could fix by caching one per element, e.g.
        # functools.lru_cache on a small helper if needs be
        return MassAttenuationCoefficient(element).func(energy)


class McMasterDatabase(CrossSectionDatabase):
    """McMaster et al. 1969/1970 coefficients (matches SSW xray's xsec.pro/det_xsec.pro)."""

    name = "mcmaster"

    def total(self, element, energy):
        return mcmaster.total_mass_attenuation(element, energy)

    def photoelectric(self, element, energy):
        return mcmaster.photoelectric_mass_attenuation(element, energy)

    def compton(self, element, energy):
        return mcmaster.compton_mass_attenuation(element, energy)

    def rayleigh(self, element, energy):
        return mcmaster.rayleigh_mass_attenuation(element, energy)


class XCOMDatabase(CrossSectionDatabase):
    """NIST XCOM data, bundled with SSW (matches ``xsec.pro``'s
    ``use_xcom`` path)."""

    name = "xcom"

    def total(self, element, energy):
        return xcom.total_mass_attenuation(element, energy)

    def photoelectric(self, element, energy):
        return xcom.photoelectric_mass_attenuation(element, energy)

    def compton(self, element, energy):
        return xcom.incoherent_mass_attenuation(element, energy)

    def rayleigh(self, element, energy):
        return xcom.coherent_mass_attenuation(element, energy)


_DATABASES = {"nist": NISTDatabase(), "mcmaster": McMasterDatabase(), "xcom": XCOMDatabase()}


def get_database(database):
    """Resolve a database name (or an existing `CrossSectionDatabase` instance) to an instance."""
    if isinstance(database, CrossSectionDatabase):
        return database
    try:
        return _DATABASES[database]
    except KeyError:
        raise ValueError(
            f"Unknown cross-section database {database!r}. Available: {list(_DATABASES)}."
        ) from None
