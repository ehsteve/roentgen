Data Sources
============

The data used in this package comes from the following data providers.

* The U.S National Institute of Standards and Technology (`NIST <https://www.nist.gov>`__)
* The Center for X-ray Optics and Advanced Light Source (`CXRO <http://cxro.lbl.gov/>`__) at the Lawrence Berkeley National Laboratory
* The Library for gamma and alpha emissions (`LARA <http://www.lnhb.fr/accueil/donnees-nucleaires/module-lara/>`_) at the Laboratoire Nationale Henri Becquerel

We would like to extend our thanks to these groups for making their data public.

The CXRO publish primarily makes their data available through the `X-ray Data Booklet <https://xdb.lbl.gov>`__.

The following list the data files contained in the data directory
along with their source and other relevant information to properly
interpret the data therein.

The elements and compounds_mixtures folders contains mass attenuation coefficients from the `NIST <https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients>`__.
More specifically, elements are from `here <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html>`__ while compounds and mixtures are from `here <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab4.html>`__.

elements.csv
------------
Source: `https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>`__
supplemented for elements z>92 by wikipedia (zovera), `NIST Atomic Spectra Database <https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html>`_ (ionization energy), and
`NIST ASTAR Material composition <https://physics.nist.gov/cgi-bin/Star/compos.pl?ap-text>`_ (mean excitation energy, density)
Provides translation between atomic number, element names, and element symbols.
Energy ranges from 1 keV to 20 MeV.

- **z**: Atomic number
- **name**: Element name
- **zovera**: Ratio of atomic number-to-mass (Z/A)
- **i**: mean excitation energy I in eV
- **density**: density in g/cm^3

compounds_mixtures.csv
----------------------
Source: `https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html>`__
Provides a list of all compounds and mixtures supported
density is in units of g/cm^3
Energy ranges from 1 keV to 20 MeV.

emission_lines.csv
------------------
Source: `CXRO X-ray data Booklet Table 1-3 <https://xdb.lbl.gov/Section1/Table_1-3.pdf>`__
Source: `AXAF Project Science: Widths for Characteristic X-Ray Lines <https://wwwastro.msfc.nasa.gov/xraycal/linewidths.html>`__

- **energy**: Energy of emission line in eV
- **transition**: Name of the transition in Siegbahn notation
- **intensity**: Relative intensity of the line
- **width**: Line width in eV

emission_energies.csv
---------------------
Source: `CXRO X-ray data Booklet Table 1-3 <https://xdb.lbl.gov/Section1/Table_1-2.pdf>`__
energy in eV

electron_binding_energies.csv
-----------------------------
Source: `CXRO X-ray data Booklet Table 1-1 <https://xdb.lbl.gov/Section1/Table_1-1.pdf>`_
energies in eV

lara
----
Source:  `Laboratoire Nationale Henri Becquerel Library for gamma and alpha emissions <http://www.lnhb.fr/accueil/donnees-nucleaires/module-lara/>`_
More specifically, data files are from `here <http://www.lnhb.fr/accueil/donnees-nucleaires/donnees-nucleaires-tableau/>`__