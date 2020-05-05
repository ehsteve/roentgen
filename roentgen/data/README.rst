Data Sources & References
=========================

The data used in this package comes from the following data providers.

* The U.S National Institute of Standards and Technology (`NIST <https://www.nist.gov>`_)
* The Center for X-ray Optics and Advanced Light Source (`CXRO <http://cxro.lbl.gov/>`_) at the Lawrence Berkeley National Laboratory

We would like to extend our thanks to these groups for making their data public.

The CXRO publish primarily make their data available through the `X-ray Data Booklet <https://xdb.lbl.gov`_.

The following list the data files contained in this data directory
along with their source and other relevant information to properly
interpret the data therein.

The elements and compounds_mixtures folders contains mass attenuation coefficients from the `NIST <https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients>`_.
More specifically, elements are from `here <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html>`_ while compounds and mixtures are from `here <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab4.html>`_.

elements.csv
------------
Source: `https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>`_
Provides translation between atomic number, element names, and element symbols.
* **z**: Atomic number
* **name**: Element name
* **zovera**: Ratio of atomic number-to-mass (Z/A)
* **i**: mean excitation energy I in eV
* **density**: density in g/cm^3

compounds_mixtures.csv
----------------------
Source: `https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html>`_
Provides a list of all compounds and mixtures supported
density is in units of g/cm^3

emission_lines.csv
------------------
Source: `CXRO X-ray data Booklet Table 1-3 <https://xdb.lbl.gov/Section1/Table_1-3.pdf>`_
* **energy**: Energy of emission line in eV
* **transition**: Name of the transition in Siegbahn notation
* **intensity**: Relative intensity of the line

emission_energies.csv
---------------------
Source: `CXRO X-ray data Booklet Table 1-3 <https://xdb.lbl.gov/Section1/Table_1-2.pdf>`_
energy in eV

electron_binding_energies.csv
-----------------------------
Source `CXRO X-ray data Booklet Table 1-1 <https://xdb.lbl.gov/Section1/Table_1-1.pdf>`_
energies in eV
