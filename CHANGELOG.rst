Changelog
=========

We follow `semantic versioning <https://semver.org>`_ practices so versions increment as

#. MAJOR version when you make incompatible API changes
#. MINOR version when you add functionality in a backward compatible manner
#. PATCH version when you make backward compatible bug fixes

Note that PATCH versions releases may not be documented below.


2.4.0 (2026-Jan)
----------------
* Moved to ruff for linting
* Moved to uv
* Added nuclides module that provides access to radioactive nuclide emission lines


2.3.0 (2024-10-14)
------------------
* Added line widths to the emission line list


2.2.0 (2024-03-14)
------------------
* Fix and improvement to gui
* General improvements to documentation and examples
* Removed poetry
* added pre-commit support and updated testing


2.1.0 (2022-01-07)
------------------
* Improvement to documentation to illustrate custom materials
* Added tutorial and documentation to illustrate interpolation
* Going outside of data range (1 keV to 2 MeV) now provides a ValueError instead of failing quietly and filling in zeros for extrapolated values.


2.0.1 (2022-10-10)
------------------
* The Material class can now represent substances with multiple and arbitrary constitutients
* Moved to poetry for dependency management and packaging
* Automated testing now also using poetry
* Added data source information into the documentation
* Updated dependencies versions
* Added a Code of Conduct
* If material (element or compounds) are not found, raises a ValueError instead of returning None
* Increased test coverage
* Added attenuator section to guide.

1.0.0 (2020-04-23)
------------------
* First stable release
* Added gui, deployed to heroku
* Added data and support for emission lines
* Added Response class
* Major update to documentation and moved to sunpy doc theme

0.1.0 (2017-03-24)
------------------
* First release