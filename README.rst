========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - build status
      - |testing| |codestyle| |coverage|
    * - package
      - |version| |downloads| |wheel|

.. |docs| image:: https://readthedocs.org/projects/roentgen/badge/?version=latest
    :target: https://roentgen.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |testing| image:: https://github.com/ehsteve/roentgen/actions/workflows/testing.yml/badge.svg
    :target: https://github.com/ehsteve/roentgen/actions/workflows/testing.yml
    :alt: Build Status

.. |codestyle| image:: https://github.com/ehsteve/roentgen/actions/workflows/codestyle.yml/badge.svg
    :target: https://github.com/ehsteve/roentgen/actions/workflows/codestyle.yml
    :alt: Black linting

.. |coverage| image:: https://codecov.io/gh/ehsteve/roentgen/branch/master/graph/badge.svg?token=feNCnYTjB3
    :alt: Test coverage on codecov
    :target: https://codecov.io/gh/ehsteve/roentgen

.. |version| image:: https://img.shields.io/pypi/v/roentgen.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/roentgen

.. |downloads| image:: https://img.shields.io/pypi/dm/roentgen.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/roentgen

.. |wheel| image:: https://img.shields.io/pypi/wheel/roentgen.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/roentgen

.. end-badges

.. image:: docs/logo/roentgen_logo.svg
    :height: 100
    :width: 100

A Python package for the quantitative analysis of the interaction of energetic x-rays with matter.
This package is named after one of the discoverers of X-rays, `Wilhelm Röntgen <https://en.wikipedia.org/wiki/Wilhelm_Röntgen>`_.

Installation
============

::

    pip install roentgen

This project makes use of `Poetry <https://python-poetry.org>`_ for dependency management. To install this project for development, clone the repository and then run the following command inside the package directory

::

    poetry install --with dev,docs,gui


Documentation
=============

http://roentgen.readthedocs.io/en/stable/

GUI
===
This package provides a gui interface to quickly investigate the absorption and transmission of x-rays through different materials.
It is based on `bokeh <https://docs.bokeh.org/en/stable/>`_. To run it locally use the following command

::

   bokeh serve --show <roengten_directory>/gui


Data Sources
============
This package includes on a number of data files which were translated and imported from a few key sources.
The package developers would like to thank the following data providers

* The U.S National Institute of Standards and Technology (NIST)
* The Center for X-ray Optics and Advanced Light Source at the Lawrence Berkeley National Laboratory

For more information see the `README <roentgen/data/README.rst>`_ in data directory.

Contributing
============

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.
Have a look at the `great guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__ from the `sunpy project <https://sunpy.org>`__ which provides advice for new contributors.

Code of Conduct
===============

When you are interacting with members of this community, you are asked to follow the SunPy `Code of Conduct <https://sunpy.org/coc>`__.
