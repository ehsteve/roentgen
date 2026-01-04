.. _coding-standards:

****************
Coding Standards
****************

The purpose of the page is to describe the standards that are expected of all the code in this repository.
All developers should read and abide by the following standards.
Code which does not follow these standards closely will generally not be accepted.

We try to closely follow the coding style and conventions proposed by `Astropy <https://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_ and `Sunpy <https://docs.sunpy.org/en/latest/dev_guide/contents/code_standards.html>`_.

The following text is mostly cribbed from Sunpy.

Language Standard
=================

* All code must be compatible with Python 3.11 and later.

* The new Python 3 formatting style should be used (i.e.
  ``f"{spam:s}"`` instead of ``"%s" % "spam"``).

* The addition of new dependencies will be considered but is discouraged.

Coding Style/Conventions
========================

* The code will follow the standard `PEP8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_.
  In particular, this includes using only 4 spaces for indentation, and never tabs.

* **Follow the existing coding style** within a file and avoid making changes that are purely stylistic.
  Please try to maintain the style when adding or modifying code.

* Following PEP8's recommendation, absolute imports are to be used in general.
  We allow relative imports within a module to avoid circular import chains.

* The ``import numpy as np``, ``import matplotlib as mpl``, and ``import matplotlib.pyplot as plt`` naming conventions should be used wherever relevant.
  ``from packagename import *`` should never be used (except in ``__init__.py``)

* Classes should either use direct variable access, or Python's property mechanism for setting object instance variables.

* Classes should use the builtin `super` function when making calls to methods in their super-class(es) unless there are specific reasons not to.
  `super` should be used consistently in all subclasses since it does not work otherwise.

* Multiple inheritance should be avoided in general without good reason.

* ``__init__.py`` files for modules should not contain any significant implementation code. ``__init__.py`` can contain docstrings and code for organizing the module layout.

* We make use of Ruff for managing code formating and import sorting. Run `ruff check --fix` and `ruff format`.

Private code
============

It is often useful to designate code as private, which means it is not part of the user facing API, only used internally by HERMES, and can be modified without a deprecation period.
Any classes, functions, or variables that are private should either:

- Have an underscore as the first character of their name, e.g., ``_my_private_function``.
- If you want to do that to entire set of functions in a file, name the file with a underscore as the first character, e.g., ``_my_private_file.py``.

Documentation and Testing
=========================

* American English is the default language for all documentation strings and inline commands.
  Variables names should also be based on English words.

* Documentation strings must be present for all public classes/methods/functions, and must follow the form outlined in the :ref:`docs_guidelines` page.
  Additionally, examples or tutorials in the package documentation are strongly recommended.

* Write usage examples in the docstrings of all classes and functions whenever possible.
  These examples should be short and simple to reproduce–users should be able to copy them verbatim and run them.
  These examples should, whenever possible, be in the :ref:`doctests` format and will be executed as part of the test suite.

* Unit tests should be provided for as many public methods and functions as possible, and should adhere to the standards set in the :ref:`testing` document.

Data and Configuration
======================

* We store test data in ``./data/test`` as long as it is less than about 100 kB.

