.. _dev_env:

*********************
Developer Environment
*********************

This package uses `uv <https://uv.readthedocs.io/en/latest/>`_ to manage the development environment.
Uv is a tool to manage virtual environments and dependencies for Python projects.
The following instructions assume you have `uv` installed.

After you've cloned the repository, navigate to the root of the repository and run::

    $ uv sync

This will create a virtual environment and install all the dependencies needed for development, including testing and building the documentation.

You can then activate the virtual environment with::

    $ uv activate

To deactivate the virtual environment, simply run::

    $ uv deactivate

Once the virtual environment is activated, you can run various commands within it using `uv run <command>`.

For example, to run the test suite, you would use::

    $ uv run pytest

To build the documentation, you would use::

    $ uv run sphinx-build -b html docs/ docs/_build/html
