# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath(".."))

now = datetime.now()

# -- Project information -----------------------------------------------------

project = "roentgen"
copyright = f"US Government (not copyrighted) 2023-{now.year}"
author = "Steven D. Christe, Nabil Freij, Shane Maloney, Daniel Ryan"

# The full version, including alpha/beta/rc tags
from roentgen import __version__

release = __version__
is_development = ".dev" in __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "matplotlib.sphinxext.plot_directive",
    "pytest_doctestplus.sphinx.doctestplus",
    "sphinx_copybutton",
]

# to exclude traditional Python prompts from your copied code
copybutton_prompt_text = ">>> "

# plot_directive default to always show source when including a plot
plot_include_source = True

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The reST default role (used for this markup: `text`) to use for all
# documents. Set to the "smart" one.
default_role = "obj"

autosummary_generate = True  # Turn on sphinx.ext.autosummary
autosummary_ignore_module_all = False
autosummary_imported_members = False

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": True,
}

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://data.astropy.org/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://docs.scipy.org/doc/numpy/",
        (None, "http://data.astropy.org/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://data.astropy.org/intersphinx/scipy.inv"),
    ),
    "matplotlib": (
        "https://matplotlib.org/",
        (None, "http://data.astropy.org/intersphinx/matplotlib.inv"),
    ),
    "astropy": ("http://docs.astropy.org/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_logo = "logo/roentgen.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
html_theme = "pydata_sphinx_theme"

html_theme_options = {
   "use_edit_page_button": True,
   "back_to_top_button": True,
}

html_context = {
    "display_github": True,
    "github_user": "ehsteve",
    "github_repo": "roentgen",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]
