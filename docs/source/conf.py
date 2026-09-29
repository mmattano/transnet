# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Document the package from the source tree, without requiring an install.
sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------

project = "TransNet"
copyright = "2026, Matthias Anagho-Mattanovich"
author = "Matthias Anagho-Mattanovich"

try:
    from transnet import __version__ as release
except Exception:                                   # pragma: no cover
    release = "0.2.0"

version = release


# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",       # NumPy-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
# The walkthroughs are jupytext scripts, not .ipynb: the .py file is the
# source of truth, so there is no notebook output to keep in step with it.
# The offline ones execute at build time (seconds, no network); the studies
# need built networks and minutes, so they render from saved output.
sys.path.insert(0, os.path.abspath("."))     # for _readers below
nb_custom_formats = {".py": ["_readers.read", {}]}
nb_execution_mode = "cache"
nb_execution_timeout = 300
# The studies need built networks and minutes of compute, so they are executed
# once by `make studies` into docs/source/studies/*.ipynb, outputs included,
# and rendered from those saved outputs rather than re-run here.
nb_execution_excludepatterns = [
    "studies/*",                    # executed by `make studies`, rendered as saved
    "notebooks/studies/*",
    "notebooks/extra/*",
    "notebooks/walkthroughs/external_annotation.py",   # needs the network
]
myst_enable_extensions = ["colon_fence", "dollarmath"]
# Interactive figures are the point of the interactive view, so prefer the
# plotly output over its static fallback when rendering to HTML.
nb_mime_priority_overrides = [
    ("html", "application/vnd.plotly.v1+json", 10),
    ("html", "text/html", 20),
    ("html", "image/png", 30),
]

exclude_patterns = [
    "_build", "**.ipynb_checkpoints",
    "_readers.py", "conf.py",
    # The studies need built networks and minutes of compute; their narrative
    # lives in the .rst pages beside them, with figures from a real run.
    "notebooks/studies/*", "notebooks/extra/*", "notebooks/README.md",
    "notebooks/exports/*",
]

# Sphinx should not import heavy optional database clients just to read a
# docstring; the package imports fine without them, but the docs build should
# not depend on them at all.
autodoc_mock_imports = [
    "bioservices", "mygene", "pyensembl", "zeep", "Bio",
]

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}


# -- Options for HTML output -------------------------------------------------

html_theme = "alabaster"
html_static_path = ["_static"]
html_title = "TransNet: trans-omics network analysis"
