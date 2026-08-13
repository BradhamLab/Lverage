"""Configuration for the Lverage documentation."""

from importlib.metadata import version as distribution_version


project = "Lverage"
author = "Bradham Lab"
copyright = "2026, Bradham Lab"
release = distribution_version("lverage")
version = release

extensions = [
    "numpydoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

autosummary_generate = True
autodoc_typehints = "description"
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

numpydoc_show_class_members = False
numpydoc_validation_checks = {
    "GL06",
    "GL07",
    "GL08",
    "PR01",
    "PR02",
    "PR03",
    "PR04",
    "PR07",
    "RT01",
    "RT02",
    "RT03",
    "SS01",
    "YD01",
}
numpydoc_validation_exclude = {
    r"lverage\.pipeline\.LverageCode$",
}

nitpicky = True

html_theme = "pydata_sphinx_theme"
