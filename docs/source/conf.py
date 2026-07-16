# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import subprocess
import sys

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "libhankel"
copyright = "2026, Letizia"
author = "Letizia"
release = "v0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, "../..")  # adjust path

extensions = [
    "breathe",
    "sphinx.ext.autodoc",
]

breathe_projects = {"libhankel": "../xml"}

breathe_default_project = "libhankel"


templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"

autodoc_docstring_signature = True
breathe_show_include = False


# Hook doxygen
def run_doxygen(app):
    try:
        subprocess.check_call(["doxygen"])
    except Exception as e:
        print("Doxygen execution failed:", e)


def setup(app):
    app.connect("builder-inited", run_doxygen)


autodoc_default_options = {
    "members": True,
    "imported-members": True,
}
