# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'coolpackage'
copyright = '2026, G. Weston, M. T. Degiacomi'
author = 'G. Weston, M. T. Degiacomi'

# The full version, including alpha/beta/rc tags
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    "sphinx.ext.autosummary",
    'sphinx.ext.githubpages',
    "autoapi.extension",
]

autoapi_type = "python"
autoapi_dirs = [os.path.abspath('../../features')]

# Use our own AutoAPI templates, so that the title and the introductory text of
# the generated landing page can be edited (see _templates/autoapi/index.rst).
# Any template not overridden here falls back to the one shipped with autoapi.
autoapi_template_dir = "_templates/autoapi"


def autoapi_prepare_jinja_env(jinja_env):
    """Make a ``first_sentence`` filter available to the AutoAPI templates.

    AutoAPI's own ``summary`` attribute returns the first *line* of a docstring,
    which cuts a multi-line sentence mid-way. This filter instead returns the
    first complete sentence of the first paragraph, collapsed onto one line.
    """

    def first_sentence(docstring):
        paragraph = []
        for line in (docstring or "").strip().splitlines():
            if not line.strip():
                break
            paragraph.append(line.strip())
        text = " ".join(" ".join(paragraph).split())
        head, separator, _ = text.partition(". ")
        return head + "." if separator else text

    jinja_env.filters["first_sentence"] = first_sentence


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinxdoc'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'coolpackagedoc'

# -- Options for Texinfo output -------------------------------------------
add_module_names = False
autoclass_content = "both"

# Modeller is a licensed package that is not necessarily installed on the
# machine building the documentation. Mocking it lets autodoc import
# protein.py, which reaches Modeller through patcher.py.
autodoc_mock_imports = ["modeller"]
