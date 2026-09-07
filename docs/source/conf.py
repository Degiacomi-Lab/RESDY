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

# Root of the repository. Needed because the modules import each other through the
# "src" namespace package (e.g. "from src.helper import ShutUp").
REPO_ROOT = os.path.abspath('../..')

# The package source itself. Needed because the autodoc directives in the .rst files
# refer to the modules by their bare name (e.g. "protein.PDB" rather than
# "src.protein.PDB").
SRC_DIR = os.path.join(REPO_ROOT, 'src')

for path in (REPO_ROOT, SRC_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

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
autoapi_dirs = [os.path.join(SRC_DIR, 'features')]

# Internal machinery of the features package: "feature.py" is the template users copy when
# adding a feature of their own, and "error_reporting.py" is the shared error log writer.
# Neither is a feature the user can request, so neither belongs on the feature pages.
autoapi_ignore = ['*/feature.py', '*/error_reporting.py']

# Use our own AutoAPI templates, so that the title and the introductory text of
# the generated landing page can be edited (see _templates/autoapi/index.rst).
# Any template not overridden here falls back to the one shipped with autoapi.
autoapi_template_dir = "_templates/autoapi"

# The docstrings still carry unfinished ".. todo::" notes; render them rather than
# silently dropping them, so that they stay visible to whoever picks the work up.
todo_include_todos = True


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
#
# "features" is a package, so AutoAPI generates a page for the package itself
# carrying nothing but a toctree of its submodules. Dropping it lets the
# "Feature measurements" landing page (see _templates/autoapi/index.rst) link
# straight to each feature, rather than through an intermediate list.
exclude_patterns = ['autoapi/features/index.rst']

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

# Third-party packages that are either licensed (Modeller), optional, or only
# installed on the machines that run the corresponding feature. Mocking them lets
# autodoc import every module in src/ without having them present.
autodoc_mock_imports = [
    "modeller",
    "nglview",
    "torchani",
    "esm",
    "melodia_py",
    "frustratometer",
    "pkaani",
]
