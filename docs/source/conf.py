"""Configuration file for the EPFT-UP Sphinx documentation builder.

Pattern copied from IOPtics (docs/source/conf.py @ develop) with one addition:
MyST, so the Markdown documents under reports/ and docs/design/ are published
unchanged through thin include-shim pages (see models/moana/*.md).

See https://www.sphinx-doc.org/en/master/usage/configuration.html for the full
list of built-in configuration values.
"""

import os
import sys

# Make the epft_up package importable for autodoc. conf.py lives at
# docs/source/, so the repository root is two levels up.
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'EPFT-UP'
copyright = '2026, J. Xavier Prochaska and Claude'
author = 'J. Xavier Prochaska and Claude'

try:                                   # single-source the version from the package
    import epft_up
    release = epft_up.__version__
except Exception:
    release = '0.0.dev0'
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',                 # Markdown sources (reports, design docs)
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',         # NumPy / Google style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx_design',
]

source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
master_doc = 'index'

# The included Markdown documents were written for GitHub rendering: they use
# non-consecutive heading levels in places and cross-reference each other by
# repo path. Neither is an error for the site.
suppress_warnings = ['myst.header', 'myst.xref_missing']

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = [
    'colon_fence',        # ::: admonitions in Markdown shims
    'deflist',
    'attrs_inline',
]
myst_heading_anchors = 3  # GitHub-style #anchors for the reports' internal links

# -- Autodoc / autosummary ---------------------------------------------------

autosummary_generate = True

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'undoc-members': False,
    'show-inheritance': True,
}
autodoc_typehints = 'description'      # built-in (no sphinx-autodoc-typehints dep)

# `earthaccess` touches the network at import-time paths we never exercise in
# a docs build; everything else EPFT-UP imports is pip-installed by RTD.
autodoc_mock_imports = ['earthaccess']

# Napoleon settings
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# -- HTML output -------------------------------------------------------------

# Furo — a clean, modern theme (light/dark), as on the IOPtics site.
html_theme = 'furo'
html_static_path = ['_static']
html_title = 'EPFT-UP'
html_css_files = ['custom.css']       # ocean-colour theme accents

html_theme_options = {
    'sidebar_hide_name': False,
    'navigation_with_keys': True,
    'source_repository': 'https://github.com/ocean-colour/EPFT-UP/',
    'source_branch': 'main',
    'source_directory': 'docs/source/',
    'light_css_variables': {
        'color-brand-primary': '#1a6ea8',
        'color-brand-content': '#1a6ea8',
    },
    'dark_css_variables': {
        'color-brand-primary': '#6cb6e8',
        'color-brand-content': '#6cb6e8',
    },
}
