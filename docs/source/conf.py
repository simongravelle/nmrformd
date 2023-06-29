# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'NMRforMD'
copyright = 'All source code is available under the GNU General Public License v3.0'
author = 'Simon Gravelle'

version = '0.1'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'nbsphinx',
    'sphinx_favicon',
]

# Execute the notebooks
nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
exclude_patterns = ['_build', '**.ipynb_checkpoints']

templates_path = ['_templates']

exclude_patterns = ['_build', '**.ipynb_checkpoints']

pygments_style = 'tango'

html_theme = 'furo'
html_title = "    "

html_static_path = ['_static']
html_css_files = ["custom.css"]

favicons = [
    {"href": "favicon-32x32.png"},
]

html_logo = "figures/logo/logo.png"
