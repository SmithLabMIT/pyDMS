# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import sys
import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

# -- Project information -----------------------------------------------------

project = 'pyDMS'
copyright = f'{datetime.date.today().year}, Massachusetts Institute of Technology'
author = 'Massachusetts Institute of Technology'
release = '0.1.0'

language = "en"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
]


bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "super"

templates_path = ['_templates']
exclude_patterns = [
    "_build",
    "_templates",
]

show_authors = True
pygments_style = "sphinx"
todo_include_todos = False



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_theme_options = {
    "sidebar_hide_name": True,
}

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_remove_prompts = True
copybutton_line_continuation_character = "\\"

html_logo = "images/pyDMS_logo.png"
html_favicon = "images/pyDMS_favicon.png"
html_static_path = ['_static']
html_show_sourcelink = True


