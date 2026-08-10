Installation
============

pyDMS can be downloaded through PIP, Anaconda, or GitHub. It is developed to be OS-agnostic and has been tested and validated on Windows 11, Ubuntu 24, and MacOS (Intel and ARM).

To check if Python is already installed, open the "Terminal" (Windows, Mac, and Linux) or "Command prompt" (Windows) program and run the command

.. code:: shell-session

   python --version

Python is installed if the output says something like

.. code:: shell-session

   Python 3.9.X

If the returned version is anything lower than ``3.9``, an updated version of Python is required. If the returned version is ``3.14`` or higher, a lower version may be required (as not all pyDMS dependencies may support it yet). If an error occurs after running the above command, try

.. code:: shell-session

   python3 --version

instead. If an error occurs here, Python must be installed. Python can be installed from `www.python.org/downloads <https://www.python.org/downloads/>`__, or if you prefer to use Anaconda (see below), from `www.anaconda.com/download <https://www.anaconda.com/download>`__.

PIP
---

The easiest method to get started with pyDMS is by downloading from PyPi, using pip. As long as Python (version ``3.9``--``3.13``) is available on your system, you should be able to download pyDMS by opening up a "Terminal" or "Command prompt" and running

.. code:: shell-session

   pip install pyDMS-sorption

To check that it was successfully installed, in the terminal you should be able to run

.. code:: shell-session

   python -c "import pyDMS; print(pyDMS.__version__)"

which will return the version of pyDMS that was installed.

Anaconda
--------

pyDMS is available in the Conda-Forge Anaconda channel. Once an Anaconda installer is downloaded to the system (e.g., Anaconda, Miniconda, or Miniforge), pyDMS can be installed with:

.. code:: shell-session

   conda install conda-forge::pyDMS-sorption

GitHub
------

The source code for pyDMS is available on GitHub. Installation from GitHub can be done from the terminal as:

.. code:: shell-session

   pip install git+https://github.com/SmithLabMIT/pyDMS

This is helpful if there are new updates to pyDMS that have not yet been incorporated into the stable version on PyPI and Conda-Forge. Additionally, users can asks questions, report bugs, and contribute to the development of pyDMS on the GitHub page.

To download the source code from GitHub and install it locally in editable mode, installation can be done as follows:

.. code:: shell-session

   git clone https://github.com/SmithLabMIT/pyDMS
   cd pyDMS
   pip install . -e ".[dev,docs]"

where the ``-e`` flag allows for editable mode, and the ``.[dev,docs]`` allows for installation of the development and documentation dependencies. This is useful if you want to contribute to the development of pyDMS or build the documentation locally.
