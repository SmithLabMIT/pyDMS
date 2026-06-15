pyDMS
=====

Python package for dual-mode sorption (DMS) analysis
-----------------------------------------------------

.. image:: ./images/pyDMS_logo.png
   :width: 300
   :align: center

pyDMS is a Python program for the computation of dual-mode sorption (DMS) parameters using linear free energy relationship (LFER) and Van't Hoff constraints.

pyDMS is developed in the `Smith Lab <https://smithlab.mit.edu/>`_ in the `Department of Chemical Engineering <https://cheme.mit.edu/>`_ at the `Massachusetts Institute of Technology <https://mit.edu/>`_.

The paper describing pyDMS can be found `here <LINK>`_.

Features
--------

- Compute reproducible DMS parameters: :math:`C_\mathrm{H}^\prime`, :math:`k_\mathrm{D}`, and :math:`b` with uncertainty via LFERs and van't Hoff constraints.
- Access computed results in Python objects or via the automatically output PDF.
- Compute pure and mixed-gas sorption isotherms from the DMS parameters with uncertainty values.
- Calculate ideal and mixed-gas sorption selectivities with uncertainty values.
- Calculate fugacity from pressure via Virial and Peng--Robinson implementation.

Installation
------------

pyDMS can be installed using:

- conda: ``conda install package-not-yet-available``
- pip: ``pip install package-not-yet-available``

Alternatively, users can install the latest version from GitHub, either directly:

.. code-block:: bash

   python -m pip install git+https://github.com/SmithLabMIT/pyDMS.git

or by downloading the source code:

.. code-block:: bash

   git clone https://github.com/SmithLabMIT/pyDMS.git
   cd pyDMS
   python -m pip install .

For developers, it is recommended to download the source code and install in editable mode by replacing ``python -m pip install .`` with ``python -m pip install -e .``.

Citation(s)
-----------

**If you used pyDMS, please cite:**

Tapia, B. C.; Dean, P. A.; Yeo, J. Y.; Wu, A. X., Smith, Z. P. pyDMS: A Python package for the determination of physics-informed, reproducible dual-mode sorption (DMS) parameters. *In prep.* **2026**.

**We also recommend citing the following works which provide relevant background and theory on van't Hoff sorption energetics, LFERs, and their application to constrained DMS optimization, respectively:**

Koros, W. J.; Paul, D. R.; Huvard, G. S. Energetics of Gas Sorption in Glassy Polymers. *Polymer* **1979**, *20* (8), 956–960. https://doi.org/10.1016/0032-3861(79)90192-7.

Freeman, B. D. Basis of Permeability/Selectivity Tradeoff Relations in Polymeric Gas Separation Membranes. *Macromolecules* **1999**, *32* (2), 375–380. https://doi.org/10.1021/ma9814548.

Wu, A. X.; Drayton, J. A.; Mizrahi Rodriguez, K.; Benedetti, F. M.; Qian, Q.; Lin, S.; Smith, Z. P. Elucidating the Role of Fluorine Content on Gas Sorption Properties of Fluorinated Polyimides. *Macromolecules* **2021**, *54* (1), 22–34. https://doi.org/10.1021/acs.macromol.0c01746.

License
-------

Copyright © 2026 The Massachusetts Institute of Technology

This work is licensed under the 3-clause BSD license (see ``LICENSE``).

Acknowledgements
----------------

This work was supported by a MathWorks Fellowship, NSF CAREER Award (no. 2146422), and the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Separation Science Program under Award DE-SC0023252.