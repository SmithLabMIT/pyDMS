# pyDMS
### Python package for dual-mode sorption (DMS) analysis

<p align="center">
<img src="./images/pyDMS_logo.png" width="300" class="center">
</p>

pyDMS is a Python program for the computation of dual-mode sorption (DMS) parameters using linear free energy relationship (LFER) and Van't Hoff constraints. 

pyDMS was developed in the [Smith Lab](https://smithlab.mit.edu/) in the Department of Chemical Engineering at the Massachusetts Institute of Technology.

The paper describing pyDMS can be found [here](LINK)

### Features
- Compute reproducible DMS parameters: $C_\mathrm{H}^\prime$, $k_\mathrm{D}$, and $\mathrm{b}$ *with* uncertainty via LFERs and van't Hoff constraints.
- Access computed results in Python objects or via the automatically output PDF.
- Compute pure and mixed-gas sorption isotherms from the DMS parameters *with* error propogation
- Calculate ideal and mixed-gas sorption selectivities *with* error uncertainty
- Calculate fugacity from pressure via Virial and Peng-Robinson implementation


## Citation(s)
#### If you used pyDMS, please cite:
Tapia, B. C.; Dean, P. A.; Yeo, J. Y.; Smith, Z. P. pyDMS: A Python package for the determination of physics-informed, reproducible dual-mode sorption (DMS) parameters. *AIChE J.* **2026**.

#### We also recommend citing the following works which provide relevant background and theory on van't Hoff sorption energetics, LFERs, and their application to constraied DMS optimization, respectively:
Koros, W. J.; Paul, D. R.; Huvard, G. S. Energetics of Gas Sorption in Glassy Polymers. *Polymer* **1979**, *20* (8), 956–960. https://doi.org/10.1016/0032-3861(79)90192-7.

Freeman, B. D. Basis of Permeability/Selectivity Tradeoff Relations in Polymeric Gas Separation Membranes. *Macromolecules* **1999**, *32* (2), 375–380. https://doi.org/10.1021/ma9814548.

Wu, A. X.; Drayton, J. A.; Mizrahi Rodriguez, K.; Benedetti, F. M.; Qian, Q.; Lin, S.; Smith, Z. P. Elucidating the Role of Fluorine Content on Gas Sorption Properties of Fluorinated Polyimides. *Macromolecules* **2021**, *54* (1), 22–34. https://doi.org/10.1021/acs.macromol.0c01746.

## License
Copyright (C) 2026 The Massachusetts Institute of Technology

This work is licensed under the MIT License (see "LICENSE")



## Acknowledgements
This work was supported by a MathWorks Fellowship, NSF CAREER Award (no. 2146422), and the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Separation Science Program under Award DE-SC0023252.
