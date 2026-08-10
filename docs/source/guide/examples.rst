Examples
========

The best way to understand how to use pyDMS is to see it applied. We provide examples that cover all of the main aspects we envision someone using pyDMS to require. The hope is that these examples can serve as starting points for users to modify with their data. All example code is also provided as Python files within the ``examples`` directory of pyDMS.

Example 1: Fitting DMS Parameters
---------------------------------

The main purpose of pyDMS is to fit DMS parameters. pyDMS is configured to have built-in solvers and settings that work for experimental sorption isotherms for a variety of different gases. To demonstrate, we fit data from Dean et al. :cite:`dean_elucidating_2024` for sorption in PIM-1 at temperatures of 35 °C, 45 °C, 55 °C, and 65 °C. It is important to note that pyDMS requires sorption at multiple temperatures to achieve its optimization routines. While in theory only two temperatures are needed, additional temperatures naturally improve fit.

To run the optimization, create a file called ``example1_pyDMS.py`` and copy the following code:

.. literalinclude:: ../../../examples/example_1/example1_pyDMS.py
   :language: python

Save the file and in your terminal, run the command:

.. code:: shell-session

   python3 example1_pyDMS.py

You should begin to see text appearing in your terminal, detailing what step pyDMS is on. This optimization will probably take a few minutes to run. If all is successful, the following information should be printed to the terminal:

.. code:: shell-session

   ==============================================================
                 ___  __  _______
      ___  __ __/ _ \/  |/  / __/      Copyright (C) 2026
     / _ \/ // / // / /|_/ /\ \        Massachusetts Institute
    / .__/\_, /____/_/  /_/___/        of Technology
   /_/   /___/

   Authors: B.C. Tapia, P.A. Dean, J.Y. Yeo, A.X. Wu, Z.P. Smith
   Web: https://smithlab.mit.edu
   ==============================================================
   ------------------------Fugacity Check------------------------
   Fugacity data supplied in Gas.f by user
   ---------------------LFER Initial Guesses---------------------
   kd0_0: [0.001, 0.01]
   dHd0_0: [-1, -30]
   b0_0: [0.0001, 0.005]
   dHb0_0: [-1, -30]
   C_H'0: [  0 100]
   C_H'1: [  0 100]
   C_H'2: [  0 100]
   C_H'3: [  0 100]
   ---------------------LFER Solver Bounds-----------------------
   kd0_0: [0, None]
   dHd0_0: [-50, 0]
   b0_0: [0, None]
   dHb0_0: [-50, 0]
   C_H'0: [  0 150]
   C_H'1: [  0 150]
   C_H'2: [  0 150]
   C_H'3: [  0 150]
   --------------------------------------------------------------
   LFER trial: 0/1000
   ...
   LFER trial: 950/1000
   ------------------van't Hoff Initial Guesses------------------
   dHD0_0: [-1, -30]
   dHb0_0: [-1, -30]
   C_H'0: [  0 100]
   C_H'1: [  0 100]
   C_H'2: [  0 100]
   C_H'3: [  0 100]
   -------------------van't Hoff Solver Bounds-------------------
   dHD0_0: [-50, 0]
   dHb0_0: [-50, 0]
   C_H'0: [  0 150]
   C_H'1: [  0 150]
   C_H'2: [  0 150]
   C_H'3: [  0 150]
   --------------------------------------------------------------
   van't Hoff trial: 0/1000
   ...
   van't Hoff trial: 950/1000
   --------------------------------------------------------------
   Optimization successful; starting error propagation
   Error propagation successful
   --------------------------------------------------------------
   Computing sorption energetics
   Calculating Z using built-in Virial EoS
   --------------------------------------------------------------
   Pickling example1_pyDMS.pkl
   Pickling successful
   --------------------------------------------------------------
   Generating report example1_pyDMS.pdf
   pyDMS successful!
   =======================END OF PROGRAM========================

Example 2: Understanding pyDMS Results
--------------------------------------

In the folder where you ran ``example1_pyDMS.py``, you should see two new files: ``example1_pyDMS.pdf`` and ``example1_pyDMS.pkl``. The results of the optimization are shown in ``example1_pyDMS.pdf``. These results should look similar to the example report provided in the following pages, which walks you through how to interpret the results.

If you open ``example1_pyDMS.pkl``, you will see non-human-readable results. That is because ``example1_pyDMS.pkl`` serializes and stores the results of the optimization such that they can be read back into Python from "Pickle" format. :cite:`pilgrim_serializing_2009` As an example, copy the following code into a new file called ``example2_pyDMS.py`` within the same folder as ``example1_pyDMS.pkl``:

.. literalinclude:: ../../../examples/example_2/example2_pyDMS.py
   :language: python

This should print something like the following to the Terminal:

.. code:: shell-session

   The value of b at 308 K is 0.6487 atm^-1

   These results come from the stored arrays of:
       temperature: [308 318 328 338]
       b: [0.64867561 0.48218935 0.364974   0.28084285]

To perform further calculations in Python, you don't have to rerun your optimization every time. Instead, you can quickly load a previously run optimization and have all the results automatically populated. A full list of all the parameters stored by the optimized gas class is available in :numref:`table-gas-class`.

Example 3: Converting from Pressure to Fugacity
-----------------------------------------------

If you examine Example 1, you will see that the fugacity was specified as the independent variable. Inspecting :numref:`table-gas-class`, however, you will see that both ``Gas.f`` and ``Gas.p`` can be supplied. This is because pyDMS can be run in three separate ways, listed in order of preference:

#. The user supplies ``Gas.f`` as the independent variable (fugacity-based DMS parameters).

#. The user supplies ``Gas.p`` and pyDMS converts the pressure to fugacity to use as the independent variable (fugacity-based DMS parameters).

#. The user supplies ``Gas.p`` as the independent variable (pressure-based DMS parameters).

Method 1 was already shown in Example 1. Method 2 is indentical to Method 1, however, ``Gas.f`` data is replaced with ``Gas.p`` data. Method 3 can be run by removing the ``Gas.formula`` variable as pyDMS will not compute fugacity without knowledge of the penetrant. pyDMS calculates fugacity either via the Virial :cite:`dymond_virial_2002` or the Peng--Robinson :cite:`peng_new_1976` equations of state (EoS), depending on the penetrant specified. :numref:`table-fugacity` shows the built-in gases within pyDMS. If the gas of interest is not available, the user can supply custom parameters for either EoS, shown below. An in-depth explanation of how Virial and Peng--Robinson coefficients must be supplied is available in :numref:`sec-fugacity`.

To show how we can convert from pressure to fugacity, create a file called ``example3_pyDMS.py`` and copy the following code:

.. literalinclude:: ../../../examples/example_3/example3_pyDMS.py
   :language: python

Example 4: Computing Pure- and Mixed-Gas Isotherms and Sorption Selectivities
-----------------------------------------------------------------------------

The dual-mode sorption parameters that pyDMS outputs can be used to calculate both the pure- and mixed-gas isotherms from the DMS model directly within pyDMS. Additionally, the sorption selectivity can be calculated from the mixed-gas isotherms. To show how we can compute pure- and mixed-gas isotherms and sorption selectivity, create a file called ``example4_pyDMS.py`` and copy the code below. The required .pkl files can be found in the ``pyDMS/Examples/Example4`` folder perhaps most easily accessible on GitHub.

.. literalinclude:: ../../../examples/example_4/example4_pyDMS.py
   :language: python

Example 5: Adjusting Default Parameters
---------------------------------------

While the default settings within pyDMS successfully optimize different gases without needing to manually change them in microporous polymers, there are instances when a user may wish to modify the defaults, such as when the parameter ranges are too large and ``pyDMS`` cannot converge. The easiest method of performing these parameter changes, is to use the helper function, ``pyDMS.dms.parameter_assist()``, which will automatically modify parameter ranges based on the isotherms supplied. Below, we demonstrate how to incorporate this helper into your code, using sorption in PIM-1 as an example:

.. literalinclude:: ../../../examples/example_5/example5_pyDMS.py
   :language: python

Alternatively, rather than use the helper, settings can be modified manually through the ``Gas.settings`` attribute. All settings that can be changed, and their default parameters, are shown in :numref:`table-settings`. Below, we provide an example demonstrating how to modify some settings:

.. code:: python

   import numpy as np
   import pyDMS.dms as dms


   # Data from:
   # Dean, P. A.; Mizrahi Rodriguez, K.; Guo, S.; Roy, N.; Swager, T. M.; Smith,
   # Z. P. Elucidating the role of micropore-generating backbone motifs and amine
   # functionality on H2S, CO2, CH4 and N2 sorption. Journal of Membrane Science
   # 2024, 696, 122465.

   n2 = dms.Gas()

   n2.formula = "N2"

   n2.f = [np.array([0.3, 0.7, 1.2, 1.9, 2.8, 3.7, 4.8, 6.3, 9.2, 12.4, 18.3, 24.8,
                     31.5, 38.2, 44.9]),
           np.array([0.3, 0.7, 1.2, 1.9, 2.8, 3.8, 4.8, 6.3, 9.2, 12.5, 18.3, 24.8,
                     31.6, 38.3, 45.1]),
           np.array([0.3, 0.6, 1.2, 1.9, 2.8, 3.8, 4.8, 6.3, 9.2, 12.4, 18.3, 24.9,
                     31.6, 38.4, 45.2]),
           np.array([0.2, 0.6, 1.3, 1.9, 2.8, 3.8, 4.8, 6.3, 9.2, 12.5, 18.3, 24.9,
                     31.7, 38.5, 45.4]),
   ]

   n2.c = [np.array([0.7, 1.4, 2.5, 3.9, 5.5, 7.3, 8.9, 11.2, 14.9, 18.7, 24.2, 29.2,
                     33.4, 37.0, 40.1]),
           np.array([0.5, 1.2, 2.1, 3.2, 4.6, 6.2, 7.6, 9.5, 13.0, 16.3, 21.4, 25.9,
                     29.9, 33.3, 36.2]),
           np.array([0.4, 0.9, 1.8, 2.6, 3.9, 5.2, 6.5, 8.1, 11.1, 14.1, 18.7, 22.9,
                     26.6, 29.8, 32.6]),
           np.array([0.3, 0.8, 1.5, 2.3, 3.3, 4.4, 5.5, 6.9, 9.6, 12.2, 16.2, 20.1,
                     23.5, 26.3, 29.1])]

   n2.c_err = [np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3,
                         0.3, 0.3, 0.4]),
               np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3,
                         0.3, 0.3, 0.3]),
               np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                         0.3, 0.3, 0.3]),
               np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2,
                         0.3, 0.3, 0.3])]

   n2.temp = np.array([308, 318, 328, 338])

   n2.settings = {'CH_guess': [[30, 70], [30, 60], [30, 50], [30, 40]],
                  'CH_bounds': [[20, 80], [20, 70], [20, 60], [20, 50]],
                  'dHD_guess': [-15, -1],
                  'dHD_bounds': [-20, 0],
                  'dHb_guess': [-15, -1],
                  'dHb_bounds': [-20, 0],
                  'trials': 1000,
                  'seed': 34564}

   dms.compute(n2, 'n2')

As the paper describing ``pyDMS`` analyzed seven different glassy polymers, there are different settings that can be tried if pyDMS is struggling to converge. These different settings are described in the paper's supplemental information with the Python code also available in the ``examples/paper`` directory in pyDMS, which can be found on GitHub. If you suspect there may be convergence challenges, we suggest the following strategy:

#. Run pyDMS with the default parameters and see if the optimizations converge.

#. If you have a non-microporous polymer or a lower sorbing penetrant, reduce the upper bounds of ``CH_guess`` and ``CH_bounds``.

#. Inspect the LFERs (report page 3) and adjust the :math:`\Delta H_\mathrm{D}` and :math:`\Delta H_\mathrm{b}` initial guesses and bounds (``deltaHD_guess``, ``deltaHD_bounds``, ``deltaHb_guess``, and ``deltaHb_bounds``) such that the ranges considered encompass the successful optimizations.

#. Inspect the LFERs and adjust the :math:`k_\mathrm{D,0}` and :math:`b_\mathrm{0}` initial guesses and bounds (``kD0_guess``, ``kD0_bounds``, ``b0_guess``, and ``b0_bounds``) such that the ranges considered encompass the successful optimizations.

#. Increase ``trials`` to a value higher than 1000.

After each step, we recommend rerunning the optimization to see if convergence is achieved.
