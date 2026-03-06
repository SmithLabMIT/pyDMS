import numpy as np
import pyDMS.dms as dms

# The name you wish this optimization to be called.
GAS_NAME = dms.Gas()

# Use if you have pressure-based isotherms.
# The pressure points (atm) of each isotherm go within each np.array([]).
# Add more np.array([])'s if you have more than two isotherms
GAS_NAME.p = [np.array([2, 4, 6, 8]), np.array([1, 3, 7, 9])]

# Use if you have fugacity-based isotherms.
# The fugacity points (atm) of each isotherm go within each np.array([]).
# Add more np.array([])'s if you have more than two isotherms
GAS_NAME.f = [np.array([2, 4, 6, 8]), np.array([1, 3, 7, 9])]

# The concentration points (cm^3 cm^-3) of each isotherm go within each np.array([]).
# Add more np.array([])'s if you have more than two isotherms
GAS_NAME.c = [np.array([5, 10, 15, 20]), np.array([2, 7, 10, 13])]

# The std. dev. of the conc. points (cm^3 cm^-3) of each isotherm go within np.array([]).
# Add more np.array([])'s if you have more than two isotherms
GAS_NAME.c_err = [np.array([0.5, 0.7, 0.9, 1.0]), np.array([0.2, 0.3, 0.5, 0.7])]

# The temperatures (K) each isotherm was run at.
GAS_NAME.temp = np.array([308, 328])

# The formula of the gas (e.g., "CO2").
# Use this if you wish to convert provided pressures to fugacities.
# If you provide fugacities knowledge of the compound is not needed.
# If you provide pressures and want pressure-based DMS params., do not provide a formula
GAS_NAME.formula = "FORMULA_HERE"

# These are all the settings that you can change within pyDMS.
# If you do not provide any settings, these settings will automatically be applied
# The manual explains the role of each settings.
GAS_NAME.settings = {"dHD_guess":      [-30, -1],
                     "dHD_bounds":     [-50,  0],
                     "dHb_guess":      [-30, -1],
                     "dHb_bounds":     [-50,  0],
                     "kD0_guess":      [0, None],
                     "kD0_bounds":     [0.001, 0.01],
                     "b0_guess":       [0, None],
                     "b0_bounds":      [0.001, 0.01],
                     "CH_guess":       [[0, 100], [0, 100]],  # One subarray per temp
                     "CH_bounds":      [[0, 150], [0, 150]],  # One subarray per temp
                     "trials":         1000,
                     "solver_LFER":    "SLSQP",
                     "solver_vH":      "SLSQP",
                     "maxiter_LFER":   1000,
                     "maxiter_vH":     1000,
                     "ftol":           1E-7,
                     "xtol":           1E-7,
                     "gtol":           1E-7,
                     "verbose":        True,
                     "solver_verbose": False,
                     "seed":           None}

# If you supply pressure and want to compute fugacity using a Virial expansion but
# your gas is not He, H2, N2, O2, CH4, CO2, C2H6, C2H4, C3H8, or C3H6, then
# uncomment and add virial coefficients here.
# GAS_NAME.virial_coeff = {"B0":0, "B1":0, "B2":0, "B3":0, "B4":0, # B [=] cm^3 mol^-1
#                          "C0":0, "C1":0, "C2":0, "C3":0, "C4":0} # C [=] cm^6 mol^-2

# If you supply pressure and want to compute fugacity using the Peng-Robinson EoS but
# your gas is not H2S, then
# uncomment and add virial coefficients here.
# GAS_NAME.pr_coeff = {"Tc": 0, "Pc": 0, "omega": } Tc [=] K, Pc [=] MPa

# Run the code! (will save FILE_NAME.pdf and FILE_NAME.pkl)
dms.gas(GAS_NAME, "FILE_NAME")
