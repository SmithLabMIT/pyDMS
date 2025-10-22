import numpy as np
import pyDMS.dms as dms
import pyDMS.multigas as mg
import pyDMS.evaluate as ev


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

# The std. dev. of the conc. points (cm^3 cm^-3) of each isotherm go within each np.array([]).
# Add more np.array([])'s if you have more than two isotherms
GAS_NAME.c_err = [np.array([0.5, 0.7, 0.9, 1.0]), np.array([0.2, 0.3, 0.5, 0.7])]

# The temperatures (K) each isotherm was run at.
GAS_NAME.temp = [308, 328]

# The formula of the gas (e.g., "CO2").
# Use this if you wish to convert provided pressures to fugacities.
# If you provide fugacities knowledge of the compound is not needed.
# If you provide pressures and want pressure-based DMS parameters, do not provide a formula.
GAS_NAME.formula = "FORMULA_HERE"

# These are all the settings that you can change within pyDMS.
# If you do not provide any settings, the settings shown here will automatically be applied
# The manual explains the role of each settings.
# The manual shows settings used for CO2 sorption in different polymers.
GAS_NAME.settings = {"dHD_guess":      [-30, -1],
                     "dHD_bounds":     [-50,  0],
                     "dHb_guess":      [-30, -1],
                     "dHb_bounds":     [-50,  0],
                     "kD0_guess":      [0, None],
                     "kD0_bounds":     [0.001, 0.01],
                     "b0_guess":       [0, None],
                     "b0_bounds":      [0.001, 0.01],
                     "CH_guess":       [0, 100],
                     "CH_bounds":      [0, 100],
                     "trials":         1000,
                     "solver_LFER":    "SLSQP",
                     "solver_vH":      "SLSQP",
                     "maxiter_LFER":   1000,
                     "maxiter_vH":     1000,
                     "ftol":           1E-12,
                     "xtol":           1E-12,
                     "gtol":           1E-12,
                     "verbose":        True,
                     "solver_verbose": False}

dms.gas(GAS_NAME, "FILE_NAME")
