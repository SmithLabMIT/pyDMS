"""
pyDMS.evaluate

Copyright 2026 Massachusetts Institute of Technology
Licensed under the 3-clause BSD license
"""

import numpy as np
import statsmodels.api as sm
import warnings
#import copy
import pyDMS
from scipy.optimize import minimize_scalar


def isotherm(gas, temp):
    """Computes an isotherm from the optimized DMS parameters at the specified temperature.

    Args:
        gas: An instance of the `Gas` class.
        temp: The temperature of the desired isotherm.

    Returns: An array of the isotherm [pressure or fugacity, concentration, concentration error].
    """

    idx = np.where(gas.temp == temp)[0]

    if idx.size:
        index = idx[0]
    else:
        index = int(np.argmin(np.abs(gas.temp - temp)))
        pyDMS.warning_in_orange(f"Requested temperature {temp} not found. Using nearest value {gas.temp[index]} instead.")

    ch = gas.CH[index]
    kd = gas.kD[index]
    b = gas.b[index]
    ch_err = gas.CH_err[index] if gas.CH_err is not None else 0.0
    kd_err = gas.kD_err[index] if gas.kD_err is not None else 0.0
    b_err = gas.b_err[index] if gas.b_err is not None else 0.0

    if gas.f is not None:
        p_vec = gas.f
    elif gas.p is not None:
        p_vec = gas.p
    else:
        pyDMS.warning_in_orange("Neither p or f is defined, using isotherm range 0-40 atm")
        p_vec = np.tile([0.0, 40.0], (len(gas.temp), 1))

    p_val = np.linspace(1e-6, np.max(p_vec[index]), 3000)

    c = kd * p_val + ch * b * p_val / (1 + b * p_val)
    c_err = np.sqrt(
        ch**2 * p_val**2 * b_err**2 / (1 + b * p_val) ** 4
        + b**2 * p_val**2 * ch_err**2 / (1 + b * p_val) ** 2
        + p_val**2 * kd_err**2
    )

    return p_val, c, c_err


def S_inf(gas):
    """Computes the sorption coefficient at infinite dilution along with the corresponding
    uncertainty. Results can be found in `Gas.analysis.S_inf` and `Gas.analysis.S_inf_err`.

    Args:
        gas: an instance of the `Gas` class

    Returns: None
    """

    k_D = gas.kD
    b = gas.b
    ch = gas.CH

    k_D_err = gas.kD_err
    b_err = gas.b_err
    ch_err = gas.CH_err

    S_inf_calc = k_D + ch * b
    S_inf_err_calc = np.sqrt(ch**2 * b_err**2 + ch_err**2 * b**2 + k_D_err)

    gas.analysis.S_inf = S_inf_calc
    gas.analysis.S_inf_err = S_inf_err_calc


def heat_of_sorption(gas, method="all"):
    """Computes heats of sorption. Results can be found in `Gas.analysis.deltaH_S_inf`,
    `Gas.analysis.deltaH_D`, and `Gas.analysis.deltaH_b` for the infinite dilution, Henry, and
    Langmuir heats of sorption, respectively. Results are stored as lists of the form
    [slope, intercept] where slope is the heat of sorption and intercept is the pre-exponential
    factor. Uncertainties are stored in `Gas.analysis.deltaH_S_inf_err`,
    `Gas.analysis.deltaH_D_err`, and `Gas.analysis.deltaH_b_err` as lists of the form
    [slope error, intercept error].

    Args:
        gas: an instance of the `Gas` class

    Returns: an instance of the `Gas` class
    """

    temp = gas.temp

    if len(gas.temp) == 2:
        # ignore variance warning when DOF=2-2=0
        warnings.filterwarnings(
            "ignore", message="divide by zero encountered in scalar divide", category=RuntimeWarning
        )

    inv_RT = 1 / (0.008314 * temp)

    def infinite_dilution_heat(gas):
        """Computes the infinite dilution heat of sorption. Results are stored in
        `Gas.analysis.deltaH_S_inf` and `Gas.analysis.deltaH_S_inf_err

        Args:
            gas: an instance of the `Gas` class

        Returns:
            an instance of the `Gas` class
        """

        S_inf(gas)

        S_inf_calc = gas.analysis.S_inf
        # S_inf_err_calc = gas.analysis.S_inf_err
        # print(S_inf_calc)
        ln_S_inf = np.log(S_inf_calc)

        inv_RT_with_const = sm.add_constant(inv_RT)
        # print(ln_S_inf)
        deltaH_S_inf_model = sm.OLS(ln_S_inf, inv_RT_with_const).fit()
        # print(deltaH_S_inf_model.params)
        int_deltaH_S_inf, slope_deltaH_S_inf = deltaH_S_inf_model.params
        int_deltaH_S_inf_err, slope_deltaH_S_inf_err = deltaH_S_inf_model.bse
        # print(int_deltaH_S_inf, slope_deltaH_S_inf)
        gas.analysis.deltaH_S_inf = [-slope_deltaH_S_inf, np.exp(int_deltaH_S_inf)]
        gas.analysis.deltaH_S_inf_err = [
            slope_deltaH_S_inf_err,
            (np.exp(int_deltaH_S_inf) * int_deltaH_S_inf_err),
        ]

        return gas

    def henry_heat(gas):
        """Computes the Henry heat of sorption. Results are stored in `Gas.analysis.deltaH_D` and
        `Gas.analysis.deltaH_D_err`

        Args:
            gas: an instance of the `Gas` class

        Returns:
            an instance of the `Gas` class
        """

        k_D = gas.kD

        ln_k_D = np.log(k_D)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltak_D_model = sm.OLS(ln_k_D, inv_RT_with_const).fit()

        int_deltak_D, slope_deltak_D = deltak_D_model.params
        int_deltak_D_err, slope_deltak_D_err = deltak_D_model.bse

        if gas.analysis.deltaH_D:
            gas.analysis.deltaH_D[1] = np.exp(int_deltak_D)
            gas.analysis.deltaH_D_err[1] = np.exp(int_deltak_D) * int_deltak_D_err
        else:
            # print("pyDMS optimization not performed, using provided DMS parameters for Henry Energetics")
            gas.analysis.deltaH_D = [-slope_deltak_D, np.exp(int_deltak_D)]
            gas.analysis.deltaH_D_err = [
                slope_deltak_D_err,
                np.exp(int_deltak_D) * int_deltak_D_err,
            ]

        return gas

    def langmuir_heat(gas):
        """Computes the Langmuir heat of sorption. Results are stored in `Gas.analysis.deltaH_b` and `Gas.analysis.deltaH_b_err`

        Args:
            gas: an instance of the `Gas` class

        Returns:
            an instance of the `Gas` class
        """

        b = gas.b

        ln_b = np.log(b)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltab_model = sm.OLS(ln_b, inv_RT_with_const).fit()

        int_deltab, slope_deltab = deltab_model.params
        int_deltab_err, slope_deltab_err = deltab_model.bse

        # gas.analysis.deltaH_b = [-slope_deltab, np.exp(int_deltab)]
        # gas.analysis.deltaH_b_err = [slope_deltab_err,
        #   np.exp(int_deltab)*int_deltab_err]
        if gas.analysis.deltaH_b:
            gas.analysis.deltaH_b[1] = np.exp(int_deltab)
            gas.analysis.deltaH_b_err[1] = np.exp(int_deltab) * int_deltab_err
        else:
            # print("pyDMS optimization not performed, using provided DMS parameters for Langmuir Energetics")
            gas.analysis.deltaH_b = [-slope_deltab, np.exp(int_deltab)]
            gas.analysis.deltaH_b_err = [slope_deltab_err, np.exp(int_deltab) * int_deltab_err]

        return gas

    if method == "all":

        infinite_dilution_heat(gas)
        henry_heat(gas)
        langmuir_heat(gas)

    elif method == "infinite_dilution":
        infinite_dilution_heat(gas)

    elif method == "henry":
        henry_heat(gas)

    elif method == "langmuir":
        langmuir_heat(gas)

    else:
        pyDMS.error_in_red(f"Unknown method {method} for heat of sorption")

    return gas


def isosteric_heat(gas, n_points=1000):
    """Computes the isosteric heats of sorption. Results are stored in `Gas.analysis.deltaH_iso`
    and `Gas.analysis.deltaH_iso_err` as arrays of the same length as `Gas.analysis.c_iso` which
    contains the corresponding concentrations. The isosteric heat of sorption is computed by first
    determining the pressure/fugacity at which a given concentration is achieved at each
    temperature, then performing a van't Hoff analysis to determine the isosteric heat of sorption
    at that concentration. Uncertainties are determined from the standard error of the slope of the
    van't Hoff analysis.

    Args:
        gas: An instance of the `Gas` class
        n_points: The number of points to use in the isosteric heat calculation

    Returns:
        None
    """

    def minimize_obj(p, C_target, kD, CH, b):
        C_DMS = kD * p + CH * b * p / (1 + b * p)
        return (C_target - C_DMS) ** 2

    max_C = max(np.max(np.asarray(ci, dtype=float)) for ci in gas.c)

    C_target_low = np.linspace(0.1, 1, n_points)
    C_target_high = np.linspace(1.01, max_C, n_points)
    C_target = np.union1d(C_target_low, C_target_high)

    nT = len(gas.temp)
    nC = len(C_target)

    # storing p(T,C), C_iso(C), dH_iso(C)
    # p_iso = np.zeros((nT, nC))
    C_iso = np.zeros(nC)
    deltaH_iso = np.zeros(nC)
    deltaH_iso_err = np.zeros(nC)

    inv_temp = 1 / np.asarray(gas.temp)

    # can we compute Z?
    if gas.virial_coeff:
        print("Calculating Z using user-supplied Virial coefficients")
        calc_z = True

    elif gas.pr_coeff:
        print("Calculating Z using user-supplied Peng-Robinson " "coefficients")
        calc_z = True

    elif gas.formula in pyDMS.fugacity.virial_coeff:
        print("Calculating Z using built-in Virial EoS")
        calc_z = True

    elif gas.formula in pyDMS.fugacity.pr_coeff:
        print("Calculating Z using built-in Peng-Robinson EoS")
        calc_z = True

    elif gas.formula is None:
        pyDMS.warning_in_orange(
            "Gas.formula is not specfied. Assuming Z=1"
        )
        calc_z = False

    else:
        pyDMS.warning_in_orange(
            "Something went wrong trying to calculate Z.\n Will use Z=1"
        )
        calc_z = False

    for i, C_val in enumerate(C_target):

        p_results = np.zeros(nT)
        z_vals = np.zeros(nT)

        for j, T in enumerate(gas.temp):
            kD = gas.kD[j]
            CH = gas.CH[j]
            b = gas.b[j]

            result = minimize_scalar(
                minimize_obj, args=(C_val, kD, CH, b), bounds=(1e-6, 1000), method="bounded"
            )

            p = result.x
            p_results[j] = p

            # computing Z
            if calc_z:
                if gas.virial_coeff:
                    z_vals[j] = pyDMS.fugacity.virial_eos_z(gas, p, T)
                elif gas.pr_coeff:
                    z_vals[j] = pyDMS.fugacity.peng_robinson_eos_z(gas, p, T)
                elif gas.formula in pyDMS.fugacity.virial_coeff:
                    z_vals[j] = pyDMS.fugacity.virial_eos_z(gas, p, T)
                elif gas.formula in pyDMS.fugacity.pr_coeff:
                    z_vals[j] = pyDMS.fugacity.peng_robinson_eos_z(gas, p, T)

        ln_p = np.log(p_results)
        inv_temp_with_const = sm.add_constant(inv_temp)
        isosteric_model = sm.OLS(ln_p, inv_temp_with_const).fit()

        if calc_z is True:
            #print(z_vals)
            z = np.mean(z_vals)
            #print(z)
            #print(np.std(z_vals))
        else:
            z = 1.0

        int_isosteric, slope_isosteric = isosteric_model.params
        int_isosteric_err, slope_isosteric_err = isosteric_model.bse
        deltaH_isosteric = slope_isosteric * 8.314 * 10**-3 * z
        deltaH_isosteric_err = slope_isosteric_err * 8.314 * 10**-3 * z


        #p_iso[:, ] = p_avg
        C_iso[i] = C_val
        deltaH_iso[i] = deltaH_isosteric
        deltaH_iso_err[i] = deltaH_isosteric_err
    #from matplotlib import pyplot as plt
    #plt.errorbar(C_iso, deltaH_iso, yerr=deltaH_iso_err)
    #plt.show()
    gas.analysis.c_iso = C_iso
    gas.analysis.deltaH_iso = deltaH_iso
    gas.analysis.deltaH_iso_err = deltaH_iso_err

    return gas
