r"""
pyDMS.evaluate

Copyright 2025 Massachusetts Institute of Technology
Licensed under the MIT License
"""

import numpy as np
import statsmodels.api as sm
import warnings
import copy
import pyDMS
from scipy.optimize import minimize_scalar


def isotherm(gas, index):
    """Computes an isotherm from the optimized DMS parameters at
            the specified temperature.

    Args:
        gas: An instance of the Gas class.
        index: A float of the index of the temperature to compute
            the sorption isotherm at.

    Returns:
        An array of the isotherm
            [pressure or fugacity, concentration, concentration error].
    """

    # TODO: determine whether we should switch from index to temp <- user input needed

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
    """Computes the sorption coefficient at infinite dilution along with the corresponding uncertainty
    Results can be found in Gas.analysis.S_inf and Gas.analysis.S_inf_err

    Args:
        gas: an instance of the Gas class

    Returns:
        None
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
    """Computes heats of sorption
    
    Args:
        gas: an instance of the Gas class

    Returns:
        an instance of the Gas class
    """

    temp = gas.temp

    if len(gas.temp) == 2:
        # ignore variance warning when DOF=2-2=0
        warnings.filterwarnings(
            "ignore", message="divide by zero encountered in scalar divide", category=RuntimeWarning
        )

    inv_RT = 1 / (0.008314 * temp)

    def infinite_dilution_heat(gas):
        """Computes the infinite dilution heat of sorption
        
        Args:
            gas: an instance of the Gas class

        Returns:
            an instance of the Gas class
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
        """Computes the Henry heat of sorption
        
        Args:
            gas: an instance of the Gas class

        Returns:
            an instance of the Gas class
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
        """Computes the Langmuir heat of sorption
        
        Args:
            gas: an instance of the Gas class

        Returns:
            an instance of the Gas class
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


def isosteric_heat(gas, n_points):
    """Computes the isosteric heats of sorption
    
    Args:
        gas: An instance of the gas class

    Returns:
        None
    """

    def minimize(p, C_target, kD, CH, b):
        C_DMS = kD * p + CH * b * p / (1 + b * p)
        return (C_target - C_DMS) ** 2

    max_C = np.max(gas.C)

    C_target_low = np.linspace(0.001, 1, n_points)
    C_target_high = np.linspace(1.01, max_C, n_points)
    C_target = np.union1d(C_target_low, C_target_high)

    # p_guess = copy.deepcopy(C_target)

    p_iso = np.zeros((len(gas.temp), len(C_target)))
    C_iso = np.zeros(len(gas.temp), len(C_target))
    deltaH_iso = np.zeros(len(gas.temp), len(C_target))
    deltaH_iso_err = np.zeros(len(gas.temp), len(C_target))

    for i, C_val in enumerate(C_target):

        p_results = np.zeros(len(gas.temp))
        for j, temp_val in enumerate(gas.temp):
            kD = gas.kD[j]
            CH = gas.CH[j]
            b = gas.b[j]
            result = minimize_scalar(
                minimize, args=(C_val, kD, CH, b), bounds=(1e-6, 1000), method="bounded"
            )
            p_results[j] = result.x

        ln_p = np.log(p_results)
        inv_temp = 1 / gas.temp
        inv_temp_with_const = sm.add_constant(inv_temp)
        isosteric_model = sm.OLS(ln_p, inv_temp_with_const).fit()

        int_isosteric, slope_isosteric = isosteric_model.params
        int_isosteric_err, slope_isosteric_err = isosteric_model.bse
        deltaH_isosteric = -slope_isosteric * 8.313 * 10**-3  # *z
        deltaH_isosteric_err = slope_isosteric_err * 8.314 * 10**-3  # *z
        p_avg = np.mean(p_results)

        p_iso[i, j] = p_avg
        C_iso[i, j] = C_val
        deltaH_iso[i, j] = deltaH_isosteric
        deltaH_iso_err[i, j] = deltaH_isosteric_err

    gas.analysis.C_iso = C_iso
    gas.analysis.deltaH_iso = deltaH_iso
    gas.analysis.deltaH_iso_err = deltaH_iso_err
