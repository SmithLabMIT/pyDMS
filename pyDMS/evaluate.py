r'''
pyDMS.evaluate

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import numpy as np
import statsmodels.api as sm
import warnings

import pyDMS


def isotherm(gas, index, eos=None):
    '''Computes an isotherm from the optimized DMS parameters at
            the specified temperature.

    Args:
        gas: An instance of the Gas class.
        index: A float of the index of the temperature to compute
            the sorption isotherm at.

    Returns:
        An array of the isotherm
            [pressure or fugacity, concentration, concentration error].
    '''

    ch = gas.CH[index]
    kd = gas.kD[index]
    b = gas.b[index]
    ch_err = gas.CH_err[index]
    kd_err = gas.kD_err[index]
    b_err = gas.b_err[index]
    if gas.f is not None:
        p_vec = gas.f
    else:
        p_vec = gas.p
    p_val = np.linspace(1E-6, np.max(p_vec[index]), 1000)

    c = kd*p_val+ch*b*p_val/(1+b*p_val)
    c_err = (
        np.sqrt(ch**2*p_val**2*b_err**2/(1+b*p_val)**4+b**2*p_val**2*ch_err**2
                / (1+b*p_val)**2+p_val**2*kd_err**2))

    return p_val, c, c_err


def S_inf(gas):
    '''Computes the sorption coefficient at infinite dilution

    Args:
        gas: an instance of the Gas class

    Returns:
        None
    '''

    k_D = gas.kD
    b = gas.b
    ch = gas.CH
    S_inf_calc = k_D+ch*b

    gas.analysis.S_inf = S_inf_calc
    gas.analysis.S_inf_err = 'Not yet implemented*'


def heat_of_sorption(gas, method='all'):
    '''Computes heats of sorption

    '''

    temp = gas.temp

    if len(gas.temp) == 2:
        # ignore variance warning when DOF=2-2=0
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide", category=RuntimeWarning)

    inv_RT = 1/(0.008314*temp)

    def infinite_dilution_heat(gas):
        '''

        '''
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
        gas.analysis.deltaH_S_inf = [-slope_deltaH_S_inf,
                                     np.exp(int_deltaH_S_inf)]
        gas.analysis.deltaH_S_inf_err = [slope_deltaH_S_inf_err,
                                         (np.exp(int_deltaH_S_inf)
                                          * int_deltaH_S_inf_err)]

        return gas

    def henry_heat(gas):
        '''
        *
        '''

        k_D = gas.kD

        ln_k_D = np.log(k_D)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltak_D_model = sm.OLS(ln_k_D, inv_RT_with_const).fit()

        int_deltak_D, slope_deltak_D = deltak_D_model.params
        int_deltak_D_err, slope_deltak_D_err = deltak_D_model.bse

        # gas.analysis.deltaH_D = [-slope_deltak_D, np.exp(int_deltak_D)]
        # gas.analysis.deltaH_D_err = [slope_deltak_D_err,
        #   np.exp(int_deltak_D)*int_deltak_D_err]

        gas.analysis.deltaH_D[1] = np.exp(int_deltak_D)
        gas.analysis.deltaH_D_err[1] = np.exp(int_deltak_D)*int_deltak_D_err

        return gas

    def langmuir_heat(gas):
        '''
        '''

        b = gas.b

        ln_b = np.log(b)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltab_model = sm.OLS(ln_b, inv_RT_with_const).fit()

        int_deltab, slope_deltab = deltab_model.params
        int_deltab_err, slope_deltab_err = deltab_model.bse

        # gas.analysis.deltaH_b = [-slope_deltab, np.exp(int_deltab)]
        # gas.analysis.deltaH_b_err = [slope_deltab_err,
        #   np.exp(int_deltab)*int_deltab_err]

        gas.analysis.deltaH_b[1] = np.exp(int_deltab)
        gas.analysis.deltaH_b_err[1] = np.exp(int_deltab)*int_deltab_err

        return gas

    if method == 'all':

        infinite_dilution_heat(gas)
        henry_heat(gas)
        langmuir_heat(gas)

    elif method == 'infinite_dilution':
        infinite_dilution_heat(gas)

    elif method == 'henry':
        henry_heat(gas)

    elif method == 'langmuir':
        langmuir_heat(gas)

    else:
        pyDMS.error_in_red(f'Unknown method {method} for heat of sorption')

    return gas
