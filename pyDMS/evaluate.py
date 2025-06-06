r'''
pyDMS.evaluate            

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import numpy as np
import statsmodels.api as sm
import pyDMS
from . import fugacity

def isotherm(gas, index, eos=None):
    '''
    Computes the DMS isotherm
    p: pressure/fugacity (atm) to compute DMS at
    params: DMS parameters in a dictionary format
    eos: Compute fugacity from provided p if desired using a virial approach or Peng-Robinson EOS
    T: temperature to compute fugacity at
    gas: gas to compute fugacity of

    Returns: [pressure/fugacity, concentration, concentration error]
    '''

    # * update parameter locations
    ch = gas.CH[index]
    kd = gas.kD[index]
    b = gas.b[index]
    ch_err = gas.CH_err[index]
    kd_err = gas.kD_err[index]
    b_err = gas.b_err[index]

    p = np.linspace(1E-6,np.max(gas.p[index]),100)

    # making pylint happy
    p_val = None

    if eos=='virial':
        p_val = fugacity.virial_eos(gas)
    elif eos=='pr':
        p_val = fugacity.peng_robinson_eos(gas)
    elif eos is None:
        p_val = p
    else:
        print("eos not recognized")

    c = kd*p_val+ch*b*p_val/(1+b*p_val)
    c_err = np.sqrt(ch**2*p_val**2*b_err**2/(1+b*p_val)**4+b**2*p_val**2*ch_err**2/(1+b*p_val)**2+p_val**2*kd_err**2)

    return p_val, c, c_err

def isosteric_heat(gas, calc_fug=False):
    
    #name = gas.name
    c_vec = gas.c
    cerr_vec = gas.c_err
    p_vec = gas.p
    T_vec = gas.temp

    #if calc_fug:
        # calculate fugacity here

    # extract dms parameters
    kd_f = gas.kD
    b_f = gas.b
    ch_vec_f = gas.CH


    # calculate inverse temperatures
    T_inv = 1/T_vec

    # determine how many points to use
    #def minimize_dms(gas, C_target):

    # create an array of target concentrations (make it more detailed in beginning)
    # create an array of pressure guesses

    # create a function to minimize

    # use fmincon to minimize to find pressure

    # fit inverse T vs log(p_calc)

    # extract slope

    # find the average pressure

    # z-value of the average pressure

    # calculate the isosteric heat
    #return

def S_inf(gas):
    '''
    *
    '''
    
    k_D = gas.kD
    b = gas.b
    ch = gas.CH
    S_inf_calc = k_D+ch*b

    gas.analysis.S_inf = S_inf_calc
    gas.analysis.S_inf_err = 'Not yet implemented*'

    return gas

def heat_of_sorption(gas, method='all'):
    '''
    *
    '''
    
    temp = gas.temp

    inv_RT = 1/(0.008314*temp)

    def infinite_dilution_heat(gas):
        '''

        '''
        S_inf(gas)
    
        S_inf_calc = gas.analysis.S_inf
        S_inf_err_calc = gas.analysis.S_inf_err
        #print(S_inf_calc)
        ln_S_inf = np.log(S_inf_calc)

        inv_RT_with_const = sm.add_constant(inv_RT)
        #print(ln_S_inf)
        deltaH_S_inf_model = sm.OLS(ln_S_inf, inv_RT_with_const).fit()
        #print(deltaH_S_inf_model.params)
        int_deltaH_S_inf, slope_deltaH_S_inf = deltaH_S_inf_model.params
        int_deltaH_S_inf_err, slope_deltaH_S_inf_err = deltaH_S_inf_model.bse
        #print(int_deltaH_S_inf, slope_deltaH_S_inf)
        gas.analysis.deltaH_S_inf = [-slope_deltaH_S_inf, np.exp(int_deltaH_S_inf)]
        gas.analysis.deltaH_S_inf_err = [slope_deltaH_S_inf_err,np.exp(int_deltaH_S_inf)*int_deltaH_S_inf_err]

        return gas


    def henry_heat(gas):
        '''
        '''

        k_D = gas.kD
        
        ln_k_D = np.log(k_D)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltak_D_model = sm.OLS(ln_k_D,inv_RT_with_const).fit()

        int_deltak_D, slope_deltak_D = deltak_D_model.params
        int_deltak_D_err, slope_deltak_D_err = deltak_D_model.bse

        #gas.analysis.deltaH_D = [-slope_deltak_D, np.exp(int_deltak_D)]
        #gas.analysis.deltaH_D_err = [slope_deltak_D_err,np.exp(int_deltak_D)*int_deltak_D_err]

        gas.analysis.deltaH_D[1] = np.exp(int_deltak_D)
        gas.analysis.deltaH_D_err[1] = np.exp(int_deltak_D)*int_deltak_D_err

        return gas

    def langmuir_heat(gas):
        '''
        '''
        
        b = gas.b

        ln_b = np.log(b)

        inv_RT_with_const = sm.add_constant(inv_RT)

        deltab_model = sm.OLS(ln_b,inv_RT_with_const).fit()

        int_deltab, slope_deltab = deltab_model.params
        int_deltab_err, slope_deltab_err = deltab_model.bse
        
        #gas.analysis.deltaH_b = [-slope_deltab, np.exp(int_deltab)]
        #gas.analysis.deltaH_b_err = [slope_deltab_err, np.exp(int_deltab)*int_deltab_err]

        gas.analysis.deltaH_b[1] = np.exp(int_deltab)
        gas.analysis.deltaH_b_err[1] = np.exp(int_deltab)*int_deltab_err

        return gas

    if method=='all':

        infinite_dilution_heat(gas)
        henry_heat(gas)
        langmuir_heat(gas)

    elif method=='infinite_dilution':
        infinite_dilution_heat(gas)

    elif method=='henry':
        henry_heat(gas)

    elif method=='langmuir':
        langmuir_heat(gas)

    else:
        pyDMS.error_in_red(f'Unknown method {method} for heat of sorption')

    return gas