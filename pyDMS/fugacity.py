r'''
pyDMS.fugacity

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import numpy as np
import pyDMS

# Virial Coefficients tabulated from:
# Virial Coefficients of Pure Gases; Frenkel, M., Marsh, K. N., Eds.;
#   Landolt-Börnstein - Group IV Physical Chemistry;
#   Springer-Verlag: Berlin/Heidelberg, 2002; Vol. 21A. https://doi.org/10.1007/b71692.

virial_coeff = {
        'He':{
            'B0':9.2479,
            'B1':1.0876*10**3,
            'B2': 1.088*10**5,
            'B3': 2.3869*10**6,
            'B4': 0,
            'C0': 0.13876333*1000,
            'C1': -0.000131193*1000},
        'H2':{  
            'B0':1.7472*10,
            'B1':-1.2926*10**2,
            'B2':-2.6988*10**5,
            'B3':8.0282*10**6,
            'B4': 0,
            'C0': -0.0893351*1000,
            'C1': 0.001954*1000},
        'N2':{
            'B0':4.0286*10,
            'B1':-9.3378*10**3,
            'B2': -1.4164*10**6,
            'B3': 6.1253*10**7,
            'B4': -2.7198*10**9,
            'C0': 2.433098621*1000,
            'C1': -0.003424843*1000},
        'O2':{
            'B0':4.2859*10,
            'B1':-1.7696*10**4,
            'B2': 5.2007*10**5,
            'B3': -1.6393*10**8,
            'B4': 5.0855*10**9,
            'C0': 1.199412183*1000,
            'C1': -0.000178252*1000},
        'CH4':{
            'B0':4.4344*10,
            'B1': -1.6608*10**4,
            'B2': -3.5430*10**6,
            'B3': 2.9832*10**8,
            'B4': -2.3448*10**10,
            'C0': 4.930564307*1000,
            'C1': -0.008472446*1000},
        'CO2':{
            'B0':5.74*10,
            'B1':-3.8829*10**4,
            'B2': 4.2899*10**5,
            'B3': -1.4661*10**9,
            'B4': 0,
            'C0': 12.53856029*1000,
            'C1': -0.025514329*1000},
        'C2H6':{
            'B0':None,
            'B1':None,
            'B2': None,
            'B3': None,
            'B4': None,
            'C0': None,
            'C1': None},
        'C2H4':{
            'B0':None,
            'B1':None,
            'B2': None,
            'B3': None,
            'B4': None,
            'C0': None,
            'C1': None},
        'C3H8':{
            'B0':1.0971*10**2,
            'B1':-8.4673*10**4,
            'B2': 8.1215*10**6,
            'B3': -3.4382*10**9,
            'B4': 0,
            'C0': 26283.52145,
            'C1': 0},
        'C3H6':{
            'B0':1.0101*10**2,
            'B1':-7.5735*10**4,
            'B2': -7.9502*10**6,
            'B3': -2.7987*10**9,
            'B4': 0,
            'C0': 0,
            'C1': 0}}

# Peng-Robinson Parameters from NIST Chemistry Webbook (https://webbook.nist.gov/chemistry/):
# Linstrom, P. J.; Mallard, W. G. The NIST Chemistry WebBook:
#   A Chemical Data Resource on the Internet.
#   J. Chem. Eng. Data 2001, 46 (5), 1059–1063. https://doi.org/10.1021/je000236i.

pr_coeff = {
    'H2S': {
        'Tc':373.1, 
        'Pc':9, #MPa
        'omega': 0.1}
}

def virial_eos(gas, params=None):
    '''Computes the fugacity of a gas at a given temperature and pressure using a Virial expansion

    Built-in Virial coefficients are from: Virial Coefficients of Pure Gases; Frenkel, M., Marsh, K. N., Eds.;
    Landolt-Börnstein - Group IV Physical Chemistry;
    Springer-Verlag: Berlin/Heidelberg, 2002; Vol. 21A. https://doi.org/10.1007/b71692.

        
    Args:
        gas: An instance of the Gas class with Gas.formula, Gas.T, Gas.p populated
        params: A dictionary containing the Virial coefficients of a custom gas, for example:
            'params = {
            'B0':1*10**2,
            'B1':-7*10**4,
            'B2': -8*10**6,
            'B3': -3*10**9,
            'B4': 0,
            'C0': 0,
            'C1': 0}

    Returns:
        A numpy array in gas.f with the fugacities at the corresponding pressures provided in gas.p
    '''

    gas_name = gas.formula
    T = gas.T
    p = gas.p

    R = 8.314e6  # cm^3·Pa/(mol·K)
    atm_to_pa = 1.01325e5  # Conversion factor from atm to Pa
    p = p * atm_to_pa

    if params is not None:
        B0 = params.get('B0',0)
        B1 = params.get('B1',0)
        B2 = params.get('B2',0)
        B3 = params.get('B3',0)
        B4 = params.get('B4',0)
        C0 = params.get('C0',0)
        C1 = params.get('C1',0)
    else:
        B0 = gas_name.get('B0',0)
        B1 = gas_name.get('B1',0)
        B2 = gas_name.get('B2',0)
        B3 = gas_name.get('B3',0)
        B4 = gas_name.get('B4',0)
        C0 = gas_name.get('C0',0)
        C1 = gas_name.get('C1',0)

    if B0==0:
        pyDMS.warning_in_orange("WARNING: B0 = 0")
    if B1==0:
        pyDMS.warning_in_orange("WARNING: B1 = 0")

    B =  B0 + B1/T + B2/T**2 + B3/T**3 + B4/T**4
    C = C0+C1*T
    vm = R*T/p
    ln_phi = B/vm + (C+B**2)/(2*vm**2)
    phi = np.exp(ln_phi)
    f = p * phi
    f_atm = f / atm_to_pa

    gas.f = np.array(f_atm) 
    
    return gas

def peng_robinson_eos(gas, params=None):
    '''Computes the fugacity of a gas at a given temperature and pressure using the Peng-Robinson EoS

    Built-in Virial coefficients are from the NIST Chemistry Webbook (https://webbook.nist.gov/chemistry/):
        Linstrom, P. J.; Mallard, W. G. The NIST Chemistry WebBook:
        A Chemical Data Resource on the Internet.
        J. Chem. Eng. Data 2001, 46 (5), 1059–1063. https://doi.org/10.1021/je000236i.

        
    Args:
        gas: An instance of the Gas class with Gas.formula, Gas.T, Gas.p populated
        params: A dictionary containing the Peng-Robinson parameters of a custom gas, for example:
            'params = {
            'Tc':373.1, 
            'Pc':9, #MPa
            'omega': 0.1}

    Returns:
        A numpy array in gas.f with the fugacities at the corresponding pressures provided in gas.p
    '''
    gas_name = gas.formula
    T = gas.T
    p = gas.p
    
    if params is not None:
        omega = params.get('omega',0)
        Tc = params.get('Tc',0)
        pc = params.get('Pc',0)
    else:
        omega = gas_name.get('omega',0)
        Tc = gas_name.get('Tc',0)
        pc = gas_name.get('Pc',0)

    if omega==0:
        pyDMS.warning_in_orange("WARNING: The accentric factor (omega) = 0")

    if Tc==0:
        pyDMS.error_in_red("Error: Critical Temperature (Tc) cannot be 0")

    if Pc==0:
        pyDMS.error_in_red("Error: Critical Pressure (Pc) cannot be 0")

    mpa_to_atm = 9.86923 # atm/MPa
    p_mpa = p/mpa_to_atm
    R = 8.314 # J/(mol·K)

    k = 0.375 + 1.542*omega - 0.270*omega**2
    a = 0.457*(1+k*(1-np.sqrt(T/Tc)))**2*R**2*Tc**2/pc
    b = 0.0778*R*Tc/pc
    A = a * p_mpa / (R**2 * T**2)
    B = b * p_mpa / (R * T)

    f_atm_list = []  # Store fugacity values for each pressure
    
    for i, p_val in enumerate(p):
        coeffs = [1, -(1-B[i]), A[i]-3*B[i]**2-2*B[i], -A[i]*B[i]+B[i]**2+B[i]**3]
        z_factors = np.roots(coeffs)  # Solve for Z
        real_roots = np.real(z_factors[np.isreal(z_factors)])  # Keep only real solutions

        phi_array = []

        for zf in real_roots:
            if zf > B[i]:  # Ensure the log argument is valid
                phi = np.exp(
                    zf - 1 - np.log(zf - B[i]) - A[i] / (2 * np.sqrt(2) * B[i]) *
                    np.log((zf + (1 + np.sqrt(2)) * B[i]) / (zf + (1 - np.sqrt(2)) * B[i]))
                )
                phi_array.append(phi)

        if phi_array:
            phi_use = np.min(phi_array)  # Ensure proper array handling
            f_atm = phi_use * p_val
            f_atm_list.append(f_atm)
        else:
            f_atm_list.append(np.nan)  # If no valid phi, store NaN
    gas.f = np.array(f_atm_list)  # Return an array of fugacities
    
    return gas
