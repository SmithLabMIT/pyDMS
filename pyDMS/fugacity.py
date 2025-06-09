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

# the units for the virial coefficients are:
# B [=] cm^3/mol
# C [=] cm^6/mol^2

virial_coeff = {
        'He':{
            'B0': 9.2479, # check
            'B1': 1.0876*10**3, # check
            'B2': 1.088*10**5, # check
            'B3': 2.3869*10**6, # check
            'B4': 0, # check
            'C0': 0.10023*1E-3,
            'C1': 0,
            'C2': 0,
            'C3': 0,
            'C4': 0},
        'H2':{  
            'B0':1.7472*10, # check
            'B1':-1.2926*10**2, # check
            'B2':-2.6988*10**5, # check
            'B3':8.0282*10**6, # check
            'B4': 0, # check
            'C0': 0.53597*1E-3,
            'C1': 0,
            'C2': 0,
            'C3': 0,
            'C4': 0},
        'N2':{
            'B0':4.0286*10, # check
            'B1':-9.3378*10**3, # check
            'B2': -1.4164*10**6, # check
            'B3': 6.1253*10**7, # check
            'B4': -2.7198*10**9, # check
            'C0': 0.45178*1E-3,
            'C1': 282.49*1E-3,
            'C2': 0,
            'C3': 0,
            'C4': 0},
        'O2':{
            'B0':4.2859*10, # check
            'B1':-1.7696*10**4, # check
            'B2': 5.2007*10**5, # check
            'B3': -1.6393*10**8, # check
            'B4': 5.0855*10**9, # check
            'C0': 0.91432*1E-3,
            'C1': -57.003*1E-3,
            'C2': 38999*1E-3,
            'C3': 0,
            'C4': 0},
        'CH4':{
            'B0':4.4344*10, # check
            'B1': -1.6608*10**4, # check
            'B2': -3.5430*10**6, # check
            'B3': 2.9832*10**8, # check
            'B4': -2.3448*10**10, # check
            'C0': 1.468*1E-3,
            'C1': -417.68*1E-3,
            'C2': 2.1133e+05*1E-3,
            'C3': 0,
            'C4': 0},
        'CO2':{
            'B0':5.74*10, # check
            'B1':-3.8829*10**4, # check
            'B2': 4.2899*10**5, # check
            'B3': -1.4661*10**9, # check
            'B4': 0, # check
            'C0': 8.2273*1E-3,
            'C1': -11176*1E-3,
            'C2': 5.2971e+06*1E-3,
            'C3': -6.7348e+08*1E-3,
            'C4': 0},
        'C2H6':{
            'B0':0, # check
            'B1':0, # check
            'B2': 0, # check
            'B3': 0, # check
            'B4': 0, # check
            'C0': -21.966*1E-3,
            'C1': 19216*1E-3,
            'C2': -2.91e+06*1E-3,
            'C3': 0,
            'C4': 0},
        'C2H4':{
            'B0':0, # check
            'B1':0, # check
            'B2': 0, # check
            'B3': 0, # check
            'B4': 0, # check
            'C0': -19.585*1E-3,
            'C1': 14199*1E-3,
            'C2': -1.879e+06*1E-3,
            'C3': 0,
            'C4': 0},
        'C3H8':{
            'B0':1.0971*10**2, # check
            'B1':-8.4673*10**4, # check
            'B2': 8.1215*10**6, # check
            'B3': -3.4382*10**9, # check
            'B4': 0,
            'C0': 161.6*1E-3,
            'C1': -2.1173e+05*1E-3,
            'C2': 9.5225e+07*1E-3,
            'C3': -1.342e+10*1E-3,
            'C4': 0},
        'C3H6':{
            'B0':1.0101*10**2, # check
            'B1':-7.5735*10**4, # check
            'B2': -7.9502*10**6, # check
            'B3': -2.7987*10**9, # check
            'B4': 0, # check
            'C0': -11.713*1E-3,
            'C1': 9511.1*1E-3,
            'C2': 0,
            'C3': 0,
            'C4': 0}}

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

def virial_eos(gas):
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
    temps = np.atleast_1d(gas.temp)  # ensure 1D array
    p_grid = np.atleast_2d(gas.p)    # ensure 2D array

    R = 8.314e6  # cm³·Pa/(mol·K)
    atm_to_pa = 1.01325e5
    p_grid_pa = p_grid * atm_to_pa

    if gas.virial_coeff:
        virial_data = gas.virial_coeff
    else:
        virial_data = virial_coeff.get(gas_name, {})

    B0 = virial_data.get('B0', 0)
    B1 = virial_data.get('B1', 0)
    B2 = virial_data.get('B2', 0)
    B3 = virial_data.get('B3', 0)
    B4 = virial_data.get('B4', 0)
    C0 = virial_data.get('C0', 0)
    C1 = virial_data.get('C1', 0)
    C2 = virial_data.get('C2', 0)
    C3 = virial_data.get('C3', 0)
    C4 = virial_data.get('C4', 0)

    if B0 == 0:
        pyDMS.warning_in_orange("B0 = 0. This is unusual.")
    if B1 == 0:
        pyDMS.warning_in_orange("B1 = 0. This is unusual.")

    fugacities = []

    for i, T in enumerate(temps):
        p_row = p_grid_pa[i]
        B = B0 + B1/T + B2/T**2 + B3/T**3 + B4/T**4
        C = C0 + C1/T + C2/T**2 + C3/T**3 + C4/T**4
        vm = R * T / p_row
        ln_phi = B/vm + (C + B**2) / (2 * vm**2)
        phi = np.exp(ln_phi)
        f_pa = p_row * phi
        f_atm = f_pa / atm_to_pa
        fugacities.append(f_atm)

    gas.f = np.array(fugacities)  # shape: (n_conditions, n_pressures)
    return gas

def peng_robinson_eos(gas):
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
    temps = np.atleast_1d(gas.temp)   # (n_temps,)
    p_grid = np.atleast_2d(gas.p)     # (n_temps, n_pressures)

    if gas.pr_coeff:
        pr_data = gas.pr_coeff
    else:
        pr_data = pr_coeff.get(gas_name, {})

    omega = pr_data.get('omega', 0)
    Tc = pr_data.get('Tc', 0)
    Pc = pr_data.get('Pc', 0)

    if omega == 0:
        pyDMS.warning_in_orange("WARNING: The accentric factor (omega) = 0")
    if Tc == 0:
        pyDMS.error_in_red("Error: Critical Temperature (Tc) cannot be 0")
    if Pc == 0:
        pyDMS.error_in_red("Error: Critical Pressure (Pc) cannot be 0")

    mpa_to_atm = 9.86923  # atm/MPa
    R = 8.314  # J/(mol·K)

    fugacities = []

    for i, T in enumerate(temps):
        p_row = p_grid[i] / mpa_to_atm  # convert atm → MPa for this row
        k = 0.375 + 1.542 * omega - 0.270 * omega**2
        alpha = (1 + k * (1 - np.sqrt(T / Tc)))**2
        a = 0.457 * alpha * R**2 * Tc**2 / Pc
        b = 0.0778 * R * Tc / Pc

        A = a * p_row / (R**2 * T**2)
        B = b * p_row / (R * T)

        fug_row = []

        for j, p_val in enumerate(p_row):
            Aij = A[j]
            Bij = B[j]
            coeffs = [1, -(1 - Bij), Aij - 3 * Bij**2 - 2 * Bij, -Aij * Bij + Bij**2 + Bij**3]
            z_roots = np.roots(coeffs)
            z_real = np.real(z_roots[np.isreal(z_roots)])

            if len(z_real) == 0:
                fug_row.append(np.nan)
                continue

            zf = z_real[np.argmax(z_real)]  # use largest Z (vapor phase)

            try:
                ln_phi = (
                    zf - 1 - np.log(zf - Bij)
                    - Aij / (2 * np.sqrt(2) * Bij)
                    * np.log((zf + (1 + np.sqrt(2)) * Bij) / (zf + (1 - np.sqrt(2)) * Bij))
                )
                phi = np.exp(ln_phi)
                fug_row.append(phi * p_val * mpa_to_atm)
            except (ValueError, ZeroDivisionError, FloatingPointError):
                fug_row.append(np.nan)

        fugacities.append(fug_row)

    gas.f = np.array(fugacities)
    return gas