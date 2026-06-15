"""
pyDMS.fugacity

Copyright 2026 Massachusetts Institute of Technology
Licensed under the 3-clause BSD license
"""

import numpy as np
import pyDMS

# Virial Coefficients tabulated from:
# Virial Coefficients of Pure Gases; Frenkel, M., Marsh, K. N., Eds.;
#   Landolt-Börnstein - Group IV Physical Chemistry;
#   Springer-Verlag: Berlin/Heidelberg, 2002; Vol. 21A.
#   https://doi.org/10.1007/b71692.

# the units for the virial coefficients are:
# B [=] cm^3/mol
# C [=] cm^6/mol^2

virial_coeff = {
    "He": {
        "B0": 9.2479,
        "B1": 1.0876 * 10**3,
        "B2": -1.088 * 10**5,
        "B3": 2.3869 * 10**6,
        "B4": 0,
        "C0": 0.10023 * 1e-3,
        "C1": 0,
        "C2": 0,
        "C3": 0,
        "C4": 0,
    },
    "H2": {
        "B0": 1.7472 * 10,
        "B1": -1.2926 * 10**2,
        "B2": -2.6988 * 10**5,
        "B3": 8.0282 * 10**6,
        "B4": 0,
        "C0": 0.53597 * 1e-3,
        "C1": 0,
        "C2": 0,
        "C3": 0,
        "C4": 0,
    },
    "N2": {
        "B0": 4.0286 * 10,
        "B1": -9.3378 * 10**3,
        "B2": -1.4164 * 10**6,
        "B3": 6.1253 * 10**7,
        "B4": -2.7198 * 10**9,
        "C0": 0.45178 * 1e-3,
        "C1": 282.49 * 1e-3,
        "C2": 0,
        "C3": 0,
        "C4": 0,
    },
    "O2": {
        "B0": 4.2859 * 10,
        "B1": -1.7696 * 10**4,
        "B2": 5.2007 * 10**5,
        "B3": -1.6393 * 10**8,
        "B4": 5.0855 * 10**9,
        "C0": 0.91432 * 1e-3,
        "C1": -57.003 * 1e-3,
        "C2": 38999 * 1e-3,
        "C3": 0,
        "C4": 0,
    },
    "CH4": {
        "B0": 4.4344 * 10,
        "B1": -1.6608 * 10**4,
        "B2": -3.5430 * 10**6,
        "B3": 2.9832 * 10**8,
        "B4": -2.3448 * 10**10,
        "C0": 1.468 * 1e-3,
        "C1": -417.68 * 1e-3,
        "C2": 2.1133e05 * 1e-3,
        "C3": 0,
        "C4": 0,
    },
    "CO2": {
        "B0": 5.7400 * 10,
        "B1": -3.8829 * 10**4,
        "B2": 4.2899 * 10**5,
        "B3": -1.4661 * 10**9,
        "B4": 0,
        "C0": 8.2273 * 1e-3,
        "C1": -11176 * 1e-3,
        "C2": 5.2971e06 * 1e-3,
        "C3": -6.7348e08 * 1e-3,
        "C4": 0,
    },
    "C2H6": {
        "B0": 0,  # TODO
        "B1": 0,  # TODO
        "B2": 0,  # TODO
        "B3": 0,  # TODO
        "B4": 0,  # TODO
        "C0": -21.966 * 1e-3,
        "C1": 19216 * 1e-3,
        "C2": -2.91e06 * 1e-3,
        "C3": 0,
        "C4": 0,
    },
    "C2H4": {
        "B0": 0,  # TODO
        "B1": 0,  # TODO
        "B2": 0,  # TODO
        "B3": 0,  # TODO
        "B4": 0,  # TODO
        "C0": -19.585 * 1e-3,
        "C1": 14199 * 1e-3,
        "C2": -1.879e06 * 1e-3,
        "C3": 0,
        "C4": 0,
    },
    "C3H8": {
        "B0": 1.0971 * 10**2,
        "B1": -8.4673 * 10**4,
        "B2": 8.1215 * 10**6,
        "B3": -3.4382 * 10**9,
        "B4": 0,
        "C0": 161.6 * 1e-3,
        "C1": -2.1173e05 * 1e-3,
        "C2": 9.5225e07 * 1e-3,
        "C3": -1.342e10 * 1e-3,
        "C4": 0,
    },
    "C3H6": {
        "B0": 1.0101 * 10**2,
        "B1": -7.5735 * 10**4,
        "B2": -7.9502 * 10**6,
        "B3": -2.7987 * 10**9,
        "B4": 0,
        "C0": -11.713 * 1e-3,
        "C1": 9511.1 * 1e-3,
        "C2": 0,
        "C3": 0,
        "C4": 0,
    },
}

# Peng-Robinson Parameters from NIST Chemistry Webbook
#   (https://webbook.nist.gov/chemistry/):
# Linstrom, P. J.; Mallard, W. G. The NIST Chemistry WebBook:
#   A Chemical Data Resource on the Internet.
#   J. Chem. Eng. Data 2001, 46 (5), 1059–1063.
#   https://doi.org/10.1021/je000236i.

pr_coeff = {"H2S": {"Tc": 373.1, "Pc": 9, "omega": 0.1}}  # MPa


def virial_eos(gas):
    """Computes the fugacity of a gas at a given temperature and pressure
            using a Virial expansion

    Built-in Virial coefficients are from: Virial Coefficients of Pure Gases;
        Frenkel, M., Marsh, K. N., Eds.;
        Landolt-Börnstein - Group IV Physical Chemistry;
        Springer-Verlag: Berlin/Heidelberg, 2002; Vol. 21A.
        https://doi.org/10.1007/b71692.


    Args:
        gas: An instance of the Gas class with
            Gas.formula, Gas.T, Gas.p populated
        params: A dictionary containing the Virial coefficients of
            a custom gas, for example:

            'params = {
            'B0':1*10**2,
            'B1':-7*10**4,
            'B2': -8*10**6,
            'B3': -3*10**9,
            'B4': 0,
            'C0': 0,
            'C1': 0}

    Returns:
        A numpy array in gas.f with the fugacities at the corresponding
            pressures provided in gas.p
    """

    gas_name = gas.formula

    if gas_name not in virial_coeff and not gas.virial_coeff:
        pyDMS.error_in_red(
            f"Error: No Virial coefficients found for {gas_name}. Please provide custom coefficients in Gas.virial_coeff or try peng_robinson_eos()"
        )

    temps = np.atleast_1d(gas.temp)  # ensure 1D array
    p_list = gas.p  # list of arrays

    R = 8.314e6  # cm3*Pa/(mol*K)
    atm_to_pa = 1.01325e5

    if gas.virial_coeff:
        virial_data = gas.virial_coeff
    else:
        virial_data = virial_coeff.get(gas_name, {})

    B0 = virial_data.get("B0", 0)
    B1 = virial_data.get("B1", 0)
    B2 = virial_data.get("B2", 0)
    B3 = virial_data.get("B3", 0)
    B4 = virial_data.get("B4", 0)
    C0 = virial_data.get("C0", 0)
    C1 = virial_data.get("C1", 0)
    C2 = virial_data.get("C2", 0)
    C3 = virial_data.get("C3", 0)
    C4 = virial_data.get("C4", 0)

    if len(temps) != len(p_list):
        pyDMS.error_in_red("gas.temp and gas.p must have the same number of entries")

    if B0 == 0:
        pyDMS.warning_in_orange("B0 = 0. This is unusual.")
    if B1 == 0:
        pyDMS.warning_in_orange("B1 = 0. This is unusual.")

    fugacities = []
    Z_rows = []
    for i, T in enumerate(temps):
        p_row = np.array(p_list[i]) * atm_to_pa  # convert to Pa
        B = B0 + B1 / T + B2 / T**2 + B3 / T**3 + B4 / T**4
        B_star = B / (R * T)
        C = C0 + C1 / T + C2 / T**2 + C3 / T**3 + C4 / T**4
        C_star = C / (R * T) ** 2 - B_star**2
        Z_val = 1 + B_star * p_row + C_star * p_row**2
        ln_phi = B_star * p_row + (1 / 2) * C_star * p_row**2
        phi = np.exp(ln_phi)
        f_pa = p_row * phi
        f_atm = f_pa / atm_to_pa
        Z_rows.append(Z_val)
        fugacities.append(f_atm)

    gas.Z = [np.array(row) for row in Z_rows]


    gas.f = [np.array(row) for row in fugacities]

    return gas


def peng_robinson_eos(gas):
    """Computes the fugacity of a gas at a given temperature and pressure
    using the Peng-Robinson EoS

    Built-in Virial coefficients are from the NIST Chemistry Webbook
    (https://webbook.nist.gov/chemistry/):
        Linstrom, P. J.; Mallard, W. G. The NIST Chemistry WebBook:
        A Chemical Data Resource on the Internet.
        J. Chem. Eng. Data 2001, 46 (5), 1059-1063.
        https://doi.org/10.1021/je000236i.


    Args:
        gas: An instance of the `Gas` class with `Gas.formula`, `Gas.T`, and `Gas.p` populated
        params: A dictionary containing the Peng-Robinson parameters of a custom gas, for example:
            'params = {
            'Tc':373.1,
            'Pc':9, #MPa
            'omega': 0.1}
    Returns: A numpy array in `gas.f` with the fugacities at the corresponding pressures provided
    in `gas.p`
    """
    gas_name = gas.formula

    if gas_name not in pr_coeff and not gas.pr_coeff:
        pyDMS.error_in_red(
            f"Error: No Peng-Robinson parameters found for {gas_name}. Please provide custom parameters in Gas.pr_coeff or try virial_eos()"
        )

    # ensure 1D array
    temps = np.atleast_1d(gas.temp)  # (n_temps,)
    p_list = gas.p  # list of arrays

    if gas.pr_coeff:
        pr_data = gas.pr_coeff
    else:
        pr_data = pr_coeff.get(gas_name, {})

    omega = pr_data.get("omega", 0)
    Tc = pr_data.get("Tc", 0)
    Pc = pr_data.get("Pc", 0)

    if omega == 0:
        pyDMS.warning_in_orange("WARNING: The accentric factor (omega) = 0")
    if Tc == 0:
        pyDMS.error_in_red("Error: Critical Temperature (Tc) cannot be 0")
    if Pc == 0:
        pyDMS.error_in_red("Error: Critical Pressure (Pc) cannot be 0")

    mpa_to_atm = 9.86923  # atm/MPa
    R = 8.314e-6  # MPa m^3 / mol-K

    fugacities = []
    Z_values = []

    for i, T in enumerate(temps):
        p_row = np.array(p_list[i]) / mpa_to_atm  # convert atm to MPa
        k = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
        alpha = (1 + k * (1 - np.sqrt(T / Tc))) ** 2
        a = 0.45724 * alpha * R**2 * Tc**2 / Pc
        b = 0.07780 * R * Tc / Pc

        A = a * p_row / (R**2 * T**2)
        B = b * p_row / (R * T)

        fug_row = []
        Z_row = []

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
            Z_row.append(zf)
            try:
                ln_phi = (
                    zf
                    - 1
                    - np.log(zf - Bij)
                    - Aij
                    / (2 * np.sqrt(2) * Bij)
                    * np.log((zf + (1 + np.sqrt(2)) * Bij) / (zf + (1 - np.sqrt(2)) * Bij))
                )
                phi = np.exp(ln_phi)
                fug_row.append(phi * p_val * mpa_to_atm)
            except (ValueError, ZeroDivisionError, FloatingPointError):
                fug_row.append(np.nan)

        Z_values.append(Z_row)
        fugacities.append(fug_row)

    gas.f = [np.array(row) for row in fugacities]

    gas.Z = [np.array(row) for row in Z_values]

    return gas


def virial_eos_z(gas, p, T):
    """Computes the compressibility factor of a gas at a given temperature and pressure
    using a Virial expansion

    Built-in Virial coefficients are from: Virial Coefficients of Pure Gases;
        Frenkel, M., Marsh, K. N., Eds.;
        Landolt-Börnstein - Group IV Physical Chemistry;
        Springer-Verlag: Berlin/Heidelberg, 2002; Vol. 21A.
        https://doi.org/10.1007/b71692.


    Args:
        gas: An instance of the `Gas` class with `Gas.formula` populated
        params: A dictionary containing the Virial coefficients of
            a custom gas, for example:
            'params = {
            'B0':1*10**2,
            'B1':-7*10**4,
            'B2': -8*10**6,
            'B3': -3*10**9,
            'B4': 0,
            'C0': 0,
            'C1': 0}

    Returns:
        A value for Z
    """

    gas_name = gas.formula

    if gas_name not in virial_coeff and not gas.virial_coeff:
        pyDMS.error_in_red(
            f"Error: No Virial coefficients found for {gas_name}. Please provide custom coefficients in Gas.virial_coeff or try peng_robinson_eos()"
        )

    R = 8.314e6  # cm3*Pa/(mol*K)
    atm_to_pa = 1.01325e5

    if gas.virial_coeff:
        virial_data = gas.virial_coeff
    else:
        virial_data = virial_coeff.get(gas_name, {})

    B0 = virial_data.get("B0", 0)
    B1 = virial_data.get("B1", 0)
    B2 = virial_data.get("B2", 0)
    B3 = virial_data.get("B3", 0)
    B4 = virial_data.get("B4", 0)
    C0 = virial_data.get("C0", 0)
    C1 = virial_data.get("C1", 0)
    C2 = virial_data.get("C2", 0)
    C3 = virial_data.get("C3", 0)
    C4 = virial_data.get("C4", 0)

    if B0 == 0:
        pyDMS.warning_in_orange("B0 = 0. This is unusual.")
    if B1 == 0:
        pyDMS.warning_in_orange("B1 = 0. This is unusual.")

    p_row = p * atm_to_pa  # convert to Pa
    B = B0 + B1 / T + B2 / T**2 + B3 / T**3 + B4 / T**4
    B_star = B / (R * T)
    C = C0 + C1 / T + C2 / T**2 + C3 / T**3 + C4 / T**4
    C_star = C / (R * T) ** 2 - B_star**2
    Z_val = 1 + B_star * p_row + C_star * p_row**2

    return Z_val


def peng_robinson_eos_z(gas, p, T):
    """Computes the compressibility factor of a gas at a given temperature and pressure
    using the Peng-Robinson EoS

    Built-in Virial coefficients are from the NIST Chemistry Webbook
        (https://webbook.nist.gov/chemistry/):
        Linstrom, P. J.; Mallard, W. G. The NIST Chemistry WebBook:
        A Chemical Data Resource on the Internet.
        J. Chem. Eng. Data 2001, 46 (5), 1059-1063.
        https://doi.org/10.1021/je000236i.


    Args:
        gas: An instance of the `Gas` class with `Gas.formula` populated
        params: A dictionary containing the Peng-Robinson parameters of a custom gas, for example:
            'params = {
            'Tc':373.1,
            'Pc':9, #MPa
            'omega': 0.1}

    Returns:
        A value for Z
    """
    gas_name = gas.formula

    if gas_name not in pr_coeff and not gas.pr_coeff:
        pyDMS.error_in_red(
            f"Error: No Peng-Robinson parameters found for {gas_name}. Please provide custom parameters in Gas.pr_coeff or try virial_eos()"
        )

    # ensure 1D array

    if gas.pr_coeff:
        pr_data = gas.pr_coeff
    else:
        pr_data = pr_coeff.get(gas_name, {})

    omega = pr_data.get("omega", 0)
    Tc = pr_data.get("Tc", 0)
    Pc = pr_data.get("Pc", 0)

    if omega == 0:
        pyDMS.warning_in_orange("WARNING: The accentric factor (omega) = 0")
    if Tc == 0:
        pyDMS.error_in_red("Error: Critical Temperature (Tc) cannot be 0")
    if Pc == 0:
        pyDMS.error_in_red("Error: Critical Pressure (Pc) cannot be 0")

    mpa_to_atm = 9.86923  # atm/MPa
    R = 8.314e-6  # MPa m^3 / mol-K

    p_row = p / mpa_to_atm  # convert atm to MPa
    k = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + k * (1 - np.sqrt(T / Tc))) ** 2
    a = 0.45724 * alpha * R**2 * Tc**2 / Pc
    b = 0.07780 * R * Tc / Pc

    p_row = np.array([p_row])

    A = a * p_row / (R**2 * T**2)
    B = b * p_row / (R * T)


    for j, p_val in enumerate(p_row):
        Aij = A[j]
        Bij = B[j]
        coeffs = [1, -(1 - Bij), Aij - 3 * Bij**2 - 2 * Bij, -Aij * Bij + Bij**2 + Bij**3]
        z_roots = np.roots(coeffs)
        z_real = np.real(z_roots[np.isreal(z_roots)])

        zf = z_real[np.argmax(z_real)]  # use largest Z (vapor phase)

    return zf
