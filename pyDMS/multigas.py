'''
pyDMS.multigas

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import numpy as np

import pyDMS


def mixed_isotherm(*gases, p_or_f, mol_frac, temp):
    '''
    Computes the mixed gas DMS isotherm
    p: total pressure to compute DMS at
    params: DMS parameters in a dictionary format
    x: mole fractions of mixture
    eos: Compute fugacity from provided p if desired using a virial approach
        or Peng-Robinson EOS
    T: temperature to compute fugacity at
    gases: gases to compute fugacities of

    Returns: [partial pressure/fugacity, concentration, concentration error]
    '''

    if sum(mol_frac) != 1:
        pyDMS.warning_in_orange(f'Sum of mole fractions is {sum(mol_frac)}')

    # p_f_max_array = np.zeros(len(mol_frac))

    # *currently this just overwrites which is fine
    # if every gas has the same location but we cant ensure that
    for i, gas in enumerate(gases):
        index = np.where(gas.temp == temp)[0][0]
    # this is somewhat useless but fine for now
    p_use = np.ones((len(gases), len(p_or_f)))*p_or_f
    p_partial = np.zeros_like(p_use)

    for i, val in enumerate(mol_frac):
        p_partial[i] = p_or_f*val

    p0 = p_partial[0]

    b_terms_sum = 0  # initializing b1*p1+b2*p2+...+bi*pi
    term3_err_prop = 0  # initializing the error propogation
    term2_partial_err_prop = 0  # initializing the error propogation

    for i, gas in enumerate(gases):

        if i == 0:
            ch = gas.CH[index]
            kd = gas.kD[index]
            b0 = gas.b[index]
            ch0_err = gas.CH_err[index]
            kd_err = gas.kD_err[index]
            b0_err = gas.b_err[index]

        b_terms_sum += gas.b[index] * p_partial[i]

    b_terms_sum_plus_one = 1 + b_terms_sum

    # error propogation
    for i, gas in enumerate(gases):
        if i != 0:
            b = gas.b[index]
            b_err = gas.b_err[index]
            term3_err_prop += (
                b0**2*ch**2*p0**2*p_partial[i]**2*b_err**2
                / (b_terms_sum_plus_one)**4
                )

            term2_partial_err_prop += b*p_partial[i]

    term1_err_prop = p0**2*kd_err**2
    term2_err_prop = (
        ch**2*p0**2*(1+term2_partial_err_prop)**2*b0_err**2
        / (b_terms_sum_plus_one)**4
        )
    term4_err_prop = b0**2*p0**2*ch0_err**2/(b_terms_sum_plus_one)**2
    c = kd*p0+ch*b0*p0/(b_terms_sum_plus_one)
    c_err = np.sqrt(
        term1_err_prop+term2_err_prop+term3_err_prop+term4_err_prop
        )

    return [p0, c, c_err]


def selectivity(*isotherms, calc=None):
    '''
    # calc options: 1_numerator, 2_numerator*
    '''

    if len(isotherms) not in [2, 3]:
        pyDMS.error_in_red("selectivity() requires 2 or 3 gas objects.")

    if len(isotherms) == 3 and calc is None:
        pyDMS.error_in_red('with three gases, calc must be specified '
                           'as 1_numerator or 2_numerator')

    iso1, iso2, *rest = isotherms
    iso3 = rest[0] if rest else None

    p1 = iso1[0]
    c1 = iso1[1]
    c1err = iso1[2]
    p2 = iso2[0]
    c2 = iso2[1]
    c2err = iso2[2]

    p3 = 0

    if iso3:
        p3 = iso3[0]
        c3 = iso3[1]
        c3err = iso3[2]

    p_tot = p1+p2+p3

    if len(isotherms) == 2:
        alpha = (c1/p1)/(c2/p2)
        alpha_err = np.sqrt(
            p2**2*(c2**2*c1err**2+c1**2*c2err**2)/(c2**4*p1**2)
            )

    elif calc == '2_numerator':
        alpha = ((c1/p1)+(c2/p2))/(c3*p3)
        alpha_err = np.zeros(len(c1))

    elif calc == '1_numerator':
        alpha = (c1/p1)/((c2/p2)+(c3*p3))

        alpha_err_term1 = p2**2*p3**2*c1err**2/(c3*p1*p2+c2*p1*p3)**2
        alpha_err_term2 = c1**2*p2**2*p3**4*c2err**2/(p1**2*(c3*p2+c2*p3)**4)
        alpha_err_term3 = c1**1*p2**4*p3**2*c3err**2/(p1**2*(c3*p2+c2*p3)**4)
        alpha_err = np.sqrt(alpha_err_term1+alpha_err_term2+alpha_err_term3)

    else:
        pyDMS.error_in_red('calc="1_numerator" or "2_numerator" required '
                           'for ternary selectivity')

    return [p_tot, alpha, alpha_err]
