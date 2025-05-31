'''
pyDMS.multigas
                          
Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import pyDMS


def mixed_isotherm(p, params, x, eos=None, T=None, gases=None):
    '''
    Computes the mixed gas DMS isotherm
    p: total pressure to compute DMS at
    params: DMS parameters in a dictionary format
    x: mole fractions of mixture
    eos: Compute fugacity from provided p if desired using a virial approach or Peng-Robinson EOS
    T: temperature to compute fugacity at
    gases: gases to compute fugacities of

    Returns: [partial pressure/fugacity, concentration, concentration error]
    '''

    # * update parameter getting

    if eos:
        if len(params) != len(x) != len(eos) != len(T)  != len(gases):
            pyDMS.error_in_red('Dimensions of params, x, eos, T, gases must match')
    else:
        if len(params) != len(x):
            pyDMS.error_in_red(f'dim(params)={len(params)} != dim(x)={len(x)}')
        
    if sum(x) != 1:
        pyDMS.warning_in_orange(f'Sum of mole fractions is {sum(x)}')

    
    p_val_partial = np.zeros((len(x),len(p))) # initializing partial pressures

    b_terms_sum = 0 # initializing b1*p1+b2*p2+...+bi*pi
    term3_err_prop = 0 # initializing the error propogation
    term2_partial_err_prop = 0 # initializing the error propogation

    params_dict = {f"params_{i+1}": param for i, param in enumerate(params)} # extracting params

    # computing partial pressures/fugacities
    for i, x_val in enumerate(x):

        if eos:
            eos_ind = eos[i]
        else:
            eos_ind = None

        if eos_ind=='virial':
            p_val_partial[i] = virial_eos(gases[i],T[i],p)*x_val
        elif eos_ind=='pr':
            p_val_partial[i] = peng_robinson_eos(gases[i],T[i],p)*x_val
        elif eos_ind is None:
            p_val_partial[i] = p*x_val
        else:
            dms.error_in_red('Equation of State not recognized')

    p0 = p_val_partial[0]

    # initializing to make pylint happy
    ch = None
    kd = None
    b0 = None
    ch0_err = None
    kd_err = None
    b0_err = None

    for gas, param in enumerate(params):
        if gas == 0:
            ch = param["Ch"]
            kd = param["kd"]
            b0 = param["b"]
            kd_err = param.get("kd_err", 0)
            b0_err = param.get("b_err", 0)
            ch0_err = param.get("Ch_err", 0)

        b_terms_sum += param["b"] * p_val_partial[gas]  # denominator of the mixed gas DMS

    b_terms_sum_plus_one = 1 + b_terms_sum

    # error propogation
    for gas, param in enumerate(params):
        if gas != 0:
            b = param.get("b")
            b_err = param.get("b_err")
            term3_err_prop += b0**2*ch**2*p0**2*p_val_partial[gas]**2*b_err**2/(b_terms_sum_plus_one)**4

            term2_partial_err_prop += b*p_val_partial[gas]

    term1_err_prop = p0**2*kd_err**2
    term2_err_prop = ch**2*p0**2*(1+term2_partial_err_prop)**2*b0_err**2/(b_terms_sum_plus_one)**4
    term4_err_prop = b0**2*p0**2*ch0_err**2/(b_terms_sum_plus_one)**2
    c = kd*p_val_partial[0]+ch*b0*p_val_partial[0]/(b_terms_sum_plus_one)
    c_err =  np.sqrt(term1_err_prop + term2_err_prop + term3_err_prop + term4_err_prop)
    return [p0, c, c_err]

def selectivity(params1, params2, params3=None, calc='classic', custom=None, custom_err=None):
    '''
    paramNone1: [fug/press, concentration]
    params2: [Noneug/press, concentration]
    *args: additiNonenal params as desired
    custom: custom sNonelectivity rule (e.g., if you have multiple acid gases)
    
    It is important to note that the fugacity/pressure provided should be the partial pressure/fugacity of that gas
    The length of each fug/press vector must be the same as selectivity will read over and assume each index corresponds to the same total pressure 
    '''

    # initializing to make pylint happy
    alpha = 0
    alpha_err = 0
    p1 = params1[0]
    c1 = params1[1]
    c1err = params1[2]
    p2 = params2[0]
    c2 = params2[1]
    c2err = params2[2]
    p3 = 0
    c3 = 0
    c3err = 0

    if params3:
        p3 = params3[0]
        c3 = params3[1]
        c3err = params3[2]
        

    # this is here for when custom selectivities are implemented
    #for idx, arg in enumerate(args):        
    # extract the required columns from each arg (arg[1])

    if calc == 'classic':
        alpha = (c1/p1)/(c2/p2)
        alpha_err = np.sqrt(p2**2*(c2**2*c1err**2+c1**2*c2err**2)/(c2**4*p1**2))
    elif calc == 'sum_numer':
        alpha = ((c1/p1)+(c2/p2))/(c3*p3)
        alpha_err = np.zeros(len(c1))
        
    elif calc == 'sum_denom':
        alpha = (c1/p1)/((c2/p2)+(c3*p3))
        
        alpha_err_term1 = p2**2*p3**2*c1err**2/(c3*p1*p2+c2*p1*p3)**2
        alpha_err_term2 = c1**2*p2**2*p3**4*c2err**2/(p1**2*(c3*p2+c2*p3)**4)
        alpha_err_term3 = c1**1*p2**4*p3**2*c3err**2/(p1**2*(c3*p2+c2*p3)**4)
        alpha_err = np.sqrt(alpha_err_term1+alpha_err_term2+alpha_err_term3) 

    elif calc == 'custom':
        dms.warning_in_orange('custom selectivity not yet implemented')
        #result = custom(c1,c2)
    else:
         dms.error_in_red('Unknown calculation style')

    if custom_err:
        dms.warning_in_orange('custom selectivity errors not yet implemented')
        alpha_err = np.zeros(len(c1))

    return [alpha, alpha_err]