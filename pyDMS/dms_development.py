'''
pyDMS.dms
                          
Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import warnings
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.optimize import minimize
from scipy.special import erfcinv
from scipy.stats import median_abs_deviation, chi2
import statsmodels.api as sm

import pyDMS
from . import report
from . import evaluate

class Gas:
    '''
    *
    Description of each parameter:
    gas.name: formula of gas (e.g., CO2)
    gas.c: arrays of concentrations. Defined as np.array([[gas 1],[gas 2],[etc...]])
    gas.p: arrays of pressures/fugacities. Defined as np.array([[gas 1],[gas 2],[etc...]])
    gas.cerr: arrays of concentration errors. Defined as np.array([[gas 1],[gas 2],[etc...]])
    gas.T: array of temperatures in K. Defined as np.array([gas 1,gas 2,etc...])
    c_f:
    kd_f:
    b_f
    'ch_vec_f'
    S_inf:
    '''


    __slots__ = ['formula', 'c', 'p', 'f', 'cerr', 'T', 'b0',
                 'kd_f', 'b_f', 'ch_vec_f', 
                 'kd_err', 'b_err', 'ch_err',
                 'LFER', 'vH', 'analysis']

    def __init__(self):

        # initializing __slots__
        self.formula = None
        self.c = None
        self.p = None
        self.f = None
        self.cerr = None
        self.T = None

        self.ch_vec_f = None
        self.kd_f = None
        self.b_f = None
        self.ch_err = None
        self.kd_err = None
        self.b_err = None

        self.LFER = LFER()

        self.vH = vH()

        self.analysis = analysis()

    def __repr__(self):
        return 'str'

class LFER:
    '''
    LFER.out gas.out: linear fit results from LFER fitting: [slope_kd, int_kd, slope_b, int_b]
    LFER.SE gas.SE: error from LFER fitting: [slope_kd_err, int_kd_err, slope_b_err, int_b_err]
    LFER.pars gas.pars: LFER data: [log_kd0, deltaHd, log_b0, deltaHb]
    LFER.kd0 gas.kd0: kd0 derived from LFER results
    LFER.b0 gas.b0: b0 derived from LFER results
    '''
    
    __slots__ = ['out', 'SE', 'pars', 'pars_outliers','settings']
    
    def __init__(self):
        self.out = None
        self.SE = None
        self.pars = None
        self.pars_outliers = None
        self.settings = None

    def __repr__(self):
        return 'str'

class vH:
    '''
    dms:
    dms_no_out:
    avg_dms:
    plusminus_1:
    plusminus_2:
    plusminus_3:
    plusminus_4:
    plusminus_ch:
    'c_f':
    res:
    'hessian_average'
    '''
    __slots__ = ['dms', 'dms_no_out', 'avg_dms', 'plusminus_1', 'plusminus_2','plusminus_ch','c_f','res', 'settings', 'hessian_matrix']

    def __init__(self):
        self.dms = None
        self.dms_no_out = None
        self.avg_dms = None
        self.plusminus_1 = None
        self.plusminus_2 = None
        self.plusminus_ch = None
        self.c_f = None # no reason to store this, just calculate it in vis
        self.res = None
        self.hessian_matrix = None
        self.settings = None

    def __repr__(self):
        return 'str'

class analysis:

    __slots__ = ['S_inf', 'S_inf_err', 'deltaH_S_inf', 'deltaH_S_inf_err', 'deltaH_D', 'deltaH_D_err', 'deltaH_b', 'deltaH_b_err']

    def __init__(self):
        self.S_inf = None
        self.S_inf_err = None
        self.deltaH_S_inf = None
        self.deltaH_S_inf_err = None
        self.deltaH_D = None
        self.deltaH_D_err = None
        self.deltaH_b = None
        self.deltaH_b_err = None
        
    def __repr__(self):
        return 'str'

def LFER_loss(x, gas, func='chi2'):
    '''Defines the loss function for the LFER optimization.

    Args:
        x: An array of the variables to be optimized
        gas: An instance of the Gas class

    Returns:
        The loss function result as a number
    '''
    settings = gas.LFER.settings

    # finding settings
    if settings is None:
        settings = {}

    settings.setdefault('ch_constraint', 'decrease')
    ch_constraint = settings.get('ch_constraint')

    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.cerr
    p_vec = gas.p
    T_vec = gas.T

    # initializing error array
    ssr = np.zeros(len(c_vec))

    # initializing variables to solve for
    kd0 = x[0]
    dHD = x[1]
    b0 = x[2]
    dHb = x[3]

    if ch_constraint == 'decrease':
        ch_vec = x[4:]
    elif ch_constraint == 'linear':
        slope_ch = x[4]
        int_ch = x[5]
    else:
        pyDMS.error_in_red('Ch constraint not recognized. Use "decrease" or "linear"')

    # looping through each concentration vector
    for i,_ in enumerate(c_vec):
        
        # extracting parameters for relevant concentrations
        T = T_vec[i]
        p = p_vec[i]
        c = c_vec[i]
        cerr = cerr_vec[i]

        # defining Linear Free Energy Relationships (LFERS)
        kd = kd0 * np.exp(-dHD * 1000 / (8.314 * T))
        b = b0 * np.exp(-dHb * 1000 / (8.314 * T))

        if ch_constraint == 'decrease':
            ch = ch_vec[i]
        elif ch_constraint == 'linear':
            ch = slope_ch * T_vec[i] + int_ch
        else:
            pyDMS.error_in_red('Ch constraint not recognized. Use "decrease" or "linear"')
        
        # solving the DMS model
        calc = kd * p + ch * b * p / (1 + b * p)

        if func == 'chi2':
        # using the chi-squared error as the loss function
            ssr[i] = np.sum(((calc - c) ** 2) / (cerr**2))
            random_thing = None
        elif func == 'least-sq':
            ssr[i] = np.sum((calc - c) ** 2)
            #random_thing = None
        else:
            return pyDMS.error_in_red('Loss function not recognized')
    # Sum of all the errors as metric
    out = np.sum(ssr)

    # * see if scaling is really necessary
    return out/1000

def vH_loss(x, gas, func='chi2'):
    '''Defines the loss function for the van't Hoff optimization.

    Args:
        x: An array of the variables to be optimized
        gas: An instance of the Gas class with data from LFER_loss populated

    Returns:
        The loss function result as a number
    '''
    
    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.cerr
    p_vec = gas.p
    T_vec = gas.T

    # initializing error array
    ssr = np.zeros(len(c_vec))

    # initializing variables to solve for
    dHD = x[0]
    dHb = x[1]
    ch_vec = x[2:]

    # looping through each concentration vector
    for i,_ in enumerate(c_vec):

        # extracting parameters for relevant concentrations
        T = T_vec[i]
        ch = ch_vec[i]
        p = p_vec[i]
        c = c_vec[i]
        cerr = cerr_vec[i]
        LFE_params = gas.LFER.out
        a_kd0 = LFE_params[0]
        b_kd0 = LFE_params[1]
        a_b0 = LFE_params[2]
        b_b0 = LFE_params[3]

        # defining van't Hoff expressions
        kd0 = np.exp((dHD-b_kd0)/a_kd0)
        b0 = np.exp((dHb-b_b0)/a_b0)
        kd = kd0 * np.exp(-dHD * 1000 / (8.314 * T))
        b = b0 * np.exp(-dHb * 1000 / (8.314 * T))

        # solving the DMS model
        calc = kd * p + ch * b * p / (1 + b * p)

        if func == 'chi2':
            # using the chi-squared error as the loss function
            ssr[i] = np.sum(((calc - c) ** 2) / (cerr**2))
            random_thing = None
        elif func == 'least-sq':
            #ssr[i] = np.sum((calc - c) ** 2)
            random_thing = None
        else:
            return pyDMS.error_in_red('Loss function not recognized')
        
    # Sum of all the errors as metric 
    out2 = np.sum(ssr)

    # * see if scaling is really necessary
    return out2

def is_outlier(arr):
    '''Finds otuliers in an array.
    
    Utilizes the definition of isoutlier() from MATLAB
    (see https://www.mathworks.com/help/matlab/ref/isoutlier.html).

    Args:
        arr: An array to find outliers in

    Returns:
        An array with a Boolean mask for whether an entry is an outlier (TRUE) or not (FALSE).
        For example:

            if arr is [0.1,0.12,0.11,100], it will return [False, False, False, True]
    '''

    median = np.median(arr)
    mad = median_abs_deviation(arr)
    c = -1 / (np.sqrt(2) * erfcinv(3/2))  # MATLAB's scaling factor
    scaled_mad = c * mad

    return np.abs(arr - median) > 3 * scaled_mad  # Boolean mask for outliers

def hess(gas, soln):
    '''Solves the analytical Hessian for the vH_loss loss function

    Args:
        gas: An instance of the Gas class
        soln: The object from a scipy.minimize call with loss function vH_loss
    
    Returns:
        A matrix with the Hessian for the result of vH_loss
    '''
    
    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.cerr
    p_vec = gas.p
    T_vec = gas.T

    LFE_params = gas.LFER.out

    a_kd0 = LFE_params[0]
    b_kd0 = LFE_params[1]
    a_b0 = LFE_params[2]
    b_b0 = LFE_params[3]

    dHD = soln.x[0]
    dHb = soln.x[1]
    
    ch_vec = soln.x[2:]

    # initializing arrays for the Hessian
    dCdHD = np.zeros_like(c_vec)
    dCdHb = np.zeros_like(c_vec)
    dCdCH = np.zeros_like(c_vec)
    hessian = np.zeros((len(c_vec)+2,len(c_vec)+2))/1e6
    
    # looping through each concentration vector
    for i,_ in enumerate(c_vec):

        # extracting parameters for relevant concentrations
        T = T_vec[i]
        ch = ch_vec[i]
        p = p_vec[i]
        c = c_vec[i]
        cerr = cerr_vec[i]

        # solving derivatives for the Hessian
        dCdHD[i] = np.exp(-120.279*dHD/T+(-b_kd0+dHD)/a_kd0)*p*(-120.279/T+1/a_kd0)
        dCdHb[i] = (ch*np.exp(1+(120.279*dHb/T)+(b_b0+dHb)/a_b0)*p*(T-120.279*a_b0))/((np.exp(b_b0/a_b0+120.279*dHb/T)+np.exp(dHb/a_b0)*p)**2*T*a_b0)
        dCdCH[i] = (1+(np.exp((b_b0-dHb)/a_b0+120.279*dHb/T))/p)**-1

    # inserting values into Hessian
    hessian[0,0] = np.sum(2*1/cerr**2*dCdHD*dCdHD)
    hessian[0,1] = np.sum(2*1/cerr**2*dCdHD*dCdHb)
    hessian[1,0] = np.sum(2**1/cerr**2*dCdHD*dCdHb)
    hessian[1,1] = np.sum(2**1/cerr**2*dCdHb*dCdHb)

    for i, dCdCH_vec in enumerate(dCdCH):
        hessian[0,i+2] = np.sum(2**1/cerr**2*dCdHD*dCdCH_vec)
        hessian[i+2,0] = np.sum(2**1/cerr**2*dCdHD*dCdCH_vec)
        hessian[1,i+2] = np.sum(2**1/cerr**2*dCdHb*dCdCH_vec)
        hessian[i+2,1] = np.sum(2**1/cerr**2*dCdHb*dCdCH_vec)
        hessian[i+2,i+2] = np.sum(2**1/cerr**2*dCdCH_vec*dCdCH_vec)

    return hessian

def calc_LFEs(gas, settings=None):
    '''Implements the LFER_loss function to optimize the DMS model based on LFER constraints

    Args:
        gas: An instance of the Gas class
        trials: The number of optimization cycles to run with randomly chosen initial guesses
        bounds: The bounds for each guess to be chosen between
        solver: The scipy.minimize solver (SLSQP or trust-constr) to use
        verbose: Whether information should be printed to the display
        solve_verbose: Whether individual solver iterations should be printed to the display

    Returns:
        Data in the Gas.LFER subclass
    '''

    # retriving gas data
    c = gas.c
    cerr = gas.cerr
    p = gas.p
    T = gas.T

    settings = gas.LFER.settings

    # finding settings
    if settings is None:
        settings = {}

    settings.setdefault('dHD_bounds', [-1, -30])
    settings.setdefault('dHb_bounds', [-1, -30])
    settings.setdefault('kD0_bounds', [0.001, 0.01])
    settings.setdefault('b0_bounds', [-0.0001, 0.005])
    settings.setdefault('ch_bounds', np.array([[0, 150] for _ in range(len(c))]))
    settings.setdefault('ch_constraint', 'decrease')
    settings.setdefault('trials', 15)
    settings.setdefault('solver', 'SLSQP')
    settings.setdefault('verbose', True)
    settings.setdefault('solver_verbose', False)

    kd0_0_bnd = settings.get('kD0_bounds')
    dHD0_0_bnd = settings.get('dHD_bounds')
    b0_0_bnd = settings.get('b0_bounds')
    dHb_0_bnd = settings.get('dHb_bounds')
    ch_0_bnd = settings.get('ch_bounds')
    ch_constraint = settings.get('ch_constraint')

    trials = settings.get("trials")
    solver = settings.get("solver")
    verbose = settings.get("verbose")
    solver_verbose = settings.get("solver_verbose")

    # * consider allowing this to be turned off?

    if ch_constraint == 'decrease':
        # setting up the solver constraints: A*x <= b
        A_con = np.zeros((len(c)-1,len(c)+4))
        b_con = np.zeros(len(c)-1)

        # setting up the linear constraint so that each ch' must be lower than the previous
        for i in range(len(c)-1):
            first_index = i+4
            second_index = i+5
            A_con[i,first_index] = -1
            A_con[i,second_index] = 1

        def linear_constraint(x):
            return b_con - (A_con @ x)

    # checks to ensure proper data formatting
    if np.shape(c) != np.shape(cerr) != np.shape(p):
        pyDMS.error_in_red('c, cerr, and p dimensions do not match') # *

    elif len(c) != len(T):
        pyDMS.error_in_red('T does not match c')


    # setting up arrays for storing minimization results
    if ch_constraint == 'decrease':
        nOptVars = 4 + len(c) # no. of optimization variables\
    elif ch_constraint == 'linear':
        nOptVars = 6
    else:
        pyDMS.error_in_red('Ch constraint not recognized. Use "decrease" or "linear"')

    dms = np.zeros((trials,nOptVars)) # results *
    res = np.zeros(trials) # func(x)
    flag = np.zeros(trials) # optimization flags
    avg_dms = np.zeros(nOptVars) # averaged results
    std_dev = np.zeros(nOptVars) # standard deviations
    num_par_final = np.zeros(nOptVars) # no. of non-outlier chains


    # finding the range of each bound
    kdval = kd0_0_bnd[1] - kd0_0_bnd[0]
    dHDval = dHD0_0_bnd[1] - dHD0_0_bnd[0]
    b0val = b0_0_bnd[1] - b0_0_bnd[0]
    dHbval = dHb_0_bnd[1] - dHb_0_bnd[0]
    ch0_vals = [bnd[1] - bnd[0] for bnd in ch_0_bnd]

    ch_0 = np.zeros(len(ch0_vals)) # holding ch_0 guesses

    # printing bounds to search through
    if verbose:
        print('----------------------------LFER Initial Bounds----------------------------')
        print(f'kd0_0: {kd0_0_bnd}')
        print(f'dHd0_0: {dHD0_0_bnd}')
        print(f'b0_0: {b0_0_bnd}')
        print(f'dHb0_0: {dHb_0_bnd}')
        for i, ch_0_bnds in enumerate(ch_0_bnd):
            print(f'C_H\'{i}: {ch_0_bnds}')
        print('---------------------------------------------------------------------------')

    # initializing random number generator
    rng = np.random.default_rng()

    for j in range(trials):

        if verbose and np.mod(j,10)==0:
            print(f'LFER trial: {j}/{trials}')
                  
        # picking a random number in each bound
        kd0_0 = kd0_0_bnd[0] + kdval*rng.random() 
        dHD_0 = dHD0_0_bnd[0] + dHDval*rng.random()
        b0_0 = b0_0_bnd[0] + b0val*rng.random()
        dHb_0 = dHb_0_bnd[0] + dHbval*rng.random()

        for i, bounds in enumerate(ch_0_bnd):
            ch_0[i]=bounds[0]+ch0_vals[i]*rng.random()

        # setting up initial guess for solver
        if ch_constraint == 'decrease':
            x0 = np.concatenate([[kd0_0, dHD_0, b0_0,dHb_0],ch_0])
        elif ch_constraint == 'linear':
            slope_ch_0 = rng.uniform(-0.9, -0.1)
            int_ch_0 = rng.uniform(100, 200)
            x0 = np.array([kd0_0, dHD_0, b0_0, dHb_0, slope_ch_0, int_ch_0])
        else:
            pyDMS.error_in_red('*oops')
        # * consider storing the inital guess inside initial_g

        # * add custom bounds
        # Setting bounds for solver. Generally less rigorous than initial guesses
        #                kd0_0     dHD_0     b0_0       dHb_0
        solver_bounds = [(0,None), (-50, 0), (0, None), (-50, 0)]

        if ch_constraint == 'decrease':
            solver_bounds.extend([(ch_bnd[0], ch_bnd[1]) for ch_bnd in ch_0_bnd])
        elif ch_constraint == 'linear':
            solver_bounds.extend([(-0.9, -0.1), (100, 200)])
        else:
            pyDMS.error_in_red('oops*')

        if ch_constraint == 'decrease':
            solver_constraint = {'type': 'ineq', 'fun': linear_constraint}
        elif ch_constraint == 'linear':
            solver_constraint = None
        else:
            pyDMS.error_in_red('*oops')

        # optimizing with the 'trust-constr' algorithm
        if solver=='trust-constr':

            settings.setdefault('xtol', 1e-12)
            settings.setdefault('gtol', 1e-12)
            settings.setdefault('maxiter', 1000)
            
            xtol = settings.get('xtol')
            gtol = settings.get('gtol')
            maxiter = settings.get('maxiter')

            options = {
                    'maxiter': maxiter,        # maximum iterations
                    'xtol': xtol,              # variable tolerance
                    'gtol': gtol,              # function tolerance
                    'disp': solver_verbose,    # display optimization progress
            }

            # finding initial guess
            LFER_loss(x0, gas)

            if not verbose:
                warnings.filterwarnings("ignore", message="delta_grad == 0.0")

            # running optimization
            result = minimize(
                LFER_loss,
                x0,
                args=(gas,),
                method='trust-constr',
                bounds=solver_bounds,
                constraints={'type': 'ineq', 'fun': linear_constraint},
                options=options)

        # optimizing with the 'SLQSP' algorithm
        elif solver=='SLSQP':

            settings.setdefault('ftol', 1e-12)
            settings.setdefault('maxiter', 1000)
            
            ftol = settings.get('ftol')
            maxiter = settings.get('maxiter')
        
            options = {
                    'maxiter': maxiter,     # max iterations
                    'ftol': ftol,           # function tolerance 
                    'disp': solver_verbose, # display iteration results
            }

            # finding initial guess
            LFER_loss(x0,gas)

            # running optimization
            result = minimize(
                LFER_loss,
                x0,
                args=(gas),
                method='SLSQP',
                bounds=solver_bounds,
                constraints=solver_constraint,
                options=options)

        # collecting optimization results
        flag[j] = result.status
        res[j] = result.fun
        dms[j,:] = result.x

        transposed_dms = np.transpose(dms)


    # FINDING OUTLIERS
    outlier_track = np.zeros(nOptVars)  # store number of outliers for each parameter
    par_outliers_removed = np.copy(transposed_dms)  # Copy of the transposed results matrix
    par_with_outliers = np.copy(transposed_dms) # Copy of the transposed results matrix


    # finding outliers for each parameter
    for i in range(nOptVars):

        par = transposed_dms[i] # parameter of interest
        par_no_outliers = np.copy(par) # copying parameter of interest to separate array

        TF = is_outlier(par) # finding if iteration parameter is an outlier (TRUE if outlier)

        par_no_outliers[TF] = 0  # setting values to zero if outlier * unclear why this is needed
        par_outliers_removed[i][TF] = 0  # setting outliers to zero in complete array as well
        par = par[~TF] # removing values that are outliers from array

        outliers = np.where(TF == 1)[0]  # indices where outliers are present
        # * removed an if statement and should still work: keep an eye an it
        outlier_track[i] = len(outliers) # store number of outliers for each parameter

        avg_dms[i] = np.mean(par) # finding average of each cleaned parameter
        std_dev[i] = np.std(par) # finding standard deviation of each cleaned parameter
        num_par_final[i] = len(par) # number of remaining parameter iterations post-cleaning

    # *TO DO: NEED TO EXAMINE THIS SOME MORE
    par_outliers_removed = np.transpose(par_outliers_removed)
    par_outliers_removed = par_outliers_removed[~np.any(par_outliers_removed == 0, axis=1)]
    par_outliers_removed = np.transpose(par_outliers_removed)

    # For plotting purposes
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")  # This ensures that the warning is captured
        log_kd0_out = np.log(par_with_outliers[0])
        deltaHd_out = par_with_outliers[1]
        log_b0_out = np.log(par_with_outliers[2])
        deltaHb_out = par_with_outliers[3]


    # FINDING LFER PARAMETERS
    # deltaHD = aD*ln(kd0) + bD
    log_kd0 = np.log(par_outliers_removed[0])
    deltaHd = par_outliers_removed[1]
    log_b0 = np.log(par_outliers_removed[2])
    deltaHb = par_outliers_removed[3]

    # ensuring an intercept is fitted as well
    log_kd0_with_const = sm.add_constant(log_kd0)
    log_b0_with_const = sm.add_constant(log_b0)

    # fitting linear regressions
    deltaHd_model = sm.OLS(deltaHd, log_kd0_with_const).fit()
    deltaHb_model = sm.OLS(deltaHb, log_b0_with_const).fit()

    int_kd, slope_kd = deltaHd_model.params # int_kd = aD, slope_kD = bD
    int_kd_err, slope_kd_err = deltaHd_model.bse # errors in int_kD and slope_kD
    int_b, slope_b = deltaHb_model.params # int_b= ab, slope_b = bb
    int_b_err, slope_b_err = deltaHb_model.bse # errors in int_b and slope_b

    out = [slope_kd, int_kd, slope_b, int_b] # collecting slopes, intercepts
    SE = [slope_kd_err, int_kd_err, slope_b_err, int_b_err] # collecting errors

    # collecting data
    # * keep an eye on this. Changed from len(c) to 4
    pars = np.zeros((len(log_kd0), 4))
    pars[:,0] = log_kd0
    pars[:,1] = deltaHd
    pars[:,2] = log_b0
    pars[:,3] = deltaHb

    # collecting data with outliers for plotting
    # * keep an eye on this. Changed from len(c) to 4
    pars_outliers = np.zeros((len(log_kd0_out), len(c)))
    pars_outliers[:,0] = log_kd0_out
    pars_outliers[:,1] = deltaHd_out
    pars_outliers[:,2] = log_b0_out
    pars_outliers[:,3] = deltaHb_out

    gas.LFER.out = out
    gas.LFER.SE = SE
    gas.LFER.pars = pars
    gas.LFER.pars_outliers = pars_outliers

    gas.LFER.settings = settings
    
    return gas

def calc_params(gas, settings=None):
    '''Implements the vH_loss function to optimize the DMS model based on van't Hoff constraints

    Args:
        gas: An instance of the Gas class
        trials: The number of optimization cycles to run with randomly chosen initial guesses
        bounds: The bounds for each guess to be chosen between
        solver: The scipy.minimize solver (SLSQP or trust-constr) to use
        verbose: Whether information should be printed to the display
        solve_verbose: Whether individual solver iterations should be printed to the display

    Returns:
        Data in the Gas.vH
          subclass
    '''

    # pulling in data
    c = gas.c
    cerr = gas.cerr
    p = gas.p
    T = gas.T
    LFE_params = gas.LFER.out

    settings = gas.vH.settings

    # finding settings
    if settings is None:
        settings = {}
    
    settings.setdefault('dHD_bounds', [-1, -30])
    settings.setdefault('dHb_bounds', [-1, -30])
    settings.setdefault('ch_bounds', np.array([[0, 150] for _ in range(len(c))]))
    settings.setdefault('trials', 15)
    settings.setdefault('solver', 'SLSQP')
    settings.setdefault('verbose', True)
    settings.setdefault('solver_verbose', False)

    dHD0_0_bnd = settings.get('dHD_bounds')
    dHb_0_bnd = settings.get('dHb_bounds')
    ch_0_bnd = settings.get('ch_bounds')

    trials = settings.get("trials")
    solver = settings.get("solver")
    verbose = settings.get("verbose")
    solver_verbose = settings.get("solver_verbose")

    # setting up the solver constraints: A*x <= b
    A_con = np.zeros((len(c)-1,len(c)+2))
    b_con = np.zeros(len(c)-1)

    if np.shape(c) != np.shape(cerr) != np.shape(p):
        pyDMS.error_in_red('c, cerr, and p dimensions do not match')

    elif len(c) != len(T):
        pyDMS.error_in_red('*some error')

    nOptVars = 2 + len(c) # number of optimization variables
    dms = np.zeros((trials,nOptVars)) # creating a matrix for holding results
    res = np.zeros(trials)
    flag = np.zeros(trials)
    avg_dms = np.zeros(nOptVars)
    std_dev = np.zeros(nOptVars)
    num_par_final = np.zeros(nOptVars)
    avg_dms_raw = np.zeros(nOptVars)
    std_dev_raw = np.zeros(nOptVars)
    hessian_matrix = np.zeros((len(c)+2, len(c)+2, trials))

    # finding the range of each bound
    dHDval = dHD0_0_bnd[1] - dHD0_0_bnd[0]
    dHbval = dHb_0_bnd[1] - dHb_0_bnd[0]
    ch0_vals = [bnd[1] - bnd[0] for bnd in ch_0_bnd]

    ch_0 = np.zeros(len(ch0_vals)) # # holding ch_0 guesses

    # printing bounds to search through
    if verbose:
        print('----------------------Van\'t Hoff Initial Bounds----------------------')
        print(f'dHD0_0: {dHD0_0_bnd}')
        print(f'dHb0_0: {dHb_0_bnd}')
        for i, ch_0_bnds in enumerate(ch_0_bnd):
            print(f'C_H\'{i}: {ch_0_bnds}')
        print('----------------------------------------------------------------------')

    # initializing random number generator
    rng = np.random.default_rng()

    for j in range(trials):

        if verbose and np.mod(j,10)==0:
            print(f'van\'t Hoff trial: {j}/{trials}')

        # picking a random number in each bound
        dHD_0 = dHD0_0_bnd[0] + dHDval*rng.random()
        dHb_0 = dHb_0_bnd[0] + dHbval*rng.random()

        for i, bounds in enumerate(ch_0_bnd):
            ch_0[i]=bounds[0]+ch0_vals[i]*rng.random()

        # setting up initial guess for solver
        x1 = np.concatenate([[dHD_0, dHb_0],ch_0])

        def linear_constraint(x):
            return b_con - (A_con @ x)

        # * consider storing the inital guess inside initial_g
        for i in range(len(c)-1):
            first_index = i+2
            second_index = i+3
            A_con[i,first_index] = -1
            A_con[i,second_index] = 1

        # * add custom bounds
        # Setting bounds for solver. Generally less rigorous than initial guesses
        #                kd0_0     dHD_0     b0_0       dHb_0
        solver_bounds = [(-50, 0), (-50, 0)]

        solver_bounds.extend([(ch_bnd[0], ch_bnd[1]) for ch_bnd in ch_0_bnd])

        # optimizing with the 'trust-constr' algorithm
        if solver=='trust-constr':
            
            settings.setdefault('xtol', 1e-12)
            settings.setdefault('gtol', 1e-12)
            settings.setdefault('maxiter', 1000)
            
            xtol = settings.get('xtol')
            gtol = settings.get('gtol')
            maxiter = settings.get('maxiter')

            options = {
                    'maxiter': maxiter,     # maximum iterations
                    'xtol': xtol,           # variable tolerance
                    'gtol': gtol,           # function tolerance
                    'disp': solver_verbose, # display optimization progress
            }

            # finding initial guess
            vH_loss(x1,gas)

            #flag = np.zeros((1,j))
            if not verbose:
                warnings.filterwarnings("ignore", message="delta_grad == 0.0")

            # running optimization
            result = minimize(
                vH_loss,
                x1,
                args=(gas),
                method='trust-constr',
                bounds=solver_bounds,
                constraints={'type': 'ineq', 'fun': linear_constraint},
                options=options #,
                #hess=hessian_approx
                )

        # optimizing with the 'SLQSP' algorithm
        elif solver=='SLSQP':
            
            settings.setdefault('ftol', 1e-12)
            settings.setdefault('maxiter', 1000)
            
            ftol = settings.get('ftol')
            maxiter = settings.get('maxiter')
        
            options = {
                    'maxiter': maxiter,     # max iterations
                    'ftol': ftol,           # function tolerance 
                    'disp': solver_verbose, # display iteration results
            }

            # finding initial guess
            vH_loss(x1,gas)

            # * create flag array
            #flag = np.zeros((1,j))

            # running optimization
            result = minimize(
                vH_loss,
                x1,
                args=(gas,),
                method='SLSQP',
                bounds=solver_bounds,
                constraints={'type': 'ineq', 'fun': linear_constraint},
                options=options)
            
        hessian_matrix[:,:,j] = hess(gas, result)
        flag[j] = result.status
        res[j] = result.fun
        dms[j,:] = result.x

        transposed_dms = np.transpose(dms)
        par2_outliers_removed = np.copy(transposed_dms)  # Copy of the transposed results matrix
    
    outlier_track = np.zeros(nOptVars)  # store number of outliers for each parameter

    for i in range(len(x1)):

        par2 = transposed_dms[i]
        avg_dms_raw = np.mean(par2)
        std_dev_raw = np.std(par2)

        TF = is_outlier(par2)
        par2 = par2[~TF]

        par2_outliers_removed[i][TF] = 0  # setting outliers to zero in complete array as well

        outliers = np.where(TF == 1)[0]  # indices where outliers are True

        outlier_track[i] = len(outliers) # store number of outliers for each parameter

        avg_dms[i] = np.mean(par2) # finding average of each cleaned parameter
        std_dev[i] = np.std(par2) # finding standard deviation of each cleaned parameter
        num_par_final[i] = len(par2) # number of remaining parameter iterations post-cleaning

    par2_outliers_removed = np.transpose(par2_outliers_removed)
    par2_outliers_removed = par2_outliers_removed[~np.any(par2_outliers_removed == 0, axis=1)]

    dHD_f = avg_dms[0]
    dHb_f = avg_dms[1]
    ch_vec_f = avg_dms[2:]
    a_kd0_f = LFE_params[0]
    b_kd0_f = LFE_params[1]
    a_b0_f = LFE_params[2]
    b_b0_f = LFE_params[3]
    kd0_f = np.exp((dHD_f-b_kd0_f)/a_kd0_f)
    b0_f = np.exp((dHb_f-b_b0_f)/a_b0_f)

    kd_f = np.zeros(len(T))
    b_f = np.zeros(len(T))

    C_f = np.zeros((len(T),100))

    for i, temp in enumerate(T):
        max_press = np.max(p[i])
        press = np.linspace(1E-6,max_press,100)
        kd_f[i] = kd0_f*np.exp(-dHD_f*1000/(8.314*temp))
        b_f[i] = b0_f*np.exp(-dHb_f*1000/(8.314*temp))
        C_f[i,:] = kd_f[i]*press+ch_vec_f[i]*b_f[i]*press/(1+b_f[i]*press) 

    gas.vH.c_f = C_f
    gas.kd_f = kd_f
    gas.b_f = b_f
    gas.ch_vec_f = ch_vec_f
    gas.vH.dms = dms
    gas.vH.dms_no_out = par2_outliers_removed
    gas.vH.avg_dms = avg_dms
    gas.vH.res = res
    gas.vH.hessian_matrix = hessian_matrix

    gas.vH.settings = settings

    return gas

def chi2_error_fit(gas):
    '''
    *
    '''
    
    hessian_matrix = gas.vH.hessian_matrix
    res = gas.vH.res
    avg_dms = gas.vH.avg_dms
    ch_vec_f = gas.ch_vec_f
    
    # mean of the Hessians across all trials
    hessian_avg = np.mean(hessian_matrix, axis=2)

    # mean of the residuals across all trials
    f = np.mean(res)

    # mean of the residuals shifted by the critical chi-square value
    chi_crit = f + chi2.isf(0.05, len(avg_dms))

    # setting up arrays for error analysis
    avg_dms1 = np.copy(avg_dms)
    avg_dms2 = np.copy(avg_dms)
    avg_dms_ch = [np.copy(avg_dms) for _ in range(len(ch_vec_f))] # making sure arrays are not linked

    # finding dHD and dHb +/- 1 away from optimal value
    dHD_vec = np.linspace(avg_dms[0] - 1, avg_dms[0] + 1, 100)
    dHb_vec = np.linspace(avg_dms[1] - 1, avg_dms[1] + 1, 100)

    # finding C_H' +/- 1.5 away from optimal value
    ch_vec = [np.linspace(avg_dms_val - 1.5, avg_dms_val + 1.5, 100) for avg_dms_val in avg_dms[2:]]

    # setting up more arrays for error analysis
    chi2_1 = np.zeros(len(dHD_vec))
    chi2_2 = np.zeros(len(dHD_vec))
    chi2_ch = [np.zeros(len(dHD_vec)) for _ in range(len(ch_vec_f))]

    # perturbijng the parameters to create a chi-square distribution
    for i, dHD_val in enumerate(dHD_vec):
        avg_dms1[0] = dHD_val
        avg_dms2[1] = dHb_vec[i]

        for j in range(len(ch_vec_f)):        
            avg_dms_ch[j][j+2] = ch_vec[j][i]

        # solving for the chi-square distributions for each parameter as a Taylor series expansion
        chi2_1[i] = f + (1/2)*((avg_dms1-avg_dms) @ hessian_avg @ np.transpose((avg_dms1-avg_dms)))
        chi2_2[i] = f + (1/2)*((avg_dms2-avg_dms) @ hessian_avg @ np.transpose((avg_dms2-avg_dms)))

        for j in range(len(ch_vec_f)):
            chi2_ch[j][i] = f + (1/2)*((avg_dms_ch[j]-avg_dms) @ hessian_avg @ np.transpose((avg_dms_ch[j]-avg_dms)))

    # fitting the chi-square distributions to parabolas
    chi2_1_fit = poly.polyfit(dHD_vec, chi2_1,2)
    chi2_2_fit = poly.polyfit(dHb_vec, chi2_2,2)
    chi2_ch_fit = [np.zeros(3) for _ in range(len(ch_vec_f))]

    for j in range(len(ch_vec_f)):
        chi2_ch_fit[j] = poly.polyfit(ch_vec[j], chi2_ch[j],2)

    # shifting the fits by the chi-square minimum
    chi2_1_fit[0] = chi2_1_fit[0] - chi_crit
    chi2_2_fit[0] = chi2_2_fit[0] - chi_crit

    for j in range(len(ch_vec_f)):
        chi2_ch_fit[j][0] = chi2_ch_fit[j][0] - chi_crit

    # finding where the distributions are equal to the critical chi-square value
    roots_1 = poly.polyroots(chi2_1_fit)
    roots_2 = poly.polyroots(chi2_2_fit)

    roots_ch = [np.zeros(2) for _ in range(len(ch_vec_f))]

    for j in range(len(ch_vec_f)):
        roots_ch[j] = poly.polyroots(chi2_ch_fit[j])

    plusminus_1 = roots_1 - avg_dms[0] # error for dHD
    plusminus_2 = roots_2 - avg_dms[1] # error for dHb

    plusminus_ch = [np.zeros(2) for _ in range(len(ch_vec_f))]

    for j in range(len(ch_vec_f)):
        plusminus_ch[j] = roots_ch[j] - avg_dms[j+2]

    ch_err = np.array([arr[1] for arr in plusminus_ch]) # error for C_H'

    gas.vH.plusminus_1 = plusminus_1
    gas.vH.plusminus_2 = plusminus_2
    gas.vH.plusminus_ch = plusminus_ch
    gas.ch_err = ch_err

    return gas

def propogate_error(gas):
    '''Determines error from the van't Hoff optimization and propogates the error to further parameters

    Args:
        gas: An instance of the Gas class with Gas.vH data populated

    Returns:
        Data in the Gas.vH subclass and Gas class
    '''
    print('-----------------------------------------------------------------')
    print('Optimization successful; starting error propogation')

    LFE_params = gas.LFER.out
    LFE_error = gas.LFER.SE
    avg_dms = gas.vH.avg_dms
    T_vec = gas.T
    plusminus_1 = gas.vH.plusminus_1
    plusminus_2 = gas.vH.plusminus_2

    a_kd0 = LFE_params[0]
    b_kd0 = LFE_params[1]
    a_b0 = LFE_params[2]
    b_b0 = LFE_params[3]

    err_a_kd0 = LFE_error[0]
    err_b_kd0 = LFE_error[1]
    err_a_b0 = LFE_error[2]
    err_b_b0 = LFE_error[3]

    dHD_f = avg_dms[0]
    dHb_f = avg_dms[1]

    delkddel_a_kd0 = np.zeros(len(T_vec))
    delkddel_b_kd0 = np.zeros(len(T_vec))
    delkddelT = np.zeros(len(T_vec))
    delkddelHd = np.zeros(len(T_vec))
    kd_err = np.zeros(len(T_vec))


    for i, T in enumerate(T_vec):
        delkddelT[i] = (1000*dHD_f/(8.314*T**2))*np.exp(-1000*dHD_f/(8.314*T) + (dHD_f - b_kd0)/a_kd0)
        delkddel_a_kd0[i] = ((b_kd0 - dHD_f)/a_kd0**2)*np.exp(-1000*dHD_f/(8.314*T) + (dHD_f - b_kd0)/a_kd0)
        delkddel_b_kd0[i] = (-1/a_kd0)*np.exp(-1000*dHD_f/(8.314*T) + (dHD_f - b_kd0)/a_kd0)
        delkddelHd[i] = (1/a_kd0 - 1000/8.314/T)*np.exp(-1000*dHD_f/(8.314*T) + (dHD_f - b_kd0)/a_kd0)
        kd_err[i] = np.sqrt(0.01**2*delkddelT[i]**2 + err_a_kd0**2*delkddel_a_kd0[i]**2 + err_b_kd0**2*delkddel_b_kd0[i]**2 + plusminus_1[1]**2*delkddelHd[i]**2)

    b_err = np.zeros(len(T_vec))
    delbdel_a_b0 = np.zeros(len(T_vec))
    delbdel_b_b0 = np.zeros(len(T_vec))
    delbdelT = np.zeros(len(T_vec))
    delbdelHb = np.zeros(len(T_vec))

    for i, T in enumerate(T_vec):
        delbdelT[i] = (1000*dHb_f/(8.314*T**2))*np.exp(-1000*dHb_f/(8.314*T) + (dHb_f - b_b0)/a_b0)
        delbdel_a_b0[i] = ((b_b0 - dHb_f)/a_b0**2)*np.exp(-1000*dHb_f/(8.314*T) + (dHb_f - b_b0)/a_b0)
        delbdel_b_b0[i] = (-1/a_b0)*np.exp(-1000*dHb_f/(8.314*T) + (dHb_f - b_b0)/a_b0)
        delbdelHb[i] = (1/a_b0 - 1000/8.314/T)*np.exp(-1000*dHb_f/(8.314*T) + (dHb_f - b_b0)/a_b0)
        b_err[i] = np.sqrt(0.01**2*delbdelT[i]**2 + err_a_b0**2*delbdel_a_b0[i]**2 + err_b_b0**2*delbdel_b_b0[i]**2 + plusminus_2[1]**2*delbdelHb[i]**2)

    gas.kd_err = kd_err
    gas.b_err = b_err
    #gas.ch_err = ch_err

    print('Error propogation successful')
    print('------------------------------------------------------------------')

    return gas

def compute(gas, info = True):
    
    '''Provides a wrapper to run the entire optimization procedure with a single function call

    Args:
        gas: An instance of the Gas class
        trials: The number of optimization cycles to run with randomly chosen initial guesses
        bounds: The bounds for each guess to be chosen between
        solver: The scipy.minimize solver (SLSQP or trust-constr) to use
        verbose: Whether information should be printed to the display
        solve_verbose: Whether individual solver iterations should be printed to the display

    Returns:
        All data in the Gas.LFER and Gas.vH subclasses
    '''

    calc_LFEs(gas)
    calc_params(gas)

    chi2_error_fit(gas)
    propogate_error(gas)

    evaluate.heat_of_sorption(gas)
    #vis.heat_of_sorption(gas)

    if info:
        report.LFER(gas)
        report.LFER(gas, outliers=True)
        report.histograms(gas)
        report.isotherms(gas)
        report.heat_of_sorption(gas)
        report.generate(gas)
