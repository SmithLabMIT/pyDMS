"""
pyDMS.dms
Fitting the dual-mode sorption (DMS) model with LFER and van't Hoff constraints

Copyright 2025 Massachusetts Institute of Technology
Licensed under the MIT License
"""

import warnings
import pickle
import numpy as np
import numpy.polynomial.polynomial as poly
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.special import erfcinv
from scipy.stats import median_abs_deviation, chi2
from scipy.linalg import cho_factor, cho_solve
from sklearn.linear_model import RANSACRegressor, LinearRegression

import pyDMS
from . import report
from . import evaluate
from . import fugacity


class Gas:
    """Holds all input and output data of a pyDMS optimization run

    Attributes:
        formula: a string indicating the chemical formula of the gas
            (e.g., 'CO2').
        temp: An array of temperatures in K. Defined as
            np.array([Temp 1, Temp 2, ...]).
        p: An array of pressures. Defined as
            np.array([[Array 1], [Array 2], ...]).
        f: An array of fugacities. Defined as
        c: An array of concentrations. Defined as
            np.array([[Array 1],[Array 2], ...]).
            [[Array 1], [Array 2], ...].
        c_err: An array of uncertainties in concentration. Defined as
            np.array([[Array 1], [Array 2], ...]).
        Z: An array of compressibility factors. Defined as
            [Array 1], [Array 2], ...].
        kD:  An array of DMS parameter kD. Defined as
            np.array([Param 1, Param 2, ...]).
        b: An array of DMS parameter b. Defined as
            np.array([Param 1, Param 2, ...]).
        CH: An array of DMS parameter CH'. Defined as
            np.array([Param 1, Param 2, ...]).
        kD_err: An array of uncertainties in kD. Defined as
            np.array([Param 1, Param 2, ...]).
        b_err: An array of uncertainties in b. Defined as
            np.array([Param 1, Param 2, ...]).
        CH_err: An array of uncertainties in CH'. Defined as
            np.array([Param 1, Param 2, ...]).
        LFER: An instance of the LFER class containing results from
            the LFER fitting.
        vH: An instance of the vH class containing results from
            the van't Hoff fitting.
        analysis: An instance of the analysis class containing results
            from post-optimization analysis.
        settings: A dictionary of settings used for the optimization runs.
        virial_coeff: A dictionary containing user-supplied Virial
            coefficients.
        pr_coeff: A dictionary containing user-supplied Peng-Robinson
            coefficients.
    """

    __slots__ = [
        "formula",
        "c",
        "p",
        "f",
        "Z",
        "c_err",
        "temp",
        "kD",
        "b",
        "CH",
        "kD_err",
        "b_err",
        "CH_err",
        "LFER",
        "vH",
        "analysis",
        "settings",
        "virial_coeff",
        "pr_coeff",
    ]

    def __init__(self):

        self.formula = None
        self.c = None
        self.p = None
        self.f = None
        self.Z = None
        self.c_err = None
        self.temp = None
        self.CH = None
        self.kD = None
        self.b = None
        self.CH_err = None
        self.kD_err = None
        self.b_err = None
        self.settings = None
        self.virial_coeff = None
        self.pr_coeff = None

        self.LFER = LFER()

        self.vH = vH()

        self.analysis = analysis()

        def __str__(self):
            lines = [f"{self.__class__.__name__}("]
            for attr in self.__slots__:
                value = getattr(self, attr)
                if isinstance(value, np.ndarray):
                    array_str = np.array2string(
                        value,
                        precision=3,
                        suppress_small=True,
                        threshold=10,
                        separator=", ",
                        max_line_width=80,
                    )
                    lines.append(f"  {attr} = {array_str}")
                else:
                    lines.append(f"  {attr} = {str(value)}")
            lines.append(")")
            return "\n".join(lines)


class LFER:
    """Holds output from the LFER fitting optimization.

    Attributes:
        fit: An array of the linear fit results
            [slope_kd, int_kd, slope_b, int_b]
        fit_err: An array of the uncertainty
            [slope_kd_err, int_kd_err, slope_b_err, int_b_err]
        out: LFER inlier output
            [log_kd0, deltaHd, log_b0, deltaHb]
        out_outliers: LFER inlier and ouitlier output
            [log_kd0, deltaHd, log_b0, deltaHb]
    """

    __slots__ = ["fit", "fit_err", "out", "out_outliers"]

    def __init__(self):
        self.fit = None
        self.fit_err = None
        self.out = None
        self.out_outliers = None

    def __str__(self):
        lines = [f"{self.__class__.__name__}("]
        for attr in self.__slots__:
            value = getattr(self, attr)
            if isinstance(value, np.ndarray):
                array_str = np.array2string(
                    value,
                    precision=3,
                    suppress_small=True,
                    threshold=10,
                    separator=", ",
                    max_line_width=80,
                )
                lines.append(f"  {attr} = {array_str}")
            else:
                lines.append(f"  {attr} = {str(value)}")
        lines.append(")")
        return "\n".join(lines)


class vH:
    """Holds output from the LFER fitting optimization.
    # * need to check this
    Attributes:
        fit: An array of the linear fit results
            [slope_kd, int_kd, slope_b, int_b]
        fit_err: An array of the uncertainty
            [slope_kd_err, int_kd_err, slope_b_err, int_b_err]
        out: LFER inlier output
            [log_kd0, deltaHd, log_b0, deltaHb]
        out_outliers: LFER inlier and outlier output
            [log_kd0, deltaHd, log_b0, deltaHb]
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
    """

    __slots__ = [
        "out_outliers",
        "out",
        "avg_dms",
        "plusminus_1",
        "plusminus_2",
        "plusminus_ch",
        "residuals",
        "hessian_matrix",
    ]

    def __init__(self):
        self.out_outliers = None
        self.out = None
        self.avg_dms = None
        self.plusminus_1 = None
        self.plusminus_2 = None
        self.plusminus_ch = None
        self.residuals = None
        self.hessian_matrix = None

    def __str__(self):
        lines = [f"{self.__class__.__name__}("]
        for attr in self.__slots__:
            value = getattr(self, attr)
            if isinstance(value, np.ndarray):
                array_str = np.array2string(
                    value,
                    precision=3,
                    suppress_small=True,
                    threshold=10,
                    separator=", ",
                    max_line_width=80,
                )
                lines.append(f"  {attr} = {array_str}")
            else:
                lines.append(f"  {attr} = {str(value)}")
        lines.append(")")
        return "\n".join(lines)


class analysis:
    """TODO"""

    __slots__ = [
        "S_inf",
        "S_inf_err",
        "deltaH_S_inf",
        "deltaH_S_inf_err",
        "deltaH_D",
        "deltaH_D_err",
        "deltaH_b",
        "deltaH_b_err",
        "c_iso",
        "deltaH_iso",
        "deltaH_iso_err",
    ]

    def __init__(self):
        self.S_inf = None
        self.S_inf_err = None
        self.deltaH_S_inf = None
        self.deltaH_S_inf_err = None
        self.deltaH_D = None
        self.deltaH_D_err = None
        self.deltaH_b = None
        self.deltaH_b_err = None
        self.c_iso = None
        self.deltaH_iso = None
        self.deltaH_iso_err = None

    def __str__(self):
        lines = [f"{self.__class__.__name__}("]
        for attr in self.__slots__:
            value = getattr(self, attr)
            if isinstance(value, np.ndarray):
                array_str = np.array2string(
                    value,
                    precision=3,
                    suppress_small=True,
                    threshold=10,
                    separator=", ",
                    max_line_width=80,
                )
                lines.append(f"  {attr} = {array_str}")
            else:
                lines.append(f"  {attr} = {str(value)}")
        lines.append(")")
        return "\n".join(lines)


def save_gas_class(gas, filename):
    """Saves a Gas object to a .pkl file.

    Args:
        gas: An instance of the Gas class.
        filename: Name of the .pkl file to save the Gas instance to.

    Returns:
        None
    """

    with open(filename, "wb") as f:
        pickle.dump(gas, f, protocol=pickle.HIGHEST_PROTOCOL)
    # print("Pickling successful")
    # print("--------------------------------------------------------------")


def load_gas_class(filename):
    """Loads a Gas object from a .pkl file.

    Args:
        filename: Name of the .pkl file to load the Gas instance from.

        Returns:
            A Gas instance loaded from the .pkl file.
    """

    with open(filename, "rb") as f:
        return pickle.load(f)


def base_loss(x, p, c, cerr):
    """Defines the loss function for TODO
    Args:

    Returns:

    """
    ssr = 0

    # initializing variables to solve for
    kD = x[0]
    CH = x[1]
    b = x[2]

    calc = kD * p + CH * b * p / (1 + b * p)

    ssr = np.sum(((calc - c) ** 2) / (cerr**2))

    return ssr


def parameter_hints(gas):
    # TODO

    c = gas.c

    calculate_fugacity(gas)

    cerr_vec = gas.c_err

    if gas.f is not None:
        p_vec = gas.f
    else:
        p_vec = gas.p

    T_vec = gas.temp

    x0 = [4.0, 50.0, 1.0]

    for i, p_i in enumerate(p_vec):
        c_i = c[i]
        cerr_i = cerr_vec[i]

        base_loss(x0, p_i, c_i, cerr_i)

        result = minimize(base_loss, x0, args=(p_i, c_i, cerr_i))

        print(result.x)
    # read in individual arrays

    # optimize DMS

    # save values

    # run regressions
    return


def calculate_fugacity(gas):
    """Calculates the fugacity from Gas.p.

    If Gas.f is None, it will attempt to calculate the fugacity using
        (in order of priority):
    (1) user-supplied Virial coefficients
    (2) user-supplied Peng-Robinson coefficients
    (3) built-in Virial coefficients (Gas.formula must be specified)
    (4) built-in Peng-Robinson coefficients (Gas.formula must be specified)

    Args:
        gas: An instance of the Gas class.

    Returns:
        None
    """
    print("------------------------Fugacity Check------------------------")
    if gas.f is not None:
        print("Fugacity data supplied in Gas.f by user")

    elif gas.virial_coeff:
        print("Calculating fugacity using user-supplied Virial coefficients")
        fugacity.virial_eos(gas)

    elif gas.pr_coeff:
        print("Calculating fugacity using user-supplied Peng-Robinson " "coefficients")
        fugacity.peng_robinson_eos(gas)

    elif gas.formula in fugacity.virial_coeff:
        print("Calculating fugacity using built-in Virial EoS")
        fugacity.virial_eos(gas)

    elif gas.formula in fugacity.pr_coeff:
        print("Calculating fugacity using built-in Peng-Robinson EoS")
        fugacity.peng_robinson_eos(gas)

    elif gas.f is None and gas.formula is None:
        pyDMS.warning_in_orange(
            "Fugacity data and gas unspecified.\nFitting will be performed " "with pressure data."
        )

    elif gas.formula not in fugacity.virial_coeff and gas.formula not in fugacity.pr_coeff:
        pyDMS.warning_in_orange(
            "Fugacity data unspecified and Gas.formula not found in built-in "
            "data.\nSupply Gas.virial_coeff or Gas.pr_coeff to calculate"
            "fugacity.\nFitting will be performed with pressure data."
        )

    else:
        pyDMS.warning_in_orange(
            "Something went wrong trying to calculate fugacity.\n Will use " "gas.p"
        )


def LFER_loss(x, gas, loss="chi2"):
    """Defines the loss function for the LFER optimization.

    Args:
        x: An array of the variables to be optimized.
        gas: An instance of the Gas class.
        loss: The loss function to use. Currently, only 'chi2' is supported.

    Returns:
        The loss function result as a number.
    """

    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.c_err
    if gas.f is not None:
        p_vec = gas.f
    else:
        p_vec = gas.p
    T_vec = gas.temp

    # initializing error array
    ssr = np.zeros(len(c_vec))

    # initializing variables to solve for
    kd0 = x[0]
    dHD = x[1]
    b0 = x[2]
    dHb = x[3]

    ch_vec = x[4:]

    # looping through each concentration vector
    for i, _ in enumerate(c_vec):

        # extracting parameters for relevant concentrations
        T = T_vec[i]
        ch = ch_vec[i]
        p = p_vec[i]
        c = c_vec[i]
        cerr = cerr_vec[i]

        # defining Linear Free Energy Relationships (LFERS)
        kd = kd0 * np.exp(-dHD * 1000 / (8.314 * T))
        b = b0 * np.exp(-dHb * 1000 / (8.314 * T))

        # solving the DMS model
        calc = kd * p + ch * b * p / (1 + b * p)

        if loss == "chi2":
            # using the chi-squared error as the loss function
            ssr[i] = np.sum(((calc - c) ** 2) / (cerr**2))
        else:
            return pyDMS.error_in_red("Loss function not recognized")
    # Sum of all the errors as metric
    out = np.sum(ssr)

    # TODO: see if scaling is really necessary (yes it is I think)
    return out / 1000


def vH_loss(x, gas, loss="chi2"):
    """Defines the loss function for the van't Hoff optimization.

    Args:
        x: An array of the variables to be optimized.
        gas: An instance of the Gas class with data from LFER_loss populated.
        loss: The loss function to use. Currently, only 'chi2' is supported.

    Returns:
        The final value of the loss function.
    """

    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.c_err
    if gas.f is not None:
        p_vec = gas.f
    else:
        p_vec = gas.p
    T_vec = gas.temp

    # initializing error array
    ssr = np.zeros(len(c_vec))

    # initializing variables to solve for
    dHD = x[0]
    dHb = x[1]
    ch_vec = x[2:]

    # looping through each concentration vector
    for i, _ in enumerate(c_vec):

        # extracting parameters for relevant concentrations
        T = T_vec[i]
        ch = ch_vec[i]
        p = p_vec[i]
        c = c_vec[i]
        cerr = cerr_vec[i]
        LFE_params = gas.LFER.fit
        a_kd0 = LFE_params[0]
        b_kd0 = LFE_params[1]
        a_b0 = LFE_params[2]
        b_b0 = LFE_params[3]

        # defining van't Hoff expressions
        kd0 = np.exp((dHD - b_kd0) / a_kd0)
        b0 = np.exp((dHb - b_b0) / a_b0)
        kd = kd0 * np.exp(-dHD * 1000 / (8.314 * T))
        b = b0 * np.exp(-dHb * 1000 / (8.314 * T))

        # solving the DMS model
        calc = kd * p + ch * b * p / (1 + b * p)

        if loss == "chi2":
            # using the chi-squared error as the loss function
            ssr[i] = np.sum(((calc - c) ** 2) / (cerr**2))
        else:
            return pyDMS.error_in_red("Loss function not recognized")

    # Sum of all the errors as metric
    out2 = np.sum(ssr)

    return out2


def is_outlier(arr):
    """Finds otuliers in an array.

    Utilizes the definition of isoutlier() from MATLAB
    (see https://www.mathworks.com/help/matlab/ref/isoutlier.html).

    Args:
        arr: An array to find outliers in

    Returns:
        An array with a Boolean mask for whether an entry is an outlier
            (TRUE) or not (FALSE).

        For example:

            if arr is [0.1, 0.12, 0.11, 100], it will return
                [False, False, False, True]
    """

    median = np.median(arr)
    mad = median_abs_deviation(arr)
    c = -1 / (np.sqrt(2) * erfcinv(3 / 2))  # MATLAB's scaling factor
    scaled_mad = c * mad

    return np.abs(arr - median) > 3 * scaled_mad  # Boolean mask for outliers


def hess(gas, soln):
    """Solves the analytical Hessian for the vH_loss loss function

    Args:
        gas: An instance of the Gas class
        soln: The object from a scipy.minimize call with loss function vH_loss

    Returns:
        A matrix with the Hessian for the result of vH_loss
    """

    # retriving gas data
    c_vec = gas.c
    cerr_vec = gas.c_err
    if gas.f is not None:
        p_vec = gas.f
    else:
        p_vec = gas.p
    T_vec = gas.temp

    LFE_params = gas.LFER.fit

    a_kd0 = LFE_params[0]
    b_kd0 = LFE_params[1]
    a_b0 = LFE_params[2]
    b_b0 = LFE_params[3]

    dHD = soln.x[0]
    dHb = soln.x[1]
    ch_vec = soln.x[2:]

    # initializing arrays for the Hessian
    dCdHD = [np.zeros_like(c_i) for c_i in c_vec]
    dCdHb = [np.zeros_like(c_i) for c_i in c_vec]
    dCdCH = [np.zeros_like(c_i) for c_i in c_vec]
    hessian = np.zeros((len(c_vec) + 2, len(c_vec) + 2)) / 1e6

    # looping through each concentration vector
    for i, _ in enumerate(c_vec):

        # extracting parameters for relevant concentrations
        T = T_vec[i]
        ch = ch_vec[i]
        p = p_vec[i]
        cerr = cerr_vec[i]

        # solving derivatives for the Hessian
        dCdHD[i] = (
            np.exp(-120.279 * dHD / T + (-b_kd0 + dHD) / a_kd0) * p * (-120.279 / T + 1 / a_kd0)
        )
        dCdHb[i] = (
            ch * np.exp(1 + (120.279 * dHb / T) + (b_b0 + dHb) / a_b0) * p * (T - 120.279 * a_b0)
        ) / ((np.exp(b_b0 / a_b0 + 120.279 * dHb / T) + np.exp(dHb / a_b0) * p) ** 2 * T * a_b0)
        dCdCH[i] = (1 + (np.exp((b_b0 - dHb) / a_b0 + 120.279 * dHb / T)) / p) ** -1

    # inserting values into Hessian
    hessian[0, 0] = np.sum(
        [np.sum(2 / cerr_vec[i] ** 2 * dCdHD[i] * dCdHD[i]) for i in range(len(c_vec))]
    )
    hessian[0, 1] = np.sum(
        [np.sum(2 / cerr_vec[i] ** 2 * dCdHD[i] * dCdHb[i]) for i in range(len(c_vec))]
    )
    hessian[1, 0] = hessian[0, 1]  # symmetric
    hessian[1, 1] = np.sum(
        [np.sum(2 / cerr_vec[i] ** 2 * dCdHb[i] * dCdHb[i]) for i in range(len(c_vec))]
    )

    for i in range(len(c_vec)):
        cerr = cerr_vec[i]
        hessian[0, i + 2] = np.sum(2 / cerr**2 * dCdHD[i] * dCdCH[i])
        hessian[i + 2, 0] = hessian[0, i + 2]
        hessian[1, i + 2] = np.sum(2 / cerr**2 * dCdHb[i] * dCdCH[i])
        hessian[i + 2, 1] = hessian[1, i + 2]
        hessian[i + 2, i + 2] = np.sum(2 / cerr**2 * dCdCH[i] * dCdCH[i])

    return hessian


def calc_LFEs(gas, settings=None):
    """Implements the LFER_loss function to optimize the DMS model based on
        LFER constraints

    Args:
        gas: An instance of the Gas class
        settings: A dictionary of settings for the optimization

    Returns:
        Data in the Gas.LFER subclass
    """

    # retriving gas data
    c = gas.c

    calculate_fugacity(gas)

    # cerr = gas.c_err
    # p = gas.p
    # T = gas.temp

    settings = gas.settings

    # finding settings
    if settings is None:
        settings = {}

    settings.setdefault("dHD_guess", [-1, -30])
    settings.setdefault("dHb_guess", [-1, -30])
    settings.setdefault("kD0_guess", [0.001, 0.01])
    settings.setdefault("b0_guess", [0.0001, 0.005])
    settings.setdefault("CH_guess", np.array([[0, 100] for _ in range(len(c))]))

    settings.setdefault("dHD_bounds", [-50, 0])
    settings.setdefault("dHb_bounds", [-50, 0])
    settings.setdefault("kD0_bounds", [0, None])
    settings.setdefault("b0_bounds", [0, None])
    settings.setdefault("CH_bounds", np.array([[0, 150] for _ in range(len(c))]))

    settings.setdefault("trials", 1000)
    settings.setdefault("solver_LFER", "SLSQP")
    settings.setdefault("ftol", 1e-7)
    settings.setdefault("xtol", 1e-7)
    settings.setdefault("gtol", 1e-7)
    settings.setdefault("maxiter_LFER", 1000)
    settings.setdefault("verbose", True)
    settings.setdefault("solver_verbose", False)
    settings.setdefault("seed", None)

    kd0_0_bnd = settings.get("kD0_guess")
    dHD0_0_bnd = settings.get("dHD_guess")
    b0_0_bnd = settings.get("b0_guess")
    dHb_0_bnd = settings.get("dHb_guess")
    ch_0_bnd = settings.get("CH_guess")

    kd0_b_solver = settings.get("kD0_bounds")
    dHD0_solver = settings.get("dHD_bounds")
    b0_solver = settings.get("b0_bounds")
    dHb_solver = settings.get("dHb_bounds")
    ch_solver = settings.get("CH_bounds")

    trials = settings.get("trials")
    solver = settings.get("solver_LFER")
    ftol_val = settings.get("ftol")
    xtol_val = settings.get("xtol")
    gtol_val = settings.get("gtol")
    maxiter_LFER = settings.get("maxiter_LFER")
    verbose = settings.get("verbose")
    solver_verbose = settings.get("solver_verbose")    
    rng_seed = settings.get("seed")

    # setting up the solver constraints: A*x <= b
    # (* consider allowing this to be turned off?)
    A_con = np.zeros((len(c) - 1, len(c) + 4))
    b_con = np.zeros(len(c) - 1)

    # setting up the constraint so that each CH' is lower than the previous
    for i in range(len(c) - 1):
        first_index = i + 4
        second_index = i + 5
        A_con[i, first_index] = -1
        A_con[i, second_index] = 1

    def linear_constraint(x):
        return b_con - (A_con @ x)

    # checks to ensure proper data formatting
    # if np.shape(c) != np.shape(cerr) != np.shape(p):
    #    pyDMS.error_in_red('c, cerr, and p dimensions do not match')

    # elif len(c) != len(T):
    #    pyDMS.error_in_red('T does not match c')

    # setting up arrays for storing minimization results
    nOptVars = 4 + len(c)  # no. of optimization variables
    dms = np.zeros((trials, nOptVars))  # results *
    res = np.zeros(trials)  # func(x)
    flag = np.zeros(trials)  # optimization flags

    # finding the range of each bound
    kdval = kd0_0_bnd[1] - kd0_0_bnd[0]
    dHDval = dHD0_0_bnd[1] - dHD0_0_bnd[0]
    b0val = b0_0_bnd[1] - b0_0_bnd[0]
    dHbval = dHb_0_bnd[1] - dHb_0_bnd[0]

    ch0_vals = [bnd[1] - bnd[0] for bnd in ch_0_bnd]

    ch_0 = np.zeros(len(ch0_vals))  # holding ch_0 guesses

    # printing bounds to search through (TODO: consider making it more table-like)
    if verbose:
        print("---------------------LFER Initial Guesses---------------------")
        print(f"kd0_0: {kd0_0_bnd}")
        print(f"dHd0_0: {dHD0_0_bnd}")
        print(f"b0_0: {b0_0_bnd}")
        print(f"dHb0_0: {dHb_0_bnd}")
        for i, ch_0_bnds in enumerate(ch_0_bnd):
            print(f"C_H'{i}: {ch_0_bnds}")
        print("---------------------LFER Solver Bounds-----------------------")
        print(f"kd0_0: {kd0_b_solver}")
        print(f"dHd0_0: {dHD0_solver}")
        print(f"b0_0: {b0_solver}")
        print(f"dHb0_0: {dHb_solver}")
        for i, ch_0_bnds in enumerate(ch_solver):
            print(f"C_H'{i}: {ch_0_bnds}")
        print("--------------------------------------------------------------")

    # initializing random number generator
    rng = np.random.default_rng(rng_seed)

    for j in range(trials):

        if verbose and np.mod(j, 50) == 0:
            print(f"LFER trial: {j}/{trials}")

        # picking a random number in each bound
        kd0_0 = kd0_0_bnd[0] + kdval * rng.random()
        dHD_0 = dHD0_0_bnd[0] + dHDval * rng.random()
        b0_0 = b0_0_bnd[0] + b0val * rng.random()
        dHb_0 = dHb_0_bnd[0] + dHbval * rng.random()

        for i, bounds in enumerate(ch_0_bnd):
            ch_0[i] = bounds[0] + ch0_vals[i] * rng.random()

        # setting up initial guess for solver
        x0 = np.concatenate([[kd0_0, dHD_0, b0_0, dHb_0], ch_0])

        # Setting bounds for solver
        solver_bounds = [
            (kd0_b_solver[0], kd0_b_solver[1]),  # kd0_0
            (dHD0_solver[0], dHD0_solver[1]),  # dHD_0
            (b0_solver[0], b0_solver[1]),  # b0_0
            (dHb_solver[0], dHb_solver[1]),
        ]  # dHb_0

        solver_bounds.extend([(ch_s[0], ch_s[1]) for ch_s in ch_solver])

        # optimizing with the 'trust-constr' algorithm
        if solver == "trust-constr":

            options = {
                "maxiter": maxiter_LFER,  # maximum iterations
                "xtol": xtol_val,  # variable tolerance
                "gtol": gtol_val,  # function tolerance
                "disp": solver_verbose,  # display optimization progress
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
                method="trust-constr",
                bounds=solver_bounds,
                constraints={"type": "ineq", "fun": linear_constraint},
                options=options,
            )

        # optimizing with the 'SLQSP' algorithm
        elif solver == "SLSQP":

            options = {
                "maxiter": maxiter_LFER,  # max iterations
                "ftol": ftol_val,  # function tolerance
                "disp": solver_verbose,  # display iteration results
            }

            # finding initial guess
            LFER_loss(x0, gas)

            # running optimization
            result = minimize(
                LFER_loss,
                x0,
                args=(gas),
                method="SLSQP",
                bounds=solver_bounds,
                constraints={"type": "ineq", "fun": linear_constraint},
                options=options,
            )

        # collecting optimization results
        flag[j] = result.status
        res[j] = result.fun
        dms[j, :] = result.x

        transposed_dms = np.transpose(dms)

    # FINDING OUTLIERS
    # Copy of the transposed results matrix
    par_with_outliers = np.copy(transposed_dms)

    # For plotting purposes
    with warnings.catch_warnings(record=True) as _:
        warnings.simplefilter("always")
        log_kd0_out = np.log(par_with_outliers[0])
        deltaHd_out = par_with_outliers[1]
        log_b0_out = np.log(par_with_outliers[2])
        deltaHb_out = par_with_outliers[3]

    ####################################
    # Old method does not work well to find outliers around a regression
    # keeping in case we ever want to retry it
    # outlier_track = np.zeros(nOptVars)
    # par_outliers_removed = np.copy(transposed_dms)

    # finding outliers for each parameter
    # for i in range(nOptVars):

    #    par = transposed_dms[i] # parameter of interest
    #    par_no_outliers = np.copy(par)

    # TF = is_outlier(par) # finding if iteration parameter is an outlier

    # par_no_outliers[TF] = 0
    # par_outliers_removed[i][TF] = 0
    # par = par[~TF]

    # outliers = np.where(TF == 1)[0]  # indices where outliers are present
    # * removed an if statement and should still work: keep an eye an it
    # outlier_track[i] = len(outliers)

    # avg_dms[i] = np.mean(par) # finding average of each cleaned parameter
    # std_dev[i] = np.std(par)
    # num_par_final[i] = len(par)

    # *TO DO: NEED TO EXAMINE THIS SOME MORE
    # par_outliers_removed = np.transpose(par_outliers_removed)
    # par_outliers_removed = par_outliers_removed[~np.any(
    #     par_outliers_removed == 0, axis=1)]
    # par_outliers_removed = np.transpose(par_outliers_removed)

    # FINDING LFER PARAMETERS
    # deltaHD = aD*ln(kd0) + bD
    # log_kd0 = np.log(par_outliers_removed[0])
    # deltaHd = par_outliers_removed[1]
    # log_b0 = np.log(par_outliers_removed[2])
    # deltaHb = par_outliers_removed[3]

    # ensuring an intercept is fitted as well
    # log_kd0_with_const = sm.add_constant(log_kd0)
    # log_b0_with_const = sm.add_constant(log_b0)

    # fitting linear regressions
    # deltaHd_model = sm.OLS(deltaHd, log_kd0_with_const).fit()
    # deltaHb_model = sm.OLS(deltaHb, log_b0_with_const).fit()

    # int_kd, slope_kd = deltaHd_model.params # int_kd = aD, slope_kD = bD
    # int_kd_err, slope_kd_err = deltaHd_model.bse
    # int_b, slope_b = deltaHb_model.params # int_b= ab, slope_b = bb
    # int_b_err, slope_b_err = deltaHb_model.bse # errors in int_b and slope_b

    # out = [slope_kd, int_kd, slope_b, int_b] # collecting slopes, intercepts
    # SE = [slope_kd_err, int_kd_err, slope_b_err, int_b_err]

    # collecting data
    # * keep an eye on this. Changed from len(c) to 4
    # pars = np.zeros((len(log_kd0), 4))
    # pars[:,0] = log_kd0
    # pars[:,1] = deltaHd
    # pars[:,2] = log_b0
    # pars[:,3] = deltaHb
    #######################################################################

    valid_mask = (transposed_dms[0] > 0) & (transposed_dms[2] > 0)

    log_kd0 = np.log(transposed_dms[0][valid_mask])
    deltaHd = transposed_dms[1][valid_mask]
    log_b0 = np.log(transposed_dms[2][valid_mask])
    deltaHb = transposed_dms[3][valid_mask]

    X_kd = log_kd0.reshape(-1, 1)
    X_b = log_b0.reshape(-1, 1)
    ransac_kd = RANSACRegressor(
        LinearRegression(), residual_threshold=0.25, max_trials=1000, random_state=0
    )
    ransac_kd.fit(X_kd, deltaHd)

    ransac_b = RANSACRegressor(
        LinearRegression(), residual_threshold=0.25, max_trials=1000, random_state=0
    )
    ransac_b.fit(X_b, deltaHb)

    # Get coefficients
    slope_kd = ransac_kd.estimator_.coef_[0]
    int_kd = ransac_kd.estimator_.intercept_

    slope_b = ransac_b.estimator_.coef_[0]
    int_b = ransac_b.estimator_.intercept_

    # Store output (errors need manual calc)
    out = [slope_kd, int_kd, slope_b, int_b]

    inliers_kd = ransac_kd.inlier_mask_
    inliers_b = ransac_b.inlier_mask_

    X_kd_in = sm.add_constant(X_kd[inliers_kd])
    X_b_in = sm.add_constant(X_b[inliers_b])

    kd_fit = sm.OLS(deltaHd[inliers_kd], X_kd_in).fit()
    b_fit = sm.OLS(deltaHb[inliers_b], X_b_in).fit()

    SE = [kd_fit.bse[1], kd_fit.bse[0], b_fit.bse[1], b_fit.bse[0]]

    # Initialize full-length array with NaNs
    pars = np.full((len(log_kd0), 4), np.nan)

    # Fill in values for inliers in each model
    pars[inliers_kd, 0] = log_kd0[inliers_kd]
    pars[inliers_kd, 1] = deltaHd[inliers_kd]

    pars[inliers_b, 2] = log_b0[inliers_b]
    pars[inliers_b, 3] = deltaHb[inliers_b]

    # collecting data with outliers for plotting
    # *switched from len(c) to 4
    pars_outliers = np.zeros((len(log_kd0_out), 4))
    pars_outliers[:, 0] = log_kd0_out
    pars_outliers[:, 1] = deltaHd_out
    pars_outliers[:, 2] = log_b0_out
    pars_outliers[:, 3] = deltaHb_out

    gas.LFER.fit = out
    gas.LFER.fit_err = SE
    gas.LFER.out = pars
    gas.LFER.out_outliers = pars_outliers

    gas.settings = settings

    return gas


def calc_params(gas):
    """Implements the vH_loss function to optimize the DMS model based on
        van't Hoff constraints.
    Args:
        gas: An instance of the Gas class.

    Returns:
        Data in the Gas.vH class
    """

    # pulling in data
    c = gas.c
    # cerr = gas.c_err
    # p = gas.p
    T = gas.temp
    LFE_params = gas.LFER.fit

    settings = gas.settings

    # finding settings
    if settings is None:
        settings = {}

    settings.setdefault("dHD_guess", [-1, -30])
    settings.setdefault("dHb_guess", [-1, -30])
    settings.setdefault("CH_guess", np.array([[0, 100] for _ in range(len(c))]))

    settings.setdefault("dHD_bounds", [-50, 0])
    settings.setdefault("dHb_bounds", [-50, 0])
    settings.setdefault("CH_bounds", np.array([[0, 150] for _ in range(len(c))]))

    settings.setdefault("trials", 1000)
    settings.setdefault("solver_vH", "SLSQP")
    settings.setdefault("ftol", 1e-7)
    settings.setdefault("xtol", 1e-7)
    settings.setdefault("gtol", 1e-7)
    settings.setdefault("maxiter_vH", 1000)
    settings.setdefault("verbose", True)
    settings.setdefault("solver_verbose", False)
    settings.setdefault("seed", None)

    dHD0_0_bnd = settings.get("dHD_guess")
    dHb_0_bnd = settings.get("dHb_guess")
    ch_0_bnd = settings.get("CH_guess")

    dHD0_solver = settings.get("dHD_bounds")
    dHb_solver = settings.get("dHb_bounds")
    ch_solver = settings.get("CH_bounds")

    trials = settings.get("trials")
    solver = settings.get("solver_vH")
    ftol_val = settings.get("xtol")
    xtol_val = settings.get("xtol")
    gtol_val = settings.get("gtol")
    maxiter_vH = settings.get("maxiter_vH", 1000)
    verbose = settings.get("verbose")
    solver_verbose = settings.get("solver_verbose")
    rng_seed = settings.get("seed")

    # setting up the solver constraints: A*x <= b
    A_con = np.zeros((len(c) - 1, len(c) + 2))
    b_con = np.zeros(len(c) - 1)

    # if np.shape(c) != np.shape(cerr) != np.shape(p):
    #    pyDMS.error_in_red('c, cerr, and p dimensions do not match')

    # elif len(c) != len(T):
    #    pyDMS.error_in_red('*some error')

    nOptVars = 2 + len(c)  # number of optimization variables
    dms = np.zeros((trials, nOptVars))  # creating a matrix for holding results
    res = np.zeros(trials)
    flag = np.zeros(trials)
    avg_dms = np.zeros(nOptVars)
    std_dev = np.zeros(nOptVars)
    num_par_final = np.zeros(nOptVars)
    hessian_matrix = np.zeros((len(c) + 2, len(c) + 2, trials))

    # finding the range of each bound
    dHDval = dHD0_0_bnd[1] - dHD0_0_bnd[0]
    dHbval = dHb_0_bnd[1] - dHb_0_bnd[0]
    ch0_vals = [bnd[1] - bnd[0] for bnd in ch_0_bnd]

    ch_0 = np.zeros(len(ch0_vals))  # holding ch_0 guesses

    # printing bounds to search through
    if verbose:
        print("------------------van't Hoff Initial Guesses------------------")
        print(f"dHD0_0: {dHD0_0_bnd}")
        print(f"dHb0_0: {dHb_0_bnd}")
        for i, ch_0_bnds in enumerate(ch_0_bnd):
            print(f"C_H'{i}: {ch_0_bnds}")
        print("-------------------van't Hoff Solver Bounds-------------------")
        print(f"dHD0_0: {dHD0_solver}")
        print(f"dHb0_0: {dHb_solver}")
        for i, ch_0_bnds in enumerate(ch_solver):
            print(f"C_H'{i}: {ch_0_bnds}")
        print("--------------------------------------------------------------")

    # initializing random number generator
    rng = np.random.default_rng(rng_seed)

    for j in range(trials):

        if verbose and np.mod(j, 50) == 0:
            print(f"van't Hoff trial: {j}/{trials}")

        # picking a random number in each bound
        dHD_0 = dHD0_0_bnd[0] + dHDval * rng.random()
        dHb_0 = dHb_0_bnd[0] + dHbval * rng.random()

        for i, bounds in enumerate(ch_0_bnd):
            ch_0[i] = bounds[0] + ch0_vals[i] * rng.random()

        # setting up initial guess for solver
        x1 = np.concatenate([[dHD_0, dHb_0], ch_0])

        def linear_constraint(x):
            return b_con - (A_con @ x)

        # * consider storing the inital guess inside initial_g
        for i in range(len(c) - 1):
            first_index = i + 2
            second_index = i + 3
            A_con[i, first_index] = -1
            A_con[i, second_index] = 1

        # Setting bounds for solver. Generally less rigorous than guesses
        solver_bounds = [(dHD0_solver[0], dHD0_solver[1]), (dHb_solver[0], dHb_solver[1])]

        solver_bounds.extend([(ch_s[0], ch_s[1]) for ch_s in ch_solver])

        # optimizing with the 'trust-constr' algorithm
        if solver == "trust-constr":

            options = {
                "maxiter": maxiter_vH,  # maximum iterations
                "xtol": xtol_val,  # variable tolerance
                "gtol": gtol_val,  # function tolerance
                "disp": solver_verbose,  # display optimization progress
            }

            # finding initial guess
            vH_loss(x1, gas)

            # flag = np.zeros((1,j))
            if not verbose:
                warnings.filterwarnings("ignore", message="delta_grad == 0.0")

            # running optimization
            result = minimize(
                vH_loss,
                x1,
                args=(gas),
                method="trust-constr",
                bounds=solver_bounds,
                constraints={"type": "ineq", "fun": linear_constraint},
                options=options,
            )

        # optimizing with the 'SLQSP' algorithm
        elif solver == "SLSQP":

            options = {
                "maxiter": maxiter_vH,  # max iterations
                "ftol": ftol_val,  # function tolerance
                "disp": solver_verbose,  # display iteration results
            }

            # finding initial guess
            vH_loss(x1, gas)

            # * create flag array
            # flag = np.zeros((1,j))

            # running optimization
            result = minimize(
                vH_loss,
                x1,
                args=(gas,),
                method="SLSQP",
                bounds=solver_bounds,
                constraints={"type": "ineq", "fun": linear_constraint},
                options=options,
            )
        # print("STATUS")
        # print(result.status)
        hessian_matrix[:, :, j] = hess(gas, result)
        # print("=======HESSIAN========")
        # print(hessian_matrix[:, :, j])
        flag[j] = result.status
        res[j] = result.fun
        dms[j, :] = result.x

        transposed_dms = np.transpose(dms)
        par2_outliers_removed = np.copy(transposed_dms)

    outlier_track = np.zeros(nOptVars)

    n_trials = transposed_dms.shape[1]  # == dms.shape[0] FOR HESS
    keep = np.ones(n_trials, dtype=bool)  # start by keeping all trials FOR HESS

    for i in range(len(x1)):

        par2 = transposed_dms[i]

        TF = is_outlier(par2)

        keep &= ~TF  # drop any trial that was an outlier for parameter i THIS IS FOR HESS

        par2 = par2[~TF]

        # setting outliers to zero in complete array as well
        par2_outliers_removed[i][TF] = 0

        # indices where outliers are True
        outliers = np.where(TF == 1)[0]

        # store number of outliers for each parameter
        outlier_track[i] = len(outliers)

        # finding average of each cleaned parameter
        avg_dms[i] = np.mean(par2)
        # finding standard deviation of each cleaned parameter
        std_dev[i] = np.std(par2)
        # number of remaining parameter iterations post-cleaning
        num_par_final[i] = len(par2)

    par2_outliers_removed = np.transpose(par2_outliers_removed)
    par2_outliers_removed = par2_outliers_removed[~np.any(par2_outliers_removed == 0, axis=1)]

    valid_idx = np.where(keep)[0]

    # keep only fully valid trials (no outlier in any parameter)
    # dms_valid = dms[valid_idx, :] # HESS
    hessian_matrix_valid = hessian_matrix[:, :, valid_idx]  # HESS
    hessian_matrix_valid = 0.5 * (hessian_matrix_valid + hessian_matrix_valid.swapaxes(0, 1))
    hessian_mean = np.mean(hessian_matrix_valid, axis=2)
    # res_valid = res[valid_idx] # HESS
    # flag_valid = flag[valid_idx] # HESS

    # n_valid = int(keep.sum()) # HESS
    # valid_idx = np.flatnonzero(keep) # HESS

    # print("Number of valid optimizations:", n_valid) # HESS
    # LOOK INTO RUBIN FACTOR
    # print(hessian_sum)

    dHD_f = avg_dms[0]
    dHb_f = avg_dms[1]
    ch_vec_f = avg_dms[2:]
    a_kd0_f = LFE_params[0]
    b_kd0_f = LFE_params[1]
    a_b0_f = LFE_params[2]
    b_b0_f = LFE_params[3]
    kd0_f = np.exp((dHD_f - b_kd0_f) / a_kd0_f)
    b0_f = np.exp((dHb_f - b_b0_f) / a_b0_f)

    kd_f = np.zeros(len(T))
    b_f = np.zeros(len(T))

    for i, temp in enumerate(T):
        kd_f[i] = kd0_f * np.exp(-dHD_f * 1000 / (8.314 * temp))
        b_f[i] = b0_f * np.exp(-dHb_f * 1000 / (8.314 * temp))

    gas.kD = kd_f
    gas.b = b_f
    gas.CH = ch_vec_f
    gas.vH.out_outliers = dms
    gas.vH.out = par2_outliers_removed
    gas.vH.avg_dms = avg_dms
    gas.vH.residuals = res
    # gas.vH.hessian_matrix = hessian_matrix
    gas.vH.hessian_matrix = hessian_mean

    gas.settings = settings

    gas.analysis.deltaH_D = [avg_dms[0], None]
    gas.analysis.deltaH_b = [avg_dms[1], None]

    return gas


def chi2_error_fit(gas):
    """Calculates uncertainty via determination of the covariance materix.
    Args:
        gas: An instance of the Gas class.

    Returns:
        Data in the Gas and Gas.analysis classes
    """
    
    hessian_avg = gas.vH.hessian_matrix
    res = gas.vH.residuals

    ch_vec_f = gas.CH
    dHD_f = gas.analysis.deltaH_D[0]
    dHb_f = gas.analysis.deltaH_b[0]
    avg_dms = np.concatenate(([dHD_f, dHb_f], ch_vec_f))

    # mean of the Hessians across all trials
    # hessian_avg = np.mean(hessian_matrix, axis=2)
    # print(hessian_avg)

    hess_smoothed = 0.5 * (hessian_avg + np.transpose(hessian_avg))

    try:
        c, lower = cho_factor(hess_smoothed, check_finite=False)
    except Exception:
        lam = 1e-10 * np.trace(hess_smoothed) / hess_smoothed.shape[0]
        c, lower = c, lower = cho_factor(
            hess_smoothed + lam * np.eye(hess_smoothed.shape[0]), check_finite=False
        )

    p = hess_smoothed.shape[0]
    diag_hess_inv = np.empty(p)

    for m in range(p):
        e = np.zeros(p)
        e[m] = 1.0
        y = cho_solve((c, lower), e, check_finite=False)
        diag_hess_inv[m] = y[m]

    var = (1 / 2) * diag_hess_inv
    std_dev = np.sqrt(var)

    plusminus_1 = [std_dev[0], std_dev[0]]
    plusminus_2 = [std_dev[1], std_dev[1]]

    gas.vH.plusminus_1 = plusminus_1
    gas.vH.plusminus_2 = plusminus_2
    gas.CH_err = np.array(std_dev[2:])

    gas.analysis.deltaH_D_err = [plusminus_1[1], None]
    gas.analysis.deltaH_b_err = [plusminus_2[1], None]

    # LEGACY VERSION OF COLLECTING CONFIDENCE INTERVALS USING CHI^2 PERTURBATION
    # gas.analysis.deltaH_D_err = [plusminus_1[1], None]
    # gas.analysis.deltaH_b_err = [plusminus_2[1], None]
    # mean of the residuals across all trials
    # f = np.mean(res)

    # mean of the residuals shifted by the critical chi-square value
    # chi_crit = f + chi2.isf(0.05, len(avg_dms))

    # setting up arrays for error analysis
    # avg_dms1 = np.copy(avg_dms)
    # avg_dms2 = np.copy(avg_dms)
    # making sure arrays are not linked
    # avg_dms_ch = [np.copy(avg_dms) for _ in range(len(ch_vec_f))]

    # finding dHD and dHb +/- 1 away from optimal value
    # dHD_vec = np.linspace(avg_dms[0] - 1, avg_dms[0] + 1, 100)
    # dHb_vec = np.linspace(avg_dms[1] - 1, avg_dms[1] + 1, 100)

    # finding C_H' +/- 1 away from optimal value
    # ch_vec = [np.linspace(v - 1, v + 1, 100) for v in avg_dms[2:]]

    # setting up more arrays for error analysis
    # chi2_1 = np.zeros(len(dHD_vec))
    # chi2_2 = np.zeros(len(dHD_vec))
    # chi2_ch = [np.zeros(len(dHD_vec)) for _ in range(len(ch_vec_f))]

    # perturbing the parameters to create a chi-square distribution
    # for i, dHD_val in enumerate(dHD_vec):
    #    avg_dms1[0] = dHD_val
    #    avg_dms2[1] = dHb_vec[i]

    #    for j in range(len(ch_vec_f)):
    #        avg_dms_ch[j][j + 2] = ch_vec[j][i]

    # solving for the chi-square distributions for each parameter as a
    # Taylor series expansion
    #    chi2_1[i] = f + (1 / 2) * (
    #        (avg_dms1 - avg_dms) @ hessian_avg @ np.transpose((avg_dms1 - avg_dms))
    #    )
    #    chi2_2[i] = f + (1 / 2) * (
    #        (avg_dms2 - avg_dms) @ hessian_avg @ np.transpose((avg_dms2 - avg_dms))
    #    )

    #    for j in range(len(ch_vec_f)):
    #        chi2_ch[j][i] = f + (1 / 2) * (
    #            (avg_dms_ch[j] - avg_dms) @ hessian_avg @ np.transpose((avg_dms_ch[j] - avg_dms))
    #        )

    # fitting the chi-square distributions to parabolas
    # chi2_1_fit = poly.polyfit(dHD_vec, chi2_1, 2)
    # chi2_2_fit = poly.polyfit(dHb_vec, chi2_2, 2)
    # chi2_ch_fit = [np.zeros(3) for _ in range(len(ch_vec_f))]

    # for j in range(len(ch_vec_f)):
    #    chi2_ch_fit[j] = poly.polyfit(ch_vec[j], chi2_ch[j], 2)

    # shifting the fits by the chi-square minimum
    # chi2_1_fit[0] = chi2_1_fit[0] - chi_crit
    # chi2_2_fit[0] = chi2_2_fit[0] - chi_crit

    # for j in range(len(ch_vec_f)):
    #    chi2_ch_fit[j][0] = chi2_ch_fit[j][0] - chi_crit

    # finding where the distributions are equal to the crit. chi-square value
    # roots_1 = poly.polyroots(chi2_1_fit)
    # roots_2 = poly.polyroots(chi2_2_fit)

    # roots_ch = [np.zeros(2) for _ in range(len(ch_vec_f))]

    # for j in range(len(ch_vec_f)):
    #    roots_ch[j] = poly.polyroots(chi2_ch_fit[j])

    # plusminus_1 = roots_1 - avg_dms[0]  # error for dHD
    # plusminus_2 = roots_2 - avg_dms[1]  # error for dHb

    # plusminus_ch = [np.zeros(2) for _ in range(len(ch_vec_f))]

    # for j in range(len(ch_vec_f)):
    #    plusminus_ch[j] = roots_ch[j] - avg_dms[j + 2]

    # ch_err = np.array([arr[1] for arr in plusminus_ch])  # error for C_H'

    # gas.vH.plusminus_1 = plusminus_1
    # gas.vH.plusminus_2 = plusminus_2
    # gas.vH.plusminus_ch = plusminus_ch
    # gas.CH_err = ch_err

    # gas.analysis.deltaH_D_err = [plusminus_1[1], None]
    # gas.analysis.deltaH_b_err = [plusminus_2[1], None]

    return gas


def propogate_error(gas):
    """Determines error from the van't Hoff optimization and propagates
        the error to further parameters

    Args:
        gas: An instance of the Gas class with Gas.vH data populated

    Returns:
        Data in the Gas.vH subclass and Gas class
    """
    print("--------------------------------------------------------------")
    print("Optimization successful; starting error propagation")

    LFE_params = gas.LFER.fit
    LFE_error = gas.LFER.fit_err
    avg_dms = gas.vH.avg_dms
    T_vec = gas.temp
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
        delkddelT[i] = (1000 * dHD_f / (8.314 * T**2)) * np.exp(
            -1000 * dHD_f / (8.314 * T) + (dHD_f - b_kd0) / a_kd0
        )
        delkddel_a_kd0[i] = ((b_kd0 - dHD_f) / a_kd0**2) * np.exp(
            -1000 * dHD_f / (8.314 * T) + (dHD_f - b_kd0) / a_kd0
        )
        delkddel_b_kd0[i] = (-1 / a_kd0) * np.exp(
            -1000 * dHD_f / (8.314 * T) + (dHD_f - b_kd0) / a_kd0
        )
        delkddelHd[i] = (1 / a_kd0 - 1000 / 8.314 / T) * np.exp(
            -1000 * dHD_f / (8.314 * T) + (dHD_f - b_kd0) / a_kd0
        )
        kd_err[i] = np.sqrt(
            0.01**2 * delkddelT[i] ** 2
            + err_a_kd0**2 * delkddel_a_kd0[i] ** 2
            + err_b_kd0**2 * delkddel_b_kd0[i] ** 2
            + plusminus_1[1] ** 2 * delkddelHd[i] ** 2
        )

    b_err = np.zeros(len(T_vec))
    delbdel_a_b0 = np.zeros(len(T_vec))
    delbdel_b_b0 = np.zeros(len(T_vec))
    delbdelT = np.zeros(len(T_vec))
    delbdelHb = np.zeros(len(T_vec))

    for i, T in enumerate(T_vec):
        delbdelT[i] = (1000 * dHb_f / (8.314 * T**2)) * np.exp(
            -1000 * dHb_f / (8.314 * T) + (dHb_f - b_b0) / a_b0
        )
        delbdel_a_b0[i] = ((b_b0 - dHb_f) / a_b0**2) * np.exp(
            -1000 * dHb_f / (8.314 * T) + (dHb_f - b_b0) / a_b0
        )
        delbdel_b_b0[i] = (-1 / a_b0) * np.exp(-1000 * dHb_f / (8.314 * T) + (dHb_f - b_b0) / a_b0)
        delbdelHb[i] = (1 / a_b0 - 1000 / 8.314 / T) * np.exp(
            -1000 * dHb_f / (8.314 * T) + (dHb_f - b_b0) / a_b0
        )
        b_err[i] = np.sqrt(
            0.01**2 * delbdelT[i] ** 2
            + err_a_b0**2 * delbdel_a_b0[i] ** 2
            + err_b_b0**2 * delbdel_b_b0[i] ** 2
            + plusminus_2[1] ** 2 * delbdelHb[i] ** 2
        )

    gas.kD_err = kd_err
    gas.b_err = b_err

    return gas


def compute(gas, output="unnamed_file"):
    """Provides a wrapper to run the entire optimization procedure with a
        single function call.

    Args:
        gas: An instance of the Gas class.
        output: filenames of the PDF and .pkl files

    Returns:
        None
    """
    print(
        r"""
==============================================================
              ___  __  _______
   ___  __ __/ _ \/  |/  / __/      Copyright (C) 2025
  / _ \/ // / // / /|_/ /\ \        Massachusetts Institute
 / .__/\_, /____/_/  /_/___/        of Technology
/_/   /___/                    

Authors: B.C. Tapia, P.A. Dean, J.Y. Yeo, A.X. Wu, Z.P. Smith
Web: https://smithlab.mit.edu
=============================================================="""
    )
    calc_LFEs(gas)
    calc_params(gas)
    chi2_error_fit(gas)
    propogate_error(gas)
    print("Error propagation successful")
    print("--------------------------------------------------------------")
    print("Computing sorption energetics")
    evaluate.heat_of_sorption(gas)
    evaluate.isosteric_heat(gas)
    print("--------------------------------------------------------------")
    # vis.heat_of_sorption(gas)
    if output is not None:
        print(f"Pickling {output}.pkl")
        save_gas_class(gas, filename=output + ".pkl")
        print("Pickling successful")
        print("--------------------------------------------------------------")
        report.generate(gas, report_name=output + ".pdf")
        print("pyDMS successful!")
        print("=======================END OF PROGRAM========================")
