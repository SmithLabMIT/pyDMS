"""
pyDMS.visualize

Copyright 2025 Massachusetts Institute of Technology
Licensed under the MIT License
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from . import evaluate

small_txt = 9
medium_txt = 12
big_txt = 25

plt.rc("font", size=small_txt)
plt.rc("axes", titlesize=small_txt)
plt.rc("axes", labelsize=medium_txt)
plt.rc("xtick", labelsize=small_txt)
plt.rc("ytick", labelsize=small_txt)
plt.rc("legend", fontsize=small_txt)


def LFER(gas, outliers=False, show=False):
    """Plots the linear LFER fits

    report.LFER should be reserved specifically for printing in the report.
    For general plotting of the LFERs, use visualization.LFER

    Args:
        gas: An instance of the Gas class with Gas.LFER populated
        outliers: whether the plot should contain all outlier data
        show: whether the plot should be printed to the screen

    Returns:
        None
    """

    slope_kd, int_kd, slope_b, int_b = gas.LFER.fit
    pars = gas.LFER.out
    pars_outliers = gas.LFER.out_outliers

    mask_kd = ~np.isnan(pars[:, 0]) & ~np.isnan(pars[:, 1])
    mask_b = ~np.isnan(pars[:, 2]) & ~np.isnan(pars[:, 3])

    log_kd0 = pars[mask_kd, 0]
    deltaHd = pars[mask_kd, 1]
    log_b0 = pars[mask_b, 2]
    deltaHb = pars[mask_b, 3]

    log_kd0_out = pars_outliers[:, 0]
    deltaHd_out = pars_outliers[:, 1]
    log_b0_out = pars_outliers[:, 2]
    deltaHb_out = pars_outliers[:, 3]

    fig, ax = plt.subplots(1, 2, figsize=(6.2, 3))  # never more than two plots

    cmap = cm.viridis
    norm = mcolors.Normalize(vmin=0, vmax=10)

    def lin_fit(data, slope, ints):
        return slope * data + ints

    if outliers:
        ax[0].plot(log_kd0_out, deltaHd_out, "o", color="gray", alpha=0.7)
    ax[0].plot(log_kd0, deltaHd, "o", color="black")
    x_fit_kd = np.linspace(np.min(log_kd0), np.max(log_kd0), 100)
    ax[0].plot(
        x_fit_kd, lin_fit(x_fit_kd, slope_kd, int_kd), label="RANSAC fit", color=cmap(norm(1))
    )
    ax[0].set_xlabel(r"$\mathrm{ln}(k_{\mathrm{D},0})$")
    ax[0].set_ylabel(r"$\Delta H_\mathrm{D} \; (\mathrm{kJ \; mol^{-1}})$")

    if outliers:
        ax[1].plot(log_b0_out, deltaHb_out, "o", color="gray", alpha=0.7)
    ax[1].plot(log_b0, deltaHb, "o", color="black")
    x_fit_b = np.linspace(np.min(log_b0), np.max(log_b0), 100)
    ax[1].plot(x_fit_b, lin_fit(x_fit_b, slope_b, int_b), label="RANSAC fit", color=cmap(norm(5)))
    ax[1].set_xlabel(r"$\mathrm{ln}(\mathrm{b}_0)$")
    ax[1].set_ylabel(r"$\Delta H_\mathrm{b} \; (\mathrm{kJ \; mol^{-1}})$")

    plt.tight_layout()

    if show:
        plt.show()


def histograms(gas, show=False):
    """Plots and saves the histograms from the van't Hoff fits

    Args:
        gas: An instance of the Gas class with Gas.vH populated
        show: whether the plot should be printed to the screen

    Returns:
        None
    """

    dms = gas.vH.out_outliers
    dms_no_out = gas.vH.out
    transposed_dms = np.transpose(dms)
    transposed_dms_no_out = np.transpose(dms_no_out)
    temps = gas.temp
    n = transposed_dms.shape[0]
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, ax = plt.subplots(nrows, ncols, figsize=(10, 2.5 * nrows))
    ax = np.array(ax).reshape(nrows, ncols)

    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=n)

    for i in range(n):
        row, col = divmod(i, ncols)

        if i == 0:
            label = r"$\Delta H_{\mathrm{D}}-\overline{\Delta H_{\mathrm{D}}} \; (\mathrm{kJ \; mol^{-1}})$"
        elif i == 1:
            label = r"$\Delta H_\mathrm{b}-\overline{\Delta H_\mathrm{b}} \; (\mathrm{kJ \; mol^{-1}})$"
        else:
            tex = "$C_\\mathrm{H}^{\\prime}-\\overline{C_\\mathrm{H}^{\\prime}}$"
            label = f"({tex}) ({temps[i-2]} K)"

        ax[row, col].set_xlabel(label)
        ax[row, col].set_ylabel("Trials")
        color = cmap(norm(i))

        ax[row, col].hist(
            transposed_dms[i] - np.mean(transposed_dms[i]), color="gray", bins="sqrt", alpha=0.5
        )
        ax[row, col].hist(
            transposed_dms_no_out[i] - np.mean(transposed_dms_no_out[i]), color=color, bins="sqrt"
        )

    # Turn off unused axes
    total_axes = nrows * ncols
    for k in range(n, total_axes):
        r, c = divmod(k, ncols)
        ax[r, c].axis("off")

    plt.tight_layout()
    if show:
        plt.show()


def isotherms(gas, show=False):
    """Plots and saves the sorption isotherms with DMS parameters

    Args:
        gas: An instance of the Gas class with Gas populated
        show: Whether the plot should be printed to the screen

    Returns:
        None
    """

    T = gas.temp
    n = len(T)

    ncols = 2
    nrows = int(np.ceil(n / ncols))

    fig, ax = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3 * nrows))
    ax = np.array(ax).reshape(nrows, ncols)

    c = gas.c
    cerr = gas.c_err
    if gas.f is not None:
        p_vec = gas.f
        x_label = "f"
    else:
        p_vec = gas.p
        x_label = "P"

    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=n)

    for i, temp in enumerate(T):
        row, col = divmod(i, ncols)

        press, c_model, c_model_err = evaluate.isotherm(gas, i)
        linecolor = cmap(norm(i))

        ax[row, col].set_ylabel(r"$C \; \mathrm{(cm^3_{STP} \; cm^{-3}_{pol})}$")
        ax[row, col].set_xlabel(rf"${x_label} \; \mathrm{{(atm)}}$")
        ax[row, col].plot(press, c_model, color=linecolor, label=f"{temp} K")
        ax[row, col].fill_between(
            press, c_model - c_model_err / 2, c_model + c_model_err / 2, color=linecolor, alpha=0.3
        )

        ax[row, col].plot(p_vec[i], c[i], "o", color="black")
        ax[row, col].errorbar(p_vec[i], c[i], yerr=cerr[i], fmt="none", color="black")
        ax[row, col].legend(frameon=False)

    # Turn off unused axes
    for k in range(n, nrows * ncols):
        r, c_ = divmod(k, ncols)
        ax[r, c_].axis("off")

    plt.tight_layout()
    if show:
        plt.show()


def heat_of_sorption(gas, show=False):
    """Plots and saves the heats of sorption

    report.heat_of_sorption should be reserved specifically for printing in the
        report.
    For general plotting of the isotherms, use visualization.heat_of_sorption

    Args:
        gas: An instance of the Gas class with Gas populated

    Returns:
        The name of the file that was saved
    """

    temp = gas.temp

    inv_RT = 1 / (0.008314 * temp)

    S_inf = gas.analysis.S_inf
    k_D = gas.kD
    b = gas.b

    # S_inf_err = gas.analysis.S_inf_err
    # k_D_err = gas.kD_err
    # b_err = gas.b_err

    S_inf_0 = gas.analysis.deltaH_S_inf[1]
    deltaH_S_inf = gas.analysis.deltaH_S_inf[0]

    k_D_0 = gas.analysis.deltaH_D[1]
    deltaH_D = gas.analysis.deltaH_D[0]

    b_0 = gas.analysis.deltaH_b[1]
    deltaH_b = gas.analysis.deltaH_b[0]

    ln_S_inf = np.log(S_inf)
    ln_k_D = np.log(k_D)
    ln_b = np.log(b)

    ln_S_inf_0 = np.log(S_inf_0)
    ln_k_D_0 = np.log(k_D_0)
    ln_b_0 = np.log(b_0)

    fig, ax = plt.subplots(1, 3, figsize=(9, 3))

    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=5)

    def lin_fit(data, slope, ints):
        return slope * data + ints

    # print(deltaH_S_inf)
    ax[0].plot(inv_RT, ln_S_inf, "o", color="black")
    ax[0].plot(inv_RT, lin_fit(inv_RT, -deltaH_S_inf, ln_S_inf_0), label="OLS", color=cmap(norm(1)))
    ax[0].set_ylabel(r"$\mathrm{ln}(S_\infty)$")
    ax[0].set_xlabel(r"$(RT)^{-1} \; (\mathrm{mol \; kJ^{-1}})$")

    ax[1].plot(inv_RT, ln_k_D, "o", color="black")
    ax[1].plot(inv_RT, lin_fit(inv_RT, -deltaH_D, ln_k_D_0), label="OLS", color=cmap(norm(2)))
    ax[1].set_ylabel(r"$\mathrm{ln}(k_\mathrm{D})$")
    ax[1].set_xlabel(r"$(RT)^{-1} \; (\mathrm{mol \; kJ^{-1}})$")

    ax[2].plot(inv_RT, ln_b, "o", color="black")
    ax[2].plot(inv_RT, lin_fit(inv_RT, -deltaH_b, ln_b_0), label="OLS", color=cmap(norm(3)))
    ax[2].set_ylabel(r"$\mathrm{ln}(b)$")
    ax[2].set_xlabel(r"$(RT)^{-1} \; (\mathrm{mol \; kJ^{-1}})$")

    plt.tight_layout()

    if show:
        plt.show()
