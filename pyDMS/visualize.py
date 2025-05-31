'''
pyDMS.visualize
       

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from . import evaluate

small_txt = 9
medium_txt = 12
big_txt = 25

plt.rc('font', size=small_txt)
plt.rc('axes', titlesize=small_txt)
plt.rc('axes', labelsize=medium_txt)
plt.rc('xtick', labelsize=small_txt)
plt.rc('ytick', labelsize=small_txt)
plt.rc('legend', fontsize=small_txt)

def LFER(gas, outliers=False, show=False):
    '''Plots the linear LFER fits

    report.LFER should be reserved specifically for printing in the report.
    For general plotting of the LFERs, use visualization.LFER
    
    Args:
        gas: An instance of the Gas class with Gas.LFER populated
        outliers: whether the plot should contain all outlier data
        show: whether the plot should be printed to the screen

    Returns:
        None
    '''

    slope_kd, int_kd, slope_b, int_b = gas.LFER.out
    pars = gas.LFER.pars
    pars_outliers = gas.LFER.pars_outliers

    log_kd0 = pars[:,0]
    deltaHd = pars[:,1]
    log_b0 = pars[:,2]
    deltaHb = pars[:,3]

    log_kd0_out = pars_outliers[:,0]
    deltaHd_out = pars_outliers[:,1]
    log_b0_out = pars_outliers[:,2]
    deltaHb_out = pars_outliers[:,3]

    fig, ax = plt.subplots(1,2, figsize=(6.2,3)) # never more than two plots

    cmap = cm.viridis
    norm = mcolors.Normalize(vmin=0, vmax=10)

    def lin_fit(data, slope, ints):
        return slope*data+ints

    if outliers:
        ax[0].plot(log_kd0_out, deltaHd_out,'o', color='gray', alpha=0.7)
    ax[0].plot(log_kd0, deltaHd,'o', color='black')
    ax[0].plot(log_kd0, lin_fit(log_kd0,slope_kd,int_kd), label='OLS', color=cmap(norm(1)))
    ax[0].set_xlabel(r'$\mathrm{ln}(k_{D,0})$')
    ax[0].set_ylabel(r'$\Delta H_D$')
    if outliers:
        ax[1].plot(log_b0_out, deltaHb_out,'o', color='gray', alpha=0.7)
    ax[1].plot(log_b0, deltaHb,'o', color = 'black')
    ax[1].plot(log_b0, lin_fit(log_b0,slope_b,int_b), label='OLS', color=cmap(norm(5)))
    ax[1].set_xlabel(r'$\mathrm{ln}(b_0)$')
    ax[1].set_ylabel(r'$\Delta H_b$')

    #for axi in ax:
        #axi.legend()
        #axi.set(adjustable='box', aspect='equal')

    plt.tight_layout()

    if show:
        plt.show()

def histograms(gas, show=False):
    '''Plots and saves the histograms from the van't Hoff fits

    Args:
        gas: An instance of the Gas class with Gas.vH populated
        show: whether the plot should be printed to the screen

    Returns:
        None 
    '''

    # *TO DO MAKE THIS GENERAL
    dms = gas.vH.dms
    dms_no_out = gas.vH.dms_no_out
    transposed_dms = np.transpose(dms)
    transposed_dms_no_out = np.transpose(dms_no_out)

    fig, ax = plt.subplots(2,3, figsize=(10,5))

    ax[0,0].set_xlabel(r'$\Delta H_{D,0}-\overline{\Delta H_{D,0}}$')
    ax[0,1].set_xlabel(r'$\Delta H_b-\overline{\Delta H_b}$')

    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=len(dms[0, :]))  # Normalize colors to the number of datasets

    for i in range(len(dms[0,:])):
        if i > 1:
            tex = "$C_H^{\\prime}-\\overline{C_H^{\\prime}}$"
            #ax[i//3,i%3].set_xlabel(f"$C_H^{{\\prime}}-[#{i-1}]-\\overline{{C}}$")    
            ax[i//3,i%3].set_xlabel(f"({tex}) #{i-1}")    

        color = cmap(norm(i))

        ax[i//3,i%3].set_ylabel('# Chains')

        # bins='auto' not used here because giving known memory error
        # precision not too impoprtant because binning is for visualization only
        ax[i//3,i%3].hist(transposed_dms[i]-np.mean(transposed_dms[i]),color='gray',bins='sqrt', alpha=0.5)
        #ax[i//3,i%3].hist(transposed_dms[i],color='black')
        ax[i//3,i%3].hist(transposed_dms_no_out[i]-np.mean(transposed_dms_no_out[i]),color=color, bins='sqrt')
    #print(transposed_dms)
    #print(transposed_dms_no_out)
    plt.tight_layout()

    if show:
        plt.plot()

def isotherms(gas, show=False):
    '''Plots and saves the sorption isotherms with DMS parameters
    
    Args:
        gas: An instance of the Gas class with Gas populated
        show: Whether the plot should be printed to the screen


    Returns:
        None
    '''

    # *TO DO: MAKE THIS GERNRAL
    fig, ax = plt.subplots(2,2, figsize=(7,7))

    c = gas.c
    cerr = gas.cerr
    p = gas.p
    T = gas.T
    LFE_params = gas.LFER.out
    avg_dms = gas.vH.avg_dms

    dHD_f = avg_dms[0]
    dHb_f = avg_dms[1]
    ch_vec_f = avg_dms[2:]
    a_kd0_f = LFE_params[0]
    b_kd0_f = LFE_params[1]
    a_b0_f = LFE_params[2]
    b_b0_f = LFE_params[3]
    kd0_f = np.exp((dHD_f-b_kd0_f)/a_kd0_f)
    b0_f = np.exp((dHb_f-b_b0_f)/a_b0_f)

    b_f = gas.b_f
    kd_f = gas.kd_f
    C_f = gas.vH.c_f
    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=len(T))  # Normalize colors to the number of datasets

    for i, temp in enumerate(T):

        press, c_model, c_model_err = evaluate.isotherm(gas, i)

        linecolor = cmap(norm(i))

        ax[i//2,i%2].set_ylabel(r'$C \; \mathrm{(cm^3_{STP} \; cm^{-3}_{pol})}$')
        ax[i//2,i%2].set_xlabel(r'$p_i \; \mathrm{or} \; f \; \mathrm{(atm)}$')

        ax[i//2,i%2].plot(press, c_model, color=linecolor, label='DMS_fit')
        ax[i//2,i%2].fill_between(press, c_model-c_model_err/2, c_model+c_model_err/2, color=linecolor, label='DMS_fit', alpha=0.3)

        ax[i//2,i%2].plot(p[i],c[i],'o', color='black')
        ax[i//2,i%2].errorbar(p[i],c[i],yerr=cerr[i],xerr=None,fmt='none', color='black')

    #for axi in ax: 
    #    axi.legend()

    plt.tight_layout()

    if show:
        plt.plot

def heat_of_sorption(gas, show=False):
    '''
    *
    '''

    temp = gas.T

    inv_RT = 1/(0.008314*temp)

    S_inf = gas.analysis.S_inf
    k_D = gas.kd_f
    b = gas.b_f

    S_inf_err = gas.analysis.S_inf_err
    k_D_err = gas.kd_err
    b_err = gas.b_err

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


    fig, ax = plt.subplots(1,3, figsize=(9,3))

    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=0, vmax=5)

    def lin_fit(data, slope, ints):
        return slope*data+ints

    #print(deltaH_S_inf)
    ax[0].plot(inv_RT, ln_S_inf, 'o', color='black')
    ax[0].plot(inv_RT, lin_fit(inv_RT,-deltaH_S_inf,ln_S_inf_0), label='OLS', color=cmap(norm(1)))
    ax[0].set_ylabel(r'$\mathrm{ln}(S_\infty)$')
    ax[0].set_xlabel(r'$(RT)^{-1}$')

    #     
    ax[1].plot(inv_RT, ln_k_D, 'o', color='black')
    ax[1].plot(inv_RT, lin_fit(inv_RT,-deltaH_D,ln_k_D_0), label='OLS', color=cmap(norm(2)))
    ax[1].set_ylabel(r'$\mathrm{ln}(k_D)$')
    ax[1].set_xlabel(r'$(RT)^{-1}$')

    ax[2].plot(inv_RT, ln_b, 'o', color='black')
    ax[2].plot(inv_RT, lin_fit(inv_RT,-deltaH_b,ln_b_0), label='OLS', color=cmap(norm(3)))
    ax[2].set_ylabel(r'$\mathrm{ln}(b)$')
    ax[2].set_xlabel(r'$(RT)^{-1}$')

    plt.tight_layout()

    if show:
        plt.show()
