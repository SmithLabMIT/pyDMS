Inputs and Outputs
==================
All parameters that can either be supplied to pyDMS or collected after an
optimization run are provided in :numref:`table-gas-class`.
``Gas`` here is the variable name of the ``Gas`` class.
For example: ``co2 = dms.Gas()``.

.. note::

   Attributes marked with :sup:`†` are **user-supplied inputs only** —
   pyDMS will never calculate or overwrite these values.

.. _table-gas-class:
.. list-table:: Attributes of the ``Gas`` class.
   :header-rows: 1
   :widths: 25 55 20

   * - Attribute
     - Explanation
     - Python Object
   * - **Gas**
     -
     -
   * - ``.formula`` :sup:`†`
     - Chemical formula
     - str
   * - ``.temp`` :sup:`†`
     - Temperatures
     - numpy.ndarray
   * - ``.p`` :sup:`†`
     - Pressures
     - numpy.ndarray
   * - ``.f``
     - Fugacities
     - numpy.ndarray
   * - ``.c`` :sup:`†`
     - Concentrations
     - numpy.ndarray
   * - ``.c_err`` :sup:`†`
     - Uncertainty in concentrations
     - numpy.ndarray
   * - ``.Z``
     - Compressibility factors
     - numpy.ndarray
   * - ``.kD``
     - :math:`k_{\mathrm{D},i}`
     - numpy.ndarray
   * - ``.b``
     - :math:`b_{i}`
     - numpy.ndarray
   * - ``.CH``
     - :math:`C_{\mathrm{H},i}'`
     - numpy.ndarray
   * - ``.kD_err``
     - :math:`\sigma_{k_{\mathrm{D},i}}`
     - numpy.ndarray
   * - ``.b_err``
     - :math:`\sigma_{b_{i}}`
     - numpy.ndarray
   * - ``.CH_err``
     - :math:`\sigma_{C_{\mathrm{H},i}'}`
     - numpy.ndarray
   * - ``.settings``
     - Settings (see :numref:`table-settings`)
     - dict
   * - ``.virial_coeff``
     - User-defined Virial coefficients
     - dict
   * - ``.pr_coeff``
     - User-defined Peng-Robinson parameters
     - dict
   * - **Gas.LFER**
     -
     -
   * - ``.fit``
     - :math:`\left[\alpha_\mathrm{D}, \beta_\mathrm{D}, \alpha_\mathrm{b}, \beta_\mathrm{b} \right]`
     - numpy.ndarray
   * - ``.fit_err``
     - :math:`\left[\sigma_{\alpha_\mathrm{D}}, \sigma_{\beta_\mathrm{D}}, \sigma_{\alpha_\mathrm{b}}, \sigma_{\beta_\mathrm{b}} \right]`
     - numpy.ndarray
   * - ``.out``
     - :math:`\left[ \ln(k_{\mathrm{D},0}),\Delta H_\mathrm{D},\ln(b_{0}),\Delta H_\mathrm{b} \right]`
     - numpy.ndarray
   * - ``.out_outliers``
     - :math:`\left[ \ln(k_{\mathrm{D},0}),\Delta H_\mathrm{D},\ln(b_{0}),\Delta H_\mathrm{b} \right]`
     - numpy.ndarray
   * - **Gas.vH**
     -
     -
   * - ``.out``
     - :math:`\left[ \Delta H_\mathrm{D}, \Delta H_\mathrm{b}, C_{\mathrm{H},i}' \right]`
     - numpy.ndarray
   * - ``.out_outliers``
     - :math:`\left[ \Delta H_\mathrm{D}, \Delta H_\mathrm{b}, C_{\mathrm{H},i}' \right]`
     - numpy.ndarray
   * - ``.residuals``
     - Residuals
     - numpy.ndarray
   * - ``.hessian_matrix``
     - Hessian matrix
     - numpy.ndarray
   * - **Gas.analysis**
     -
     -
   * - ``.S_inf``
     - :math:`S_\infty`
     - float
   * - ``.S_inf_err``
     - :math:`\sigma_{S_\infty}`
     - float
   * - ``.deltaH_S_inf``
     - :math:`\left[\Delta H_{\mathrm{S},\infty},S_{\infty,0}\right]`
     - numpy.ndarray
   * - ``.deltaH_D``
     - :math:`\left[\Delta H_{\mathrm{D}}, k_{\mathrm{D},0}\right]`
     - numpy.ndarray
   * - ``.deltaH_b``
     - :math:`\left[\Delta H_{\mathrm{b}}, b_0\right]`
     - numpy.ndarray
   * - ``.deltaH_S_inf_err``
     - :math:`\left[\sigma_{\Delta S_\infty},\sigma_{S_{\infty,0}}\right]`
     - numpy.ndarray
   * - ``.deltaH_D_err``
     - :math:`\left[\sigma_{\Delta H_\mathrm{D}},\sigma_{k_{\mathrm{D},0}}\right]`
     - numpy.ndarray
   * - ``.deltaH_b_err``
     - :math:`\left[\sigma_{\Delta H_\mathrm{b}},\sigma_{b_0}\right]`
     - numpy.ndarray
   * - ``.c_iso``
     - Concentrations for :math:`\Delta H_\mathrm{iso}`
     - numpy.ndarray
   * - ``.deltaH_iso``
     - :math:`\Delta H_\mathrm{iso}`
     - numpy.ndarray
   * - ``.deltaH_iso_err``
     - :math:`\sigma_{\Delta H_\mathrm{iso}}`
     - numpy.ndarray

.. note::

   Example usage: ``Gas.formula = 'CO2'``.

   ``.out`` only includes inliers from RANSAC regression.

   ``.out_outliers`` includes inliers and outliers.

:numref:`table-settings` details all the different settings that can be provided for the optimization. If a setting is not provided, the default value(s) will be used.

.. _table-settings:
.. list-table:: Optimization settings used by ``pyDMS`` within the ``Gas.settings`` attribute.
   :header-rows: 1
   :widths: 20 60 20

   * - Parameter
     - Explanation
     - Default
   * - ``dHD_guess``
     - Range between which to randomly choose initial guesses from
     - ``[-30, -1]``
   * - ``dHD_bounds``
     - Bounds optimized parameter must reside in
     - ``[-50, 0]``
   * - ``dHb_guess``
     - Range between which to randomly choose initial guess from
     - ``[-30,-1]``
   * - ``dHb_bounds``
     - Bounds optimized parameter must reside in
     - ``[-50, 0]``
   * - ``kD0_guess``
     - Range between which to randomly choose initial guess from
     - ``[0.001, 0.01]``
   * - ``kD0_bounds``
     - Bounds optimized parameter must reside in
     - ``[0, None]``
   * - ``b0_guess``
     - Range between which to randomly choose initial guess from
     - ``[0.001, 0.01]``
   * - ``b0_bounds``
     - Bounds optimized parameter must reside in
     - ``[0, None]``
   * - ``CH_guess``
     - Range between which to randomly choose initial guess from
     - ``[0, 100]``
   * - ``CH_bounds``
     - Bounds optimized parameter must reside in
     - ``[0, 150]``
   * - ``trials``
     - Number of trials to run each optimization procedure for
     - ``1000``
   * - ``solver_LFER``
     - Solver to use for LFER optimization
     - ``'SLSQP'``
   * - ``solver_vH``
     - Solver to use for van't Hoff optimization
     - ``'SLSQP'``
   * - ``maxiter_LFER``
     - Maximum iterations for LFER solver :cite:`noauthor_minimizemethodslsqp_nodate,noauthor_minimizemethodtrust-constr_nodate`
     - ``1000``
   * - ``maxiter_vH``
     - Maximum iterations for van't Hoff solver :cite:`noauthor_minimizemethodslsqp_nodate,noauthor_minimizemethodtrust-constr_nodate`
     - ``1000``
   * - ``ftol``
     - Termination tolerance on objective function :cite:`noauthor_minimizemethodslsqp_nodate`. Only for ``SLSQP``
     - ``1E-7``
   * - ``xtol``
     - Termination tolerance on parameter change :cite:`noauthor_minimizemethodtrust-constr_nodate`. Only for ``trust-constr``
     - ``1E-7``
   * - ``gtol``
     - Termination tolerance on gradient norm :cite:`noauthor_minimizemethodtrust-constr_nodate`. Only for ``trust-constr``
     - ``1E-7``
   * - ``verbose``
     - Print general information to the terminal
     - ``True``
   * - ``solver_verbose``
     - Print information from the solver to the terminal
     - ``False``
   * - ``seed``
     - Seed number to provide to the random number generator
     - ``None``

.. note::

   Example usage: ``Gas.settings = {'dHD_guess':[-40, -1], 'dHD_bounds':[-60, 0]}``.

   Available solver options are ``'SLSQP'`` :cite:`noauthor_minimizemethodslsqp_nodate` or ``'trust-constr'`` :cite:`noauthor_minimizemethodtrust-constr_nodate` from ``scipy.optimize.minimize`` :cite:`virtanen_scipy_2020`.

:numref:`table-fugacity` details the gases for which fugacity can automatically be computed from pressure data within pyDMS.

.. _table-fugacity:
.. list-table:: Built-in fugacity parameters and their corresponding references.
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Gas
     - Chemical Formula
     - ``Gas.formula``
     - Equation of State
     - Reference
   * - Helium
     - He
     - ``"He"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Hydrogen
     - H\ :sub:`2`
     - ``"H2"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Nitrogen
     - N\ :sub:`2`
     - ``"N2"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Oxygen
     - O\ :sub:`2`
     - ``"O2"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Methane
     - CH\ :sub:`4`
     - ``"CH4"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Carbon Dioxide
     - CO\ :sub:`2`
     - ``"CO2"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Ethane
     - C\ :sub:`2`\ H\ :sub:`6`
     - ``"C2H6"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Ethylene
     - C\ :sub:`2`\ H\ :sub:`4`
     - ``"C2H4"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Propane
     - C\ :sub:`3`\ H\ :sub:`8`
     - ``"C3H8"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Propylene
     - C\ :sub:`3`\ H\ :sub:`6`
     - ``"C3H6"``
     - Virial
     - :cite:`dymond_virial_2002`
   * - Hydrogen Sulfide
     - H\ :sub:`2`\ S
     - ``"H2S"``
     - Peng-Robinson
     - :cite:`linstrom_nist_2001`