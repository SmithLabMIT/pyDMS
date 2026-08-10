Mathematical Formulation
========================

Dual-Mode Sorption Minimization
-------------------------------

The overarching goal of pyDMS is fit the dual-mode sorption (DMS) model

.. math::
      :label: eq-dms

      C=k_\mathrm{D}f+\frac{C_\mathrm{H}'bf}{1+bf}

to experimental data. Thus, pyDMS seeks values to the parameters :math:`k_\mathrm{D}`, :math:`C_\mathrm{H}'`, and :math:`b` that minimize the :math:`\chi^2` loss function:

.. math::
       :label: eq-chi2

       \chi^2= \sum \left(\frac{C_i-C_{i,exp}}{\sigma_{C_i,exp}}\right)^2

To do so, the optimization procedure applies linear free energy relationship (LFER) constraints, :cite:`freeman_basis_1999` followed by van't Hoff constraints :cite:`wu_elucidating_2021, barrer_solution_1962, koros_energetics_1979, smith_hydrogen_2013` with the aim of obtaining more reproducible results than simple nonlinear regression. Below, we work through the mathematical complexities in doing so.

pyDMS performs a sequential double-optimization procedure to reduce the number of independent variables being fit at a given time. First, the LFER section fits the DMS model :eq:`eq-dms` by minimizing :math:`\chi^2` :eq:`eq-chi2` subject to the LFER constraints,

.. math::
       :label: eq-LFER_HD

       \Delta H_\mathrm{D}=\alpha_\mathrm{D}\ln(k_{\mathrm{D},0})+\beta_\mathrm{D}

and

.. math::
       :label: eq-LFER_Hb

       \Delta H_\mathrm{b}=\alpha_\mathrm{b}\ln(b_0)+\beta_\mathrm{b},

where the fitting parameters are :math:`k_{\mathrm{D},0}`, :math:`k_{\mathrm{b},0}`, :math:`\Delta H_\mathrm{D}`, :math:`\Delta H_\mathrm{b}`, and :math:`C_{\mathrm{H},i}`. Given that :math:`k_\mathrm{D}` and :math:`b` follow van't Hoff relationships, they can be solved as

.. math::
       :label: eq-kD_vh

       k_\mathrm{D}=k_{\mathrm{D},0} \exp\left( -\frac{\Delta H_{\mathrm{D}}}{RT} \right)

and

.. math::
       :label: eq-b_vh

       b=b_0\exp\left( -\frac{\Delta H_\mathrm{b}}{RT} \right)

The optimization is performed ``trials`` number of times, each using a different uniformly random inital guesses bounded by the parameters in :numref:`table-settings`. Solver options are provided in :numref:`table-gas-class`. The resulting :math:`\alpha_i` and :math:`\beta_i` (:eq:`eq-LFER_HD` and :eq:`eq-LFER_Hb`) are then extracted via linear regression using the RANSAC algorithm :cite:`fischler_random_1981` to remove outliers.

The benefit to performing this optimization becomes clear when we perform the van't Hoff fitting. Rearranging :eq:`eq-LFER_HD` and :eq:`eq-LFER_Hb` allows for :math:`k_{\mathrm{D},0}` and :math:`b_0` to be solved for, removing the necessity of reoptimization. Now, the van't Hoff expressions (:eq:`eq-kD_vh` and :eq:`eq-b_vh`) can be constrained with only one free variable each (:math:`\Delta H_\mathrm{D}` and :math:`\Delta H_\mathrm{b}`, respectively).

Again, the optimization is performed as described above. Outliers are removed using a scaled median absolute deviation technique. :cite:`noauthor_matlab_2024` That is, a datapoint is classified as an outlier if the following criteria is met:

.. math::

       x>\frac{-\mathrm{median}[\mathrm{abs}(\mathbf{A}-\mathrm{median}[\mathbf{A}])]}{\sqrt{2}\;\mathrm{erfcinv}(3/2)}

where :math:`x` is a datapoint in the array, :math:`\mathbf{A}`, and erfcinv is the inverse complementary error function. The mean of the resulting inliers is calculated to yield the final DMS parameters, :math:`C_{\mathrm{H},i}`, :math:`k_\mathrm{D}`, and :math:`b`.

Error Propagation
~~~~~~~~~~~~~~~~~

To gather uncertainty of the DMS parameters, error within the optimized parameters is determined directly from the Hessian matrix of the loss function, :math:`\mathbf{\mathrm{H}(x)}` where :math:`\mathbf{x}` is the array of optimized parameters:

.. math::

       \mathbf{x}=\left[\Delta H_\mathrm{D}, \Delta H_\mathrm{b}, C_{\mathrm{H},i}' \right]

The Hessian is defined as

.. math::
       :label: eq-hessian

       \textbf{H}(\textbf{x}) = 
       \begin{bmatrix}
       \frac{\partial^2 \chi^2}{\partial^2 \Delta H_\mathrm{D}} & \frac{\partial^2 \chi^2}{\partial \Delta H_\mathrm{D} \partial \Delta H_\mathrm{b}} & \cdots & \frac{\partial^2 \chi^2}{\partial \Delta H_\mathrm{D} \partial C_{\mathrm{H},j}'} & \cdots \\
       \frac{\partial^2 \chi^2}{\partial \Delta H_\mathrm{b} \partial \Delta H_\mathrm{D}} & \frac{\partial^2 \chi^2}{\partial^2 \Delta H_\mathrm{b}} & \cdots & \frac{\partial^2 \chi^2}{\partial \Delta H_\mathrm{b} \partial C_{\mathrm{H},j}'} & \cdots \\
       \vdots & \vdots & \ddots & \vdots & \\
       \frac{\partial^2 \chi^2}{\partial C_{\mathrm{H},k}' \partial \Delta H_\mathrm{D}} & \frac{\partial^2 \chi^2}{\partial C_{\mathrm{H},k}' \partial \Delta H_\mathrm{b}} & \cdots & \frac{\partial^2 \chi^2}{\partial C_{\mathrm{H},k}' \partial C_{\mathrm{H},j}'} & \cdots \\
       \vdots & \vdots & & \vdots & \ddots
       \end{bmatrix}

Since there is no functional form fit to :math:`C_{\mathrm{H}, i}`, there is assumed to be negligible dependence between their values, thus,

.. math::

       \frac{\partial^2 \chi^2}{\partial C_{\mathrm{H},k}' \partial C_{\mathrm{H},j}'}=0

Taking the second derivative of :math:`\chi^2`, defined in :eq:`eq-chi2`, for arbitrary variables :math:`x_j` and :math:`x_k` results in

.. math::
       :label: eq-d2_chi2

       \frac{\partial^2 \chi^2}{\partial x_j \, \partial x_k} = 2\sum  \left[\left( \frac{1}{\sigma_{C_{i,\mathrm{exp}}}^2} \frac{\partial C_i}{\partial x_j} \frac{\partial C_i}{\partial x_k} \right) + \left( \frac{C_i - C_{i,\mathrm{exp}}}{\sigma^2_{C_{i,\mathrm{exp}}}} \right)\frac{\partial^2 C_i}{\partial x_j \, \partial x_k} \right]

Assuming that the model is able to describe the data well and thus the residuals are negligible, :eq:`eq-d2_chi2` is simplified to

.. math::
       :label: eq-d2_chi2-simplified

       \frac{\partial^2 \chi^2}{\partial x_j \, \partial x_k} \approx \sum \frac{2}{\sigma_{C_{i,\mathrm{exp}}}^2} \frac{\partial C_i}{\partial x_j} \frac{\partial C_i}{\partial x_k}

The concentration to differentiate is from the DMS model :eq:`eq-dms`, with van't Hoff and LFER Equations (:eq:`eq-LFER_HD`--:eq:`eq-b_vh`) substituted in to yield

.. math::
       :label: eq-dms-full-sub

       C=\exp\left(-\frac{\Delta H_\mathrm{D}}{RT} + \frac{-\beta_\mathrm{D} + \Delta H_\mathrm{D}}{\alpha_\mathrm{D}} \right) f +
       \frac{C_\mathrm{H}'\exp\left( -\frac{\Delta H_\mathrm{b}}{RT} + \frac{-\beta_\mathrm{b} + \Delta H_\mathrm{b}}{\alpha_\mathrm{b}} \right) f
       }{ 1 + \exp\left( -\frac{\Delta H_\mathrm{b}}{RT} + \frac{-\beta_\mathrm{b} + \Delta H_\mathrm{b}}{\alpha_\mathrm{b}} \right) f}

Taking the resulting derivatives of :eq:`eq-dms-full-sub` that appear in the Hessian yield:

.. math::
       :label: eq-dCdHD

       \frac{\partial C}{\partial H_\mathrm{D}}=\exp\left(-\frac{\Delta H_\mathrm{D}}{RT} + \frac{-\beta_{\mathrm{D}}+\Delta H_\mathrm{D}}{\alpha_{\mathrm{D}}}\right)\left(-\frac{1}{RT} + \frac{1}{\alpha_{\mathrm{D}}}\right)f,

.. math::
       :label: eq-dCdHb

       \frac{\partial C}{\partial H_\mathrm{b}} = \frac{C_\mathrm{H}'\exp\left( \frac{\Delta H_{\mathrm{b}}}{RT} + \frac{\beta_{\mathrm{b}} + \Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right)(RT - \alpha_{\mathrm{b}})f}{\left( \exp\left( \frac{\beta_{\mathrm{b}}}{\alpha_{\mathrm{b}}} + \frac{\Delta H_{\mathrm{b}}}{RT} \right) + \exp\left( \frac{\Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right) f \right)^2 RT \alpha_{\mathrm{b}}},

and

.. math::
       :label: eq-dCdCH

       \frac{\partial C}{\partial C_\mathrm{H}'} = \left(1+\frac{\exp\left(\tfrac{\beta_{\mathrm{b}}-\Delta H_\mathrm{b}}{\alpha_\mathrm{b}}+\tfrac{\Delta H_\mathrm{b}}{RT}\right)}{f}\right)^{-1}

The Hessian used to determine the standard deviations of optimized parameters is the mean Hessian of successful van't Hoff optimizations, :math:`\overline{\mathrm{\bf{H}}}\bf(x)`. The covariance matrix is computed as:

.. math::

       \mathrm{Cov(\textbf{x})}=\frac{s^2}{2}\overline{\bf{H}}(x)^{-1}

where :math:`s` is equal to 1 as the experimental uncertanties are known and the inverse Hessian, :math:`\overline{\mathbf{H}}\mathbf(x)^{-1}`, is computed via Cholesky decomposition. :cite:`cholesky_note_1924, higham_analysis_1990, bevington_data_2003` From the covariance matrix, the standard deviations of the parameters in :math:`\mathbf{x}` are found as the diagonal entries.

The uncertainties in :math:`\mathbf{x}` are propagated to find the uncertainties in the DMS parameters. Error in :math:`k_\mathrm{D}` is found as:

.. math::

       \sigma_{k_\mathrm{D}}=\sqrt{
       \sigma_T^2 \left( \frac{\partial k_\mathrm{D}}{\partial T} \right)^2 + \sigma_{\alpha_D}^2 \left( \frac{\partial k_\mathrm{D}}{\partial \alpha_\mathrm{D}} \right)^2 +\sigma_{\beta_\mathrm{D}}^2 \left( \frac{\partial k_{\mathrm{D}}}{\partial \beta_\mathrm{D}} \right)^2 +\sigma_{\Delta H_{\mathrm{D}}}^2 \left( \frac{\partial k_{\mathrm{D}}}{\partial \Delta H_{\mathrm{D}}} \right)^2}

where error in the temperature, :math:`\sigma_T^2`, is assumed to be constant at :math:`0.01\;\mathrm{K}` and the partial derivatives are equal to

.. math::
       :label: eq-dkDdT

       \frac{\partial k_\mathrm{D}}{\partial T}=\frac{\exp\left( -\frac{\Delta H_{\mathrm{D}}}{RT} + \frac{-\beta_{\mathrm{D}} + \Delta H_{\mathrm{D}}}{\alpha_{\mathrm{D}}} \right) \Delta H_{\mathrm{D}}}{RT^2},

.. math::
       :label: eq-dkdalphaD

       \frac{\partial k_\mathrm{D}}{\partial \alpha_\mathrm{D}}=\frac{\exp\left( -\frac{\Delta H_{\mathrm{D}}}{RT} + \frac{-\beta_{\mathrm{D}} + \Delta H_{\mathrm{D}}}{\alpha_{\mathrm{D}}} \right) \left(\beta_{\mathrm{D}} - \Delta H_{\mathrm{D}} \right)}{\alpha_{\mathrm{D}}^2},

.. math::
       :label: eq-dkDbetaD

        \frac{\partial k_\mathrm{D}}{\partial \beta_\mathrm{D}}=- \frac{\exp\left( -\frac{\Delta H_{\mathrm{D}}}{RT} + \frac{-\beta_{\mathrm{D}} + \Delta H_{\mathrm{D}}}{\alpha_{\mathrm{D}}} \right)}{\alpha_{\mathrm{D}}},

and

.. math::
       :label: eq-dkDdHD

       \frac{\partial k_\mathrm{D}}{\partial \Delta H_\mathrm{D}}=\exp\left( -\frac{\Delta H_{\mathrm{D}}}{RT} + \frac{-\beta_{\mathrm{D}} + \Delta H_{\mathrm{D}}}{\alpha_{\mathrm{D}}} \right)\left( -\frac{1}{RT} + \frac{1}{\alpha_{\mathrm{D}}} \right)

Likewise, a similar approach can be followed for the uncertainty in :math:`b`:

.. math::

       \sigma_{b}=\sqrt{
       \sigma_T^2 \left( \frac{\partial b}{\partial T} \right)^2 + \sigma_{\alpha_b}^2 \left( \frac{\partial b}{\partial \alpha_\mathrm{b}} \right)^2 +\sigma_{\beta_\mathrm{b}}^2 \left( \frac{\partial b}{\partial \beta_\mathrm{b}} \right)^2 +\sigma_{\Delta H_{\mathrm{b}}}^2 \left( \frac{\partial b}{\partial \Delta H_{\mathrm{b}}} \right)^2}

where error in the temperature, :math:`\sigma_T^2`, is assumed to be constant at :math:`0.01\;\mathrm{K}` and the partial derivatives are equal to

.. math::
       :label: eq-dbdT

       \frac{\partial b}{\partial T}=\frac{\exp\left( -\frac{\Delta H_{\mathrm{b}}}{RT} + \frac{-\beta_{\mathrm{b}} + \Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right) \Delta H_{\mathrm{b}}}{RT^2},

.. math::
       :label: eq-dbdalphaD

       \frac{\partial b}{\partial \alpha_\mathrm{b}}=\frac{\exp\left( -\frac{\Delta H_{\mathrm{b}}}{RT} + \frac{-\beta_{\mathrm{b}} + \Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right) \left(\beta_{\mathrm{b}} - \Delta H_{\mathrm{b}} \right)}{\alpha_{\mathrm{b}}^2},

.. math::
       :label: eq-dbbetaD

       \frac{\partial b}{\partial \beta_\mathrm{b}}=- \frac{\exp\left( -\frac{\Delta H_{\mathrm{b}}}{RT} + \frac{-\beta_{\mathrm{b}} + \Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right)}{\alpha_{\mathrm{b}}},     

and

.. math::
       :label: eq-dbdHD

       \frac{\partial b}{\partial \Delta H_\mathrm{b}}=\exp\left( -\frac{\Delta H_{\mathrm{b}}}{RT} + \frac{-\beta_{\mathrm{b}} + \Delta H_{\mathrm{b}}}{\alpha_{\mathrm{b}}} \right)\left( -\frac{1}{RT} + \frac{1}{\alpha_{\mathrm{b}}} \right),

From these errors within the DMS parameters, the error within the concentration is calculated as

.. math::

       \sigma_C'=\sqrt{\frac{C_\mathrm{H}'^2 f^2 \sigma_b^2}{\left( 1+bf\right)^4}+\frac{b^2f^2\sigma_{C_\mathrm{H}'}^2}{\left( 1+bf \right)^2}+f^2\sigma_{k_\mathrm{D}}^2}

Sorption Coefficient at Infinite Dilution
-----------------------------------------

The sorption coefficient at infinite dilution, :math:`S_\infty`, is calculated from the DMS parameters as

.. math::

       S_\infty=k_\mathrm{D}+C_\mathrm{H}'b

.. _error-propagation-1:

Error Propagation
~~~~~~~~~~~~~~~~~

Uncertainty in :math:`S_\infty` is calculated as

.. math::

   \sqrt{C_\mathrm{H}'^2\sigma_b^2+b^2\sigma_{C_\mathrm{H}'}^2+\sigma_{k_\mathrm{D}}^2}

Energetics
----------

Within the optimization, the heats (enthalpies) of Henry and Langmuir sorption are from :eq:`eq-kD_vh` and :eq:`eq-b_vh`. A third heat of sorption, the infinite dilution heat of sorption (:math:`\Delta H_{\mathrm{S},\infty}`), also defined as a van't Hoff fit:

.. math::
       :label: eq-vH_Sinf

       \ln{S_\infty}=\ln{S_{\infty,0}}-\frac{\Delta H_{\mathrm{S},\infty}}{RT}

is computed via an ordinary least squares (OLS) regression of :eq:`eq-vH_Sinf` from the parameters within the optimization.

.. _error-propagation-2:

Error Propagation
~~~~~~~~~~~~~~~~~

Uncertainties within the Henry and Langmuir heats of sorption are directly calculated from the inverse Hessian. Uncertainty within the infinite dilution heat of sorption is taken as the uncertainty within the slope of the OLS regression, where the independent variable is :math:`(RT)^{-1}`.

Isosteric Heat of Sorption
~~~~~~~~~~~~~~~~~~~~~~~~~~

The isosteric heat of sorption, :math:`\Delta H_\mathrm{iso}`, is calculated via the Clausius--Clapeyron approach :cite:`koros_energetics_1979`:

.. math::

   \left(\frac{\mathrm{d}f}{\mathrm{d}(1/T)}\right)_C=\frac{\Delta H_\mathrm{iso}}{ZR}

where :math:`f` would be substituted for :math:`P` if pressure-based DMS parameters are used and :math:`Z` is the compressibility factor (see :numref:`sec-fugacity`). :math:`\Delta H_\mathrm{iso}` is computed via OLS regression of :math:`\ln(f)` versus :math:`T^{-1}`. Fugacities at each temperature are computed from a constant concentration value calculated from the DMS parameters optimized within pyDMS. The compressibility factor at each concentration is found as the mean of the compressibility factors at each temperature.

.. _error-propagation-3:

Error Propagation
~~~~~~~~~~~~~~~~~

Uncertainty in :math:`\Delta H_\mathrm{iso}` is calculated from the uncertainty within the slope of the OLS regression as

.. math::

   \sigma_{\Delta H_\mathrm{iso}}=\sigma_m\cdot Z \cdot R

Mixed-Gas DMS and Sorption Selectivities
----------------------------------------

The mixed-gas DMS model extends the pure-gas DMS model to an :math:`N`-component mixture

.. math::
       :label: eq-mixed-dms

       C_i=k_{\mathrm{D},i}f_i+\frac{C_{\mathrm{H},i}'b_if_i}{1+b_if_i+\cdots+b_Nf_N}

where :math:`i` is the gas whose sorption isotherm is of interest. From these isotherms, binary sorption selectivity is calculated as

.. math::
       :label: eq-binary-selectivity

       \alpha_{i/j}^s=\frac{C_i/f_i}{C_j/f_j}

where :math:`f_i` and :math:`f_j` are defined by the Lewis and Randall rule :cite:`tester_thermodynamics_1997` such that

.. math::
       :label: eq-lewis-randall

       f_k=f_\mathrm{total}\cdot y_k

where :math:`f_\mathrm{total}` is the fugacity of the gas of interest at the total pressure and :math:`y_k` is the mole fraction of the component.

Calculations of sorption selectivity can be extended to any number of gases. Currently, pyDMS also has ternary selectivities built in, defined as

.. math::
       :label: eq-mixed-selec-2up

       \alpha_{(i+j)/k}^s=\frac{C_{i}/f_{i}+C_{j}/f_{j}}{C_{k}/f_{k}}

and

.. math::
       :label: eq-mixed-selec-1up

       \alpha_{i/(j+k)}^s=\frac{C_{i}/f_{i}}{C_{j}/f_{j}+C_{k}/f_{k}}

.. _error-propagation-4:

Error Propagation
~~~~~~~~~~~~~~~~~

Error propagation on the mixed-gas DMS model :eq:`eq-mixed-dms` yields the following uncertainty within concentration

.. math::
       :label: eq-mixed-dms-uncertainty

       \sigma_{C_i}=\sqrt{
       \frac{
              C_{\mathrm{H}_i}^{'2} f_i^2 \left( 1 + \sum_{j\neq i}^{N} b_{j} f_{j} \right)^2 \sigma_{b_i}^2
       }{
              \left( 1 + \sum_{j}^{N} b_{j} f_{j} \right)^4
       }
       +
       \frac{C_{\mathrm{H},i}'^2 b_i^2 f_i^2}{\left( 1 + \sum_{j}^{N} b_{j} f_{j} \right)^4} \sum_{j \neq i}^{N} f_j^2 \sigma_{b_j}^2
       +
       \frac{
              b_i^2 f_i^2 \sigma_{C_{\mathrm{H},i}}^2}{\left( 1 + \sum_{j}^{N} b_{j} f_{j} \right)^2}
       +
       f_i^2 \sigma_{k_{\mathrm{D},i}}^2
       }

Binary sorption selectivity :eq:`eq-binary-selectivity` is error propagated as

.. math::

       \sigma_{\alpha_{i/j}}=\sqrt{\frac{f_j^2 \left( C_j^2 \sigma_{C_i}^2 + C_i^2 \sigma_{C_j}^2 \right)}{C_j^4 f_i^2}}

Ternary sorption selectivities (:eq:`eq-mixed-selec-2up` and :eq:`eq-mixed-selec-1up`) are error propagated as

.. math::

       \sigma_{\alpha_{(i+j)/k}}=\sqrt{
   \frac{f_k^2 \, \sigma_{c_i}^2}{c_k^2 f_i^2} +
   \frac{f_k^2 \, \sigma_{c_j}^2}{c_k^2 f_j^2} +
   \frac{\left( \frac{c_i}{f_i} + \frac{c_j}{f_j} \right)^2 f_k^2 \, \sigma_{c_k}^2}{c_k^4}
   }\textbf{}

and

.. math::

       \sigma_{\alpha_{i/(j+k)}}=\sqrt{
   \frac{f_j^2 f_k^2 \, \sigma_{c_i}^2}{\left( c_k f_i f_j + c_j f_i f_k \right)^2}
   +
   \frac{c_i^2 f_j^2 f_k^4 \, \sigma_{c_j}^2}{f_i^2 \left( c_k f_j + c_j f_k \right)^4}
   +
   \frac{c_i^2 f_j^4 f_k^2 \, \sigma_{c_k}^2}{f_i^2 \left( c_k f_j + c_j f_k \right)^4}
   },

respectively.

.. _sec-fugacity:

Fugacity
--------

Deviation from ideal gas behavior can be represented using the compressibility factor

.. math::

       Z=\frac{p_iV_m}{RT}

where :math:`p_i` is the partial pressure, :math:`V_m` is the molar volume, :math:`T` is the absolute temperature, and :math:`R` is the ideal gas constant. The compressibility factor is related to the :math:`f` and :math:`p_i` as

.. math::

       \ln(\phi)=\ln\left(\frac{f}{p_i}\right)=\int_{0}^{p_i}\frac{Z-1}{p_i}dp_i

where :math:`\phi` is defined as the fugacity coefficient. :cite:`tester_thermodynamics_1997`

Virial Expansion
~~~~~~~~~~~~~~~~

The compressibility factor can be expressed in terms of the Virial coefficients as

.. math::

       Z=1+B^*p+C^*p^2+\cdots

where :math:`B^*` is related to the second Virial coefficient, :math:`B`, as

.. math::

       B=B^*RT

and :math:`C^*` is related to the third Virial, C, coefficient as

.. math::

       C=(B^{*2}+C^*)(RT)

In pyDMS, the Virial coefficients are expressed as temperature-dependent expansions:

.. math::

       B=B_0+B_1T^{-1}+B_2T^{-2}+B_3T^{-3}+B_4T^{-4}

and

.. math::
       :label: eq-virial-C

       C=C_0+C_1T^{-1}+C_2T^{-2}+C_3T^{-3}+C_4T^{-4}

where :math:`B` has units of :math:`\mathrm{cm^3\; mol^{-1}}`, :math:`C` has units of :math:`\mathrm{cm^6\; mol^{-2}}`, and :math:`T` is the temperature in Kelvin. Data for the Virial coefficients are taken from Dymond et al. :cite:`dymond_virial_2002` The functional form was provided for the second Virial coefficient within Dymond et al., :cite:`dymond_virial_2002` however, no functional form was fit for :math:`C`, therefore, we calculated the fit using provided experimental data at temperatures greater than 200 K and shown in `[tab:third-virial-coeffs] <#tab:third-virial-coeffs>`__.

.. list-table:: Second virial coefficients :math:`C_j` for supported gases.
   :header-rows: 1
   :widths: 20 16 16 16 16 16

   * - Gas
     - :math:`C_0`
     - :math:`C_1`
     - :math:`C_2`
     - :math:`C_3`
     - :math:`C_4`
   * - He
     - 0.10023
     - 0
     - 0
     - 0
     - 0
   * - H\ :sub:`2`
     - 0.53597
     - 0
     - 0
     - 0
     - 0
   * - N\ :sub:`2`
     - 0.45178
     - 282.49
     - 0
     - 0
     - 0
   * - O\ :sub:`2`
     - 0.91432
     - -57.00
     - 38999
     - 0
     - 0
   * - CH\ :sub:`4`
     - 1.46800
     - -417.68
     - 211330
     - 0
     - 0
   * - CO\ :sub:`2`
     - 8.22730
     - -11176
     - 5297100
     - -673400
     - 0
   * - C\ :sub:`2`\ H\ :sub:`6`
     - -21.966
     - 19216
     - -2910000
     - 0
     - 0
   * - C\ :sub:`2`\ H\ :sub:`4`
     - -19.585
     - 14199
     - -1879000
     - 0
     - 0
   * - C\ :sub:`3`\ H\ :sub:`8`
     - 161.600
     - -211730
     - 95225000
     - -13420000
     - 0
   * - C\ :sub:`3`\ H\ :sub:`6`
     - -11.713
     - 9511.10
     - 0
     - 0
     - 0

.. note::

   All coefficients are in units of :math:`\left( \mathrm{cm^6\; mol^{-2} \; K^{-j}\cdot10^{3}} \right)` where :math:`j` corresponds to :math:`C_j` to match the style of Dymond et al. :cite:`dymond_virial_2002` (i.e., each number must be multiplied by :math:`10^{-3}` before using :eq:`eq-virial-C`).

The fugacity coefficient can be expressed in terms of the Virial coefficients as

.. math::

   \ln(\phi)=B^*p_i+\frac{1}{2}C^*p_i^2+\cdots

Peng--Robinson Equation of State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fugacity coefficient using the Peng--Robinson :cite:`peng_new_1976` EoS is solved as

.. math::

       \ln(\phi)=(Z-1)-\ln(Z-B)-\frac{A}{2\sqrt{2}B}\ln \left( \frac{Z+(\sqrt{2}+1)B}{Z-(\sqrt{2}-1)B} \right)

where :math:`Z` is the compressibility factor which can be solved via the cubic equation

.. math::
       :label: eq-Z-soln
       
       Z^3+(B-1)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0

with the constants :math:`A` and :math:`B` are defined as

.. math::

       A=a\frac{P}{R^2T^2}

and

.. math::

       B=b\frac{P}{RT}

The parameters :math:`a` and :math:`b` are the Peng-Robinson coefficients, defined as

.. math::

       a=0.457 \left( 1+\kappa \left( 1-\sqrt{\frac{T}{T_c}} \right) \right)^2\frac{R^2T_c^2}{P_c}

and

.. math::

       b=0.0778\frac{R T_c}{P_c}

where

.. math::

       \kappa=0.375+1.542\omega-0.270\omega^2,

and :math:`T_c`, :math:`P_c`, and :math:`\omega` are the critical temperature, critical pressure, and acentric factor, respectively. The largest root in `eq-Z-soln`` is always selected as it corresponds to the vapor phase. :cite:`tester_thermodynamics_1997`
