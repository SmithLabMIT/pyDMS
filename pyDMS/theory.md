The overall goal of pyDMS is to perform non-linear regression on the dual mode sorption model

$$C=k_{D}f+\frac{C_{H}'b f}{1+b f}$$ 

to fit $k_D$, $C_H'$, and $b$ to an experimental sorption isotherm. The default loss function used is the $\chi^2$-function defined as

$$\chi^2=\sum{\left(\frac{C-O}{\sigma_O}\right)^2}$$

To do so, the optimization procedure applies linear free energy relationship (LFER) constraints followed by van't Hoff constraints to obtain reproducible results.

The mathematical complixities in doing so are explained below.

Linear Free Energy Relationships

The 

The LFER code solves the optimization problem

$$\Delta H_D=\alpha_D\ln(k_{D,0})+\beta_D$$
$$\Delta H_b=\alpha_b\ln(b_{0})+\beta_b$$
$$C_i=k_{D,i}f_i+\frac{C_{H,i}'b_if_i}{1+b_if_i}$$

Where k_(D,0), k_(b,0), ΔH_D, ΔH_b, and C_i are fitting parameters.

$$k_D=k_{D,0}\exp\left({-\frac{1000\Delta H_D}{8.314T}}\right)$$

$$b=b_{0}\exp\left({-\frac{1000\Delta H_b}{8.314T}}\right)$$

$$k_D = k_{D,0}=\exp\left({\frac{\Delta H_D - \beta_{k_{D,0}}}{\alpha_{k_{D,0}}}}\right)$$

$$b = b_{0}=\exp\left({\frac{\Delta H_b - \beta_{b+{0}}}{\alpha_{b_{0}}}}\right)$$

$$-\frac{\mathrm{median(|A-median(A)|)}}{\sqrt{2}\;\mathrm{erfcinv}(3/2)}$$

$$\frac{\partial C}{\partial H_D} = \exp\left(-\frac{120.279\,\Delta H_D}{T} + \frac{-b_{k_{D,0}} + \Delta H_D}{\alpha_{k_{D,0}}}\right)\;p\;\left((-\frac{120.279}{T} + \frac{1}{a_{kd0}}\right))$$

$$\frac{\partial C}{\partial H_b} = 
\frac{c_h\;\exp\!\Bigl(1 + \tfrac{120.279\,dHb}{T} + \tfrac{b_{b0} + dHb}{a_{b0}}\Bigr)\;p\;(T - 120.279\,a_{b0})}
{\Bigl(\exp\!\bigl(\tfrac{b_{b0}}{a_{b0}} + \tfrac{120.279\,dHb}{T}\bigr) + \exp\!\bigl(\tfrac{dHb}{a_{b0}}\bigr)\;p\Bigr)^2\;T\;a_{b0}}$$

$$\frac{\partial C}{\partial C_H'} = \Biggl(1 + \frac{\exp\!\bigl(\tfrac{b_{b0} - dHb}{a_{b0}} + \tfrac{120.279\,dHb}{T}\bigr)}{p}\Biggr)^{-1}$$