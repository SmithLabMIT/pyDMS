import numpy as np
import matplotlib.pyplot as plt
import pyDMS.dms as dms
import pyDMS.evaluate as ev
import pyDMS.multigas as mg

# assuming compute_dms() run previously for fugacity-based DMS
co2 = dms.load_gas_class('co2.pkl')
h2s = dms.load_gas_class('h2s.pkl')
ch4 = dms.load_gas_class('ch4.pkl')

# computing the pure isotherm for co2
f_co2_pure, c_co2_pure, cerr_co2_pure = ev.isotherm(co2, temp=308)

# specifying the total fugacity (f_CO2+f_H2S+f_CH4) of the mixture
f_tot = np.linspace(4, 20, 1000)

# computing mixed-gas isotherm of CO2
f_co2, c_co2, cerr_co2 = mg.mixed_isotherm(
                            co2, h2s, ch4,  # gas of interest always first
                            p_or_f=f_tot,  # total fugacity
                            mol_frac=[0.2, 0.2, 0.6],  # in gas order
                            temp=308)  # K

# computing mixed-gas isotherm of CH4
f_ch4, c_ch4, cerr_ch4 = mg.mixed_isotherm(
                            ch4, h2s, co2,  # gas of interest always first
                            p_or_f=f_tot,  # total fugacity
                            mol_frac=[0.6, 0.2, 0.2],  # in gas order
                            temp=308)  # K

# computing mixed-gas isotherm of H2S
f_h2s, c_h2s, cerr_h2s = mg.mixed_isotherm(
                            h2s, ch4, co2,  # gas of interest always first
                            p_or_f=f_tot,  # total fugacity
                            mol_frac=[0.2, 0.6, 0.2],  # in gas order
                            temp=308)  # K

# computing mixed-gas sorption selectivity of (CO2+H2S)/CH4
f_tot, selec, selecerr = mg.selectivity(
                        [f_co2, c_co2, cerr_co2], # numerator
                        [f_h2s, c_h2s, cerr_h2s], # numerator
                        [f_ch4, c_ch4, cerr_ch4], # denominator
                        calc='2_numerator')  # both acid gases on top

# if calc=1_numerator, the code would perform CO2/(H2S+CH4)
# for binary selectivity, calc does not have to be specified

# visualizing the selectivity, for example
plt.plot(f_tot, selec, label=f"T={co2.temp[0]} K")
plt.fill_between(f_tot, selec-selecerr, selec+selecerr, alpha=0.5)
plt.xlabel("f_CO2+f_H2S+f_CH4 (atm)")
plt.ylabel("Sorption selectivity (CO2+H2S)/CH4")
plt.legend()
plt.show()

# if you want to save the data, for example
np.savetxt("selec.txt", selec)
