import numpy as np
import pyDMS.dms as dms
import pyDMS.evaluate as ev
import pyDMS.multigas as mg

co2 = dms.load_gas_class('co2.pkl')  # assuming compute_dms() run previously
h2s = dms.load_gas_class('h2s.pkl')  # assuming compute_dms() run previously
ch4 = dms.load_gas_class('ch4.pkl')  # assuming compute_dms() run previously

# computing the pure isotherm of co2
p_co2_pure, c_co2_pure, cerr_co2_pure = ev.isotherm(co2, index=0)

# specifying the total pressure of the mixture
p_tot = np.linspace(0.01, 40, 1000)

# computing mixed-gas isotherm of CO2
p_co2, c_co2, cerr_co2 = mg.mixed_isotherm(
                            co2, h2s, ch4,  # gas of interest always first
                            p_or_f=p_tot,  # total pressure
                            mol_frac=[0.2, 0.2, 0.6],  # in gas order
                            temp=308)  # K

# computing mixed-gas isotherm of CH4
p_ch4, c_ch4, cerr_ch4 = mg.mixed_isotherm(
                            ch4, h2s, co2,  # gas of interest always first
                            p_or_f=p_tot,  # total pressure
                            mol_frac=[0.6, 0.2, 0.2],  # in gas order
                            temp=308)  # K

# computing mixed-gas isotherm of H2S
p_h2s, c_h2s, cerr_h2s = mg.mixed_isotherm(
                            h2s, ch4, co2,  # gas of interest always first
                            p_or_f=p_tot,  # total pressure
                            mol_frac=[0.2, 0.6, 0.2],  # in gas order
                            temp=308)  # K

# computing mixed-gas sorption selectivity of CH4/(CO2+H2S)
p_tot, alpha, alphaerr = mg.selectivity(
                        [p_ch4, c_ch4, cerr_ch4],
                        [p_co2, c_co2, cerr_co2],
                        [p_h2s, c_h2s, cerr_h2s],
                        calc='1_numerator')  # both acid gases on bottom