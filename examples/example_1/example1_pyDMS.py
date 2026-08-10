import numpy as np
import pyDMS.dms as dms

# Data from:
# Dean, P. A.; Mizrahi Rodriguez, K.; Guo, S.; Roy, N.; Swager, T. M.; Smith,
# Z. P. Elucidating the role of micropore-generating backbone motifs and amine
# functionality on H2S, CO2, CH4 and N2 sorption. Journal of Membrane Science
# 2024, 696, 122465.

co2 = dms.Gas()  # setting up the Class to hold all the data

co2.formula = "CO2"  # Formula of the penetrant

# fugacity data (atm)
co2.f = [np.array([0.2, 0.5, 1.1, 1.8, 2.7, 3.7, 4.7, 6.1, 8.8, 11.8, 16.8, 22.1, 27.1,
                   31.9, 36.3]),
         np.array([0.2, 0.6, 1.3, 2.0, 2.8, 3.7, 4.7, 6.2, 8.9, 11.8, 16.9, 22.3, 27.6,
                   32.5, 37.1]),
         np.array([0.2, 0.7, 1.3, 1.9, 2.8, 3.8, 4.8, 6.2, 8.9, 11.8, 17.0, 22.6, 28.0,
                   33.1, 37.9]),
         np.array([0.2, 0.6, 1.2, 1.9, 2.8, 3.8, 4.8, 6.2, 8.9, 11.9, 17.2, 22.8, 28.3,
                   33.5, 38.5])]

# concentration data (cm^3 cm^-3)
co2.c = [np.array([8.9, 18.4, 30.9, 41.0, 51.3, 59.8, 66.9, 75.6, 88.6, 101.7, 122.1,
                   140.6, 156.2, 168.5, 175.3]),
         np.array([6.5, 16.5, 27.1, 35.7, 43.6, 51.4, 57.8, 65.6, 77.4, 88.5, 105.3,
                   122.3, 135.8, 147.4, 156.6]),
         np.array([5.1, 13.4, 22.5, 28.9, 37.1, 44.5, 50.5, 57.3, 68.1, 77.7, 92.8,
                   107.3, 119.7, 131.0, 141.0]),
         np.array([4.4, 9.6, 17.4, 24.5, 31.4, 38.3, 43.2, 49.7, 60.0, 69.2, 83.1,
                   96.2, 107.5, 117.4, 126.0])]

# uncertainty in concentration data
co2.c_err = [np.array([0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.9, 1.2, 1.5, 1.7,
                       1.9, 2.0]),
             np.array([0.1, 0.2, 0.3, 0.3, 0.4, 0.5, 0.5, 0.6, 0.8, 0.9, 1.1, 1.4, 1.6,
                       1.8, 1.9]),
             np.array([0.1, 0.2, 0.3, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 1.7,
                       1.9, 2.0]),
             np.array([0.1, 0.1, 0.2, 0.2, 0.3, 0.4, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4,
                       1.5, 1.6])]

# temperature (K) of each sorption isotherm
co2.temp = np.array([308, 318, 328, 338])

# run optimization
dms.compute(co2, 'example1_pyDMS')
