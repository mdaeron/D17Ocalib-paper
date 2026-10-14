#! /usr/bin/env python3
"""
Data and functions from O'Neil et al. (1969)
https://doi.org/10.1063/1.1671982
"""

from numpy import log, sqrt

calcite_water_samples = [
	dict(T = 700, klna18 = -0.12, klna18_cl = 0.141 / sqrt(2)),
	dict(T = 500, klna18 = 1.50, klna18_cl = 0.141 / sqrt(7)),
	dict(T = 430, klna18 = 2.33, klna18_cl = 0.141 / sqrt(2)),
	dict(T = 400, klna18 = 2.75, klna18_cl = 0.141 / sqrt(5)),
	dict(T = 370, klna18 = 3.10, klna18_cl = 0.141 / sqrt(5)),
	dict(T = 350, klna18 = 3.56, klna18_cl = 0.141 / sqrt(1)),
	dict(T = 320, klna18 = 4.27, klna18_cl = 0.141 / sqrt(2)),
	dict(T = 240, klna18 = 7.18, klna18_cl = 0.141 / sqrt(7)),
	dict(T = 201, klna18 = 9.41, klna18_cl = 0.141 / sqrt(3)),
	dict(T = 25, klna18 = 27.98, klna18_cl = 0.141 / sqrt(2)),
	dict(T = 0, klna18 = 33.94, klna18_cl = 0.141 / sqrt(1)),
	]

# correct for modern fractionation factors
# 18α_CO2/H2O = 1.0412  (Kim et al., 2015)
# 18α_cc/H2O =  1.01025
for s in calcite_water_samples:
	s['klna18'] += 1000 * log(1.0412 / 1.0407 * 1.01008 / 1.01025)
