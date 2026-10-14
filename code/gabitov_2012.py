#! /usr/bin/env python3
"""
Data from Gabitov et al. (2012)
https://doi.org/10.1016/j.chemgeo.2012.02.021
"""

from numpy import log

d18Ow = -10.5

calcite_water_samples = [
	dict(run = 'DC-15', T = 15.0, D18O_slow = 31.38, D18O_center = 28.49),
	dict(run = 'DSC-4', T = 19.7, D18O_slow = 30.19, D18O_center = 28.17),
	dict(run = 'DC-3',  T = 21.7, D18O_slow = 30.12, D18O_center = 28.49),
	dict(run = 'DC-1',  T = 24.6, D18O_slow = 29.33, D18O_center = 27.59),
	]

for r in calcite_water_samples:
	r['klna18_slow'] = 1000 * log((1000 + d18Ow + r['D18O_slow']) / (1000 + d18Ow))
