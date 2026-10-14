#! /usr/bin/env python3
"""
Data and functions from Wostbrock et al. (2020a, 2020b)
https://doi.org/10.1016/j.chemgeo.2019.119432
https://doi.org/10.1016/j.gca.2020.07.045
"""

D17O_VSMOW_OF_SLAP    = -0.015
D17O_VSMOW_OF_IAEA603 = -0.100
D17O_VSMOW_OF_NBS18   = -0.048

gca_synthetic_with_ca = '''
 T, Nw,  w_d18O, w_d18O_SD, w_D17O, w_D17O_SD, Nc, c_d18O, c_d18O_SD, c_D17O, c_D17O_SD
 5,  4, -13.160,     0.241,  0.018,     0.012,  4, 20.611,         0, -0.053,     0.005
20,  3, -13.278,     0.202,  0.010,     0.007,  5, 17.312,         0, -0.051,     0.003
40,  2, -13.069,     0.096,  0.015,     0.008,  4, 13.394,         0, -0.043,     0.004
'''

gca_synthetic_with_ca = gca_synthetic_with_ca[1:-1].splitlines()
gca_synthetic_with_ca = [[s.strip() for s in l.split(',')] for l in gca_synthetic_with_ca]
gca_synthetic_with_ca = [
	{k: float(v) for k,v in zip(gca_synthetic_with_ca[0], l)}
	for l in gca_synthetic_with_ca[1:]
]

gca_synthetic_without_ca = '''
 T, Nw,  w_d18O, w_d18O_SD, w_D17O, w_D17O_SD, Nc, c_d18O, c_d18O_SD, c_D17O, c_D17O_SD
 5,  2, -13.583,     0.129,  0.018,     0.005,  6, 21.126,         0, -0.063,     0.005
20,  2, -13.293,     0.456,  0.021,     0.006,  2, 17.518,         0, -0.052,     0.001
40,  3, -13.268,     0.212,  0.020,     0.009,  9, 13.601,         0, -0.039,     0.009
'''

gca_synthetic_without_ca = gca_synthetic_without_ca[1:-1].splitlines()
gca_synthetic_without_ca = [[s.strip() for s in l.split(',')] for l in gca_synthetic_without_ca]
gca_synthetic_without_ca = [
	{k: float(v) for k,v in zip(gca_synthetic_without_ca[0], l)}
	for l in gca_synthetic_without_ca[1:]
]

gca_devils_hole = '''
   T, Nw,  w_d18O, w_d18O_SD, w_D17O, w_D17O_SD, Nc, c_d18O, c_d18O_SD, c_D17O, c_D17O_SD
33.7,  2, -13.388,     0.108,  0.025,     0.003,  5, 14.621,         0, -0.043,     0.007
'''

gca_devils_hole = gca_devils_hole[1:-1].splitlines()
gca_devils_hole = [[s.strip() for s in l.split(',')] for l in gca_devils_hole]
gca_devils_hole = [
	{k: float(v) for k,v in zip(gca_devils_hole[0], l)}
	for l in gca_devils_hole[1:]
]

gca_marine_calcite = '''
   T, Nw,  w_d18O, w_d18O_SD, w_D17O, w_D17O_SD, Nc, c_d18O, c_d18O_SD, c_D17O, c_D17O_SD
29.5,  0,    -0.1,         0, -0.005,         0,  3, 28.724,         0, -0.068,     0.006
 9.0,  0,   -0.08,         0, -0.005,         0,  3, 32.679,         0, -0.079,     0.009
-0.4,  0,   -0.87,         0, -0.005,         0,  5, 34.044,         0, -0.092,     0.003
-2.0,  0,   -1.21,         0, -0.005,         0,  8, 34.847,         0, -0.098,     0.005
'''

gca_marine_calcite = gca_marine_calcite[1:-1].splitlines()
gca_marine_calcite = [[s.strip() for s in l.split(',')] for l in gca_marine_calcite]
gca_marine_calcite = [
	{k: float(v) for k,v in zip(gca_marine_calcite[0], l)}
	for l in gca_marine_calcite[1:]
]

gca_marine_aragonite = '''
   T, Nw,  w_d18O, w_d18O_SD, w_D17O, w_D17O_SD, Nc, c_d18O, c_d18O_SD, c_D17O, c_D17O_SD
26.0,  0,     1.1,         0, -0.005,         0,  5, 29.762,         0, -0.080,     0.010
25.0,  0,     1.5,         0, -0.005,         0,  3, 30.889,         0, -0.084,     0.002
25.0,  0,     1.5,         0, -0.005,         0,  4, 30.904,         0, -0.075,     0.008
25.0,  0,    1.25,         0, -0.005,         0,  3, 31.136,         0, -0.084,     0.005
'''

gca_marine_aragonite = gca_marine_aragonite[1:-1].splitlines()
gca_marine_aragonite = [[s.strip() for s in l.split(',')] for l in gca_marine_aragonite]
gca_marine_aragonite = [
	{k: float(v) for k,v in zip(gca_marine_aragonite[0], l)}
	for l in gca_marine_aragonite[1:]
]
