#! /usr/bin/env python3
"""
Data and functions from Roy et al. (2025)
https://doi.org/10.1016/j.chemgeo.2025.122634
"""

D17O_VSMOW_OF_SLAP = 0
D17O_VSMOW_OF_NBS19_CO2 = -0.199

starting_waters = {
	'MW': {'d18O_VSMOW':  -6.707, 'D17O_VSMOW':  0.032},
	'SW': {'d18O_VSMOW': -11.589, 'D17O_VSMOW':  0.025},
	'OW': {'d18O_VSMOW':   0.020, 'D17O_VSMOW': -0.002},
}

data = '''
     Sample,  T, N, water, final_water_d18O_VSMOW, final_water_D17O_VSMOW, CO2acid_d18O_VSMOW, CO2acid_D17O_VSMOW
Carbonate-1, 10, 2,    MW,                 -6.349,                  0.045,             38.538,             -0.175
Carbonate-2, 25, 4,    MW,                 -7.974,                  0.065,             33.632,             -0.175
Carbonate-3, 25, 2,    MW,                 -6.904,                  0.036,             33.500,             -0.151
Carbonate-4, 25, 2,    MW,                 -6.747,                  0.044,             33.419,             -0.137
Carbonate-5, 25, 6,    SW,                 -12.17,                  0.029,             28.067,             -0.172
Carbonate-6, 25, 3,    OW,                 -1.986,                  0.047,             39.683,             -0.128
Carbonate-7, 35, 2,    MW,                 -6.221,                  0.035,             31.897,             -0.159
# Carbonate-8, 50, 2,    MW,                 -6.683,                  0.049,             29.025,             -0.491
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split(',')] for l in data]
data = [
	{k: v if k in ('Sample', 'water') else float(v) for k,v in zip(data[0], l)}
	for l in data[1:]
]

for r in data:
	r['starting_water_d18O_VSMOW'] = starting_waters[r['water']]['d18O_VSMOW']
	r['starting_water_D17O_VSMOW'] = starting_waters[r['water']]['D17O_VSMOW']

data = [r for r in data if not r['Sample'].startswith('#')] # discard last samples as per the original publication
