"""
Data and functions from Wostbrock et al. (2020a, 2020b)
https://doi.org/10.1016/j.chemgeo.2019.119432
https://doi.org/10.1016/j.gca.2020.07.045
"""

D17O_VSMOW_OF_SLAP        =  0.0
D17O_VSMOW_OF_IAEA603_CO2 = -0.194
D17O_VSMOW_OF_NBS18_CO2   = -0.163

data = '''
T	d18O_co2_carb	D17O_co2_carb	d18O_h2o	D17O_h2o	d18O_carb
10	44.866	-0.202	1.134	0.020	34.265
10	44.892	-0.214	0.913	0.015	33.949
10	44.649	-0.194	0.799	0.009	33.665
10	44.774	-0.215	1.21	0.020	33.784
27	41.577	-0.2	1.493	0.015	31.009
27	41.555	-0.208	1.456	0.001	30.612
35	38.341	-0.197	0.11	0.022	27.806
35	38.367	-0.195	0.229	0.041	27.8
35	38.67	-0.198	0.249	0.027	27.826
27	39.012	-0.185	0.400	0.001	28.453
27	39.683	-0.201	0.38	0.008	28.99
27	40.211	-0.185	0.493	0.031	29.27
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split('\t')] for l in data]
data = [
	{k: v if k in ('Sample', 'water') else float(v) for k,v in zip(data[0], l)}
	for l in data[1:]
]
