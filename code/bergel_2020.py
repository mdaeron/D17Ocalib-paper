"""
Data and functions from Bergel et al. (2020)
https://doi.org/10.1016/j.gca.2020.02.005
"""

D17O_VSMOW_OF_SLAP        =  0.0
D17O_VSMOW_OF_IAEA603_CO2 = -0.194
D17O_VSMOW_OF_NBS18_CO2   = -0.163

data = '''
Site	T	Sample	d18O_co2_carb	D17O_co2_carb	d18O_h2o	D17O_h2o	d18O_carb
Banias	15.3	Mb1	34.862	-0.149	-7.337	0.051	24.362
Banias	15.3	Mb2	34.862	-0.151	-7.337	0.051	24.362
Banias	15.3	Thm1	35.059	-0.158	-7.337	0.051	24.557
Banias	15.3	Thm2	35.059	-0.152	-7.337	0.051	24.557
Dan	15.4	Mb1	35.109	-0.160	-7.532	0.052	24.607
Dan	15.4	Mb2	35.070	-0.160	-7.532	0.052	24.568
Dan	15.4	Mb3	34.966	-0.157	-7.532	0.052	24.465
Dan	15.4	Thm	35.110	-0.145	-7.532	0.052	24.608
Ein Tina	21.5	Mb1	34.681	-0.142	-6.685	0.043	24.183
Ein Tina	21.5	Mb2	34.500	-0.140	-6.685	0.043	24.004
Ein Tina	21.5	Thm	34.747	-0.139	-6.685	0.043	24.248
Timsah	24.5	Mb	36.107	-0.146	-4.693	0.048	25.595
Timsah	24.5	Ml	35.800	-0.135	-4.693	0.048	25.291
Ein Shokek	26.2	Mb1	35.284	-0.135	-4.846	0.051	24.780
Ein Shokek	26.2	Mb2	35.321	-0.133	-4.846	0.051	24.817
Ein Shokek	26.2	Mb3	35.321	-0.145	-4.846	0.051	24.817
Ein Shokek	26.2	Thm1	35.706	-0.161	-4.846	0.051	25.198
Ein Shokek	26.2	Thm2	35.691	-0.166	-4.846	0.051	25.183
Ein Shokek	26.2	Thm3	35.691	-0.165	-4.846	0.051	25.183
Ein Amal	27.5	Mb	35.572	-0.145	-4.694	0.048	25.065
Ein Amal	27.5	Thm	35.541	-0.140	-4.694	0.048	25.034
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split('\t')] for l in data]
data = [
	{k: v if k in ('Sample', 'Site') else float(v) for k,v in zip(data[0], l)}
	for l in data[1:]
]
