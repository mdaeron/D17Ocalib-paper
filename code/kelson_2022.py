#! /usr/bin/env python3
"""
Data and functions from Kelson et al. (2022)
"""

D17O_VSMOW_OF_SLAP    =  0.0
D17O_VSMOW_OF_IAEA603 = -0.100

data = '''
sample	T	ar	cc	Nw	w_d18O	w_D17O	Nc	c_d18O	c_D17O
Scott Creek mussel	18	0.7	0.3	2	-0.611	-0.0081	3	31.446	-0.1040
San Simeon mussel	15	0.63	0.37	3	0.447	-0.0075	3	31.405	-0.0937
Silver Lake clam	19	1	0	3	-2.920	-0.0189	3	27.549	-0.1111
Huron  River mussel	13	1	0	3	-8.120	0.0152	2	23.520	-0.0905
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split('\t')] for l in data]
data = [
	{k: (
		float(v) if k != 'sample'
		else (
			int(v) if k.startswith('N')
			else v
		)
	)
		for k,v in zip(data[0], l)}
	for l in data[1:]
]
