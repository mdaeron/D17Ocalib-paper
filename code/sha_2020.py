#! /usr/bin/env python3
"""
Data and functions from Sha et al. (2020)
"""

D17O_VSMOW_OF_SLAP        =  0.0
D17O_VSMOW_OF_IAEA603_CO2 = -0.194
D17O_VSMOW_OF_NBS18_CO2   = -0.163

data = '''
sample	T	w_d18O	w_D17O	c_d18O	c_D17O
I-MP5	17	-7.788	0.025	33.781	-0.156
III-MP5	17	-7.717	0.023	33.629	-0.169
IV-MP5	17	-7.699	0.021	33.937	-0.162
V-MP5	17	-7.656	0.021	34.194	-0.181
III-MP9	17	-7.758	0.02	34.185	-0.171
IV-MP9	17	-7.736	0.014	34.463	-0.161
V-MP9	17	-7.697	0.02	35.297	-0.173
LH4-1	17	-7.024	0.016	35.565	-0.173
LH4-2	17	-7.060	0.014	35.341	-0.170
LH6-1	17	-7.297	0.009	35.159	-0.168
LH6-2	17	-7.200	0.012	35.946	-0.178
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split('\t')] for l in data]
data = [
	{
		k: float(v) if k != 'sample' else v
		for k,v in zip(data[0], l)
	}
	for l in data[1:]
]
