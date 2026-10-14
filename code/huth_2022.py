#! /usr/bin/env python3
"""
Data and functions from Huth et al. (2022)
"""

D17O_VSMOW_OF_SLAP    =  0.0
D17O_VSMOW_OF_IAEA603 = -0.100
TENTATIVE_D17O_VSMOW_OF_NBS18 = (
	-0.100     # IAEA-603 value used above
	- (-0.127) + (-0.130) # IAEA-603 -> NBS19 after Chaillot et al. (2025)
	- (-0.135) + (-0.098) # NBS19 -> NBS18 after Passey et al. (2014)
)

data = '''
experiment	T	water	Nw	w_d18O	w_d18O_se	w_D17O	w_D17O_se	carb	Nc	c_d18O	c_d18O_se	c_D17O	c_D17O_se
5.3	5	D17OEQ 5-3 3.17.2020	5	-6.786	0.1	0.006	0.004	EQ17O-5-3	4	26.656	0.11	-0.091	0.002
5.4	5	D17OEQ 5-4 3.17.2020	2	-6.595	0.1	0.010	0.004	EQ17O-5-4	4	26.943	0.11	-0.090	0.003
15.1	15	D17OEQ 15-1 3.2.2020	5	-6.504	0.1	0.009	0.002	EQ17O-15-1	4	24.053	0.11	-0.090	0.002
35.1	35	D17OEQ 35-1 3.4.2020	5	-6.471	0.1	0.008	0.003	EQ17O-35-1	3	20	0.11	-0.066	0.003
'''

data = data[1:-1].splitlines()
data = [[s.strip() for s in l.split('\t')] for l in data]
data = [
	{k: (
		float(v) if not (k in ('experiment', 'water', 'carb') or k.startswith('N'))
		else (
			int(v) if k.startswith('N')
			else v
		)
	)
		for k,v in zip(data[0], l)}
	for l in data[1:]
]

renormalized_passey_2014_data = '''
sample	T	w_d18O	w_D17O	c_d18O	c_D17O	c_D17O_se
JHUEGG001	39	-3.416	-0.051	27.137	-0.133	0.002
JHUEGG002	39	-4.806	-0.018	24.479	-0.120	0.010
USA2013BALchicken003	39	-3.681	-0.042	24.607	-0.113	0.003
USA2012BALchicken004	39	-4.686	-0.010	22.989	-0.092	0.003
JHUEGG007	39	-5.195	-0.006	23.927	-0.100	0.002
HPL13	24	-4.203	0.016	24.669	-0.080	0.002
JBC03coral	26	1.913	0.000	26.318	-0.087	0.008
'''

renormalized_passey_2014_data = renormalized_passey_2014_data[1:-1].splitlines()
renormalized_passey_2014_data = [[s.strip() for s in l.split('\t')] for l in renormalized_passey_2014_data]
renormalized_passey_2014_data = [
	{
		k: v if k == 'sample' else float(v)
		for k,v in zip(renormalized_passey_2014_data[0], l)
	}
	for l in renormalized_passey_2014_data[1:]
]
