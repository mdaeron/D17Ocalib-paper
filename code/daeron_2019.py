#! /usr/bin/env python3
"""
Data and functions from Daëron et al. (2019)
https://doi.org/10.1038/s41467-019-08336-5
"""

def ln_alpha18_calcite_water(T):
	return 17570 / T - 29.13

samples = {
	'LGB': {
		'T': 7.9,
		'SE_T': 0.2,
		'klna18': 33.38,
		'SE_klna18': 0.10,
	},
	'DVH': {
		'T': 33.7,
		'SE_T': 0.2,
		'klna18': 28.13,
		'SE_klna18': 0.06,
	}
}
