#! /usr/bin/env python3
"""
Data and functions from Hayles et al. (2018)
https://doi.org/10.1016/j.gca.2018.05.032
"""

from numpy import dot

ln_beta18_calcite_coefs = [
	 1.331245e-4,
	-5.393675e-1,
	 1.307870e4,
	-4.495755e5,
	-5.417531e8,
	 1.463936e11,
	-1.633009e13,
	 7.027321e14,
]

def ln_beta18_calcite(T):
	'''
	Polynomial fit of ln(18β) for calcite
	'''
	return (
		ln_beta18_calcite_coefs[0]
		+ ln_beta18_calcite_coefs[1] / T
		+ ln_beta18_calcite_coefs[2] / T**2
		+ ln_beta18_calcite_coefs[3] / T**3
		+ ln_beta18_calcite_coefs[4] / T**4
		+ ln_beta18_calcite_coefs[5] / T**5
		+ ln_beta18_calcite_coefs[6] / T**6
		+ ln_beta18_calcite_coefs[7] / T**7
	)

ln_beta18_water_coefs = [
	-7.839005e-3,
	 5.689300,
	 2.143196e4,
	-1.353889e7,
	 5.090782e9,
	-1.114055e12,
	 1.333519e14,
	-6.705843e15,
]

def ln_beta18_water(T):
	'''
	Polynomial fit of ln(18β) for water
	'''
	return (
		ln_beta18_water_coefs[0]
		+ ln_beta18_water_coefs[1] / T
		+ ln_beta18_water_coefs[2] / T**2
		+ ln_beta18_water_coefs[3] / T**3
		+ ln_beta18_water_coefs[4] / T**4
		+ ln_beta18_water_coefs[5] / T**5
		+ ln_beta18_water_coefs[6] / T**6
		+ ln_beta18_water_coefs[7] / T**7
		)

ln_beta18_water_semiempirical_coefs = [
	-1.168952e-2,
	 1.764688e1,
	 6.057270e3,
	-2.952406e6,
	 8.761135e8,
	-1.428610e11,
	 1.318367e13,
	-5.324090e14,
]

def ln_beta18_water_semiempirical(T):
	'''
	Polynomial fit of ln(18β) for water
	'''
	return (
		ln_beta18_water_semiempirical_coefs[0]
		+ ln_beta18_water_semiempirical_coefs[1] / T
		+ ln_beta18_water_semiempirical_coefs[2] / T**2
		+ ln_beta18_water_semiempirical_coefs[3] / T**3
		+ ln_beta18_water_semiempirical_coefs[4] / T**4
		+ ln_beta18_water_semiempirical_coefs[5] / T**5
		+ ln_beta18_water_semiempirical_coefs[6] / T**6
		+ ln_beta18_water_semiempirical_coefs[7] / T**7
		)

def ln_alpha18_calcite_water(T):
	return ln_beta18_calcite(T) - ln_beta18_water(T)

def ln_alpha18_calcite_water_semiempirical(T):
	return ln_beta18_calcite(T) - ln_beta18_water_semiempirical(T)

def kappa_calcite(T):
	return (
		1.019124e9 / T**5
		- 2.117501e7 / T**4
		+ 1.686453e5 / T**3
		- 5.784679e2 / T**2
		+ 1.489666e-1 / T
		+ 0.5304852
		)

def kappa_water(T):
	return (
		7.625734e6 / T**5
		+ 1.216102e6 / T**4
		- 2.135774e4 / T**3
		+ 1.323782e2 / T**2
		- 4.931630e-1 / T
		+ 0.5306551
		)

def kappa_water_semiempirical(T):
	return (
		- 2.993068e8 / T**5
		+ 6.507712e6 / T**4
		- 5.962566e4 / T**3
		+ 2.676809e2 / T**2
		- 7.097807e-1 / T
		+ 0.5308002
		)

def theta_calcite_water(T):
	return kappa_calcite(T) + (kappa_calcite(T) - kappa_water(T)) * ln_beta18_water(T) / ln_alpha18_calcite_water(T)

def theta_calcite_water_semiempirical(T):
	return kappa_calcite(T) + (kappa_calcite(T) - kappa_water_semiempirical(T)) * ln_beta18_water_semiempirical(T) / ln_alpha18_calcite_water_semiempirical(T)

def ln_alpha17_calcite_water(T):
	return theta_calcite_water(T) * ln_alpha18_calcite_water(T)

def ln_alpha17_calcite_water_semiempirical(T):
	return theta_calcite_water_semiempirical(T) * ln_alpha18_calcite_water_semiempirical(T)

if __name__ == '__main__':

	from pylab import *
	
	t = logspace(3,9)
	y = theta_calcite_water(t+273.15)
	semilogx(t, y, 'r-')
	show()

	T = linspace(0, 700, 701) + 273.15
	X = 1/T
	Y = ln_alpha18_calcite_water(T) - ln_alpha18_calcite_water_semiempirical(T)	
	plot(T-273.15, Y, 'r-')
	show()


	subplot(121)
	T = linspace(0, 700, 701) + 273.15
	X = 1/T
	Y = ln_alpha18_calcite_water(T)	
	plot(X, Y, 'k-')
	Y = ln_alpha18_calcite_water_semiempirical(T)	
	plot(X, Y, 'k-')
	
	T = array([7.9,33.7]) + 273.15
	X = 1/T
	Y = 17.57/T - 29.13e-3
	plot(X, Y, 'r+', ms = 10)

	subplot(122)
	T = linspace(0, 700, 701) + 273.15
	X = 1/T
	Y = theta_calcite_water(T)
	plot(X, Y, 'k-')

	T = array([7.9,33.7]) + 273.15
	X = 1/T
	Y = theta_calcite_water(T)
	plot(X, Y, 'r+', ms = 10)
	
	show()
	
	subplot(121)

	T = linspace(0, 250, 251) + 273.15
	X = 1/T
	Y = theta_calcite_water(T)
	plot(X, Y, 'k-')
	Y = 0.5305 - 1.39/T
	plot(X, Y, 'c-')

	T = array([7.9,33.7]) + 273.15
	X = 1/T
	Y = theta_calcite_water(T)
	plot(X, Y, 'r+', ms = 10)
	
	subplot(122)

	T = linspace(0, 1650, 1651) + 273.15
	X = 1/T
	Y = (theta_calcite_water(T) - 0.528) * 1000 * ln_alpha18_calcite_water(T)
	plot(X, Y, 'k-')
	Y = (0.5305 - 1.39/T - 0.528) * 1000 * ln_alpha18_calcite_water(T)
	plot(X, Y, 'c-')

	show()

