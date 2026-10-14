import lmfit
import config
import tomllib
import guo_zhou_2019
import huth_2022
import kelson_2022
import wostbrock_2020
import voarintsoa_2020
import bergel_2020
import roy_2025
import sha_2020

import numpy as np
import pandas as pd

from pathlib import Path
from config import vprint
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import FancyArrowPatch

CO2g_H2O_fractionations = {
	'Guo & Zhou (2019)': dict(
		alpha18_CO2g_H2O = 1.041461382,
		theta17_CO2g_H2O = 0.524606996,
	),
	'Barkan & Luz (2012)': dict(
		alpha18_CO2g_H2O = 1.041036,
		theta17_CO2g_H2O = 0.5229,
	),
	'Roy et al. (2025)': dict(
		alpha18_CO2g_H2O = 1.04122,
		theta17_CO2g_H2O = 0.5230,
	),
}

vprint(f'Starting {Path(__file__).name}')

with open('output/labs-results.toml', 'w') as f:
	f.write('[DD17O_CO2_H2O_25C]\n')
	for k,v in CO2g_H2O_fractionations.items():
		v['DD17O'] = (v['theta17_CO2g_H2O'] - 0.528) * np.log(v['alpha18_CO2g_H2O']) * 1e3
		f.write(f'"{k}" = {v["DD17O"]:.6f}\n')
	f.write('\n')

fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, config.ONE_COLUMN_WIDTH))
fig.subplots_adjust(.22, .19, .98, .95)
ax = ppl.subplot(111)
ax.axis([-0.200, -0.010, -0.230, -0.040])
ax.set_xlabel(r'$Δ^{17}O^\text{NBS18}_\text{VSMOW}$ (‰)')
ax.set_ylabel(r'$Δ^{17}O^\text{IAEA-603}_\text{VSMOW}$ (‰)')
ax.xaxis.set_major_locator(MultipleLocator(5e-2))
ax.yaxis.set_major_locator(MultipleLocator(5e-2))

gas = dict(
	ls = 'None',
	marker = 'o',
	ms = 6,
	zorder = 1000,
)

carb = dict(
	ls = 'None',
	marker = 'D',
	ms = 5,
	zorder = 1000,
)

errb = config.kw_errorbar | dict(marker = 'None')

arrstyle = dict(
	arrowstyle = '->,head_width=.15',
	mutation_scale = 12,
	shrinkA = 16,
	shrinkB = 16,
	connectionstyle = "arc3,rad=-0.4",
	lw = 0.8,
)

with open('output/DD17O-results.toml', 'rb') as f:
	_ = tomllib.load(f)
	D17O_VSMOW_OF_IAEA603 = _['bestfit']['D17O_VSMOW_OF_IAEA603']['value']
	SIGMA_D17O_VSMOW_OF_IAEA603 = _['bestfit']['D17O_VSMOW_OF_IAEA603']['se']
	D17O_VSMOW_OF_NBS18 = _['bestfit']['D17O_VSMOW_OF_NBS18']['value']
	SIGMA_D17O_VSMOW_OF_NBS18 = _['bestfit']['D17O_VSMOW_OF_NBS18']['se']
	COEFS = _['bestfit']['D17O_VSMOW_OF_IAEA603']['coefs']
	ALTD17O_VSMOW_OF_IAEA603 = _['altbestfit']['D17O_VSMOW_OF_IAEA603']['value']
	ALTSIGMA_D17O_VSMOW_OF_IAEA603 = _['altbestfit']['D17O_VSMOW_OF_IAEA603']['se']
	ALTD17O_VSMOW_OF_NBS18 = _['altbestfit']['D17O_VSMOW_OF_NBS18']['value']
	ALTSIGMA_D17O_VSMOW_OF_NBS18 = _['altbestfit']['D17O_VSMOW_OF_NBS18']['se']
	ALTCOEFS = _['altbestfit']['D17O_VSMOW_OF_IAEA603']['coefs']

xl, yl = 0.79, 0.11
dx, dy = -0.07, 0.06

X = D17O_VSMOW_OF_NBS18
Y = D17O_VSMOW_OF_IAEA603
eX = SIGMA_D17O_VSMOW_OF_NBS18 * 1.96
eY = SIGMA_D17O_VSMOW_OF_IAEA603 * 1.96
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.thisstudy.D17O_VSMOW_of_NBS18]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.thisstudy.D17O_VSMOW_of_IAEA603]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')
ref = 'this study'
color = np.array([0, 0, 0])
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = 1-(1-color)/1.5, mfc = 1-(1-color)/3, **carb, label = f'$CaCO3$, {ref}')

ax.plot(xl, yl, mec = 1-(1-color)/1.5, mfc = 1-(1-color)/3, **carb, transform = ax.transAxes)
ax.text(xl+dx/4, yl + dy/2, ref, size = 7, va = 'bottom', ha = 'left', transform = ax.transAxes, rotation = 45)
xl += dx

ref = 'Ellis & Passey (2023)'
color = np.array(config.COLORS['Ellis'])
X = config.SLAP_adjust('NBS18', ref)
Y = config.SLAP_adjust('IAEA-603', ref)
eX = 0.005 * 1.96
eY = 0.005 * 1.96
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.Ellis2023.D17O_VSMOW_of_NBS18]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.Ellis2023.D17O_VSMOW_of_IAEA603]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = color, mfc = 1-(1-color)/3, **carb, label = f'$CaCO3$, {ref}')
ax.plot(xl, yl, mec = color, mfc = 1-(1-color)/3, **carb, transform = ax.transAxes)
ax.text(xl+dx/4, yl + dy/2, ref, size = 7, va = 'bottom', ha = 'left', transform = ax.transAxes, rotation = 45)
xl += dx

ref = 'Wostbrock et al. (2020)'
color = np.array(config.COLORS['Wostbrock'])
X = config.SLAP_adjust('NBS18', ref)
eX = 0.009 * 1.96
Y = config.SLAP_adjust('IAEA-603', ref)
eY = 0.007 * 1.96
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.Wostbrock2020.D17O_VSMOW_of_NBS18]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.Wostbrock2020.D17O_VSMOW_of_IAEA603]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = color, mfc = 1-(1-color)/3, **carb, label = f'$CaCO3$, {ref}')
ax.plot(xl, yl, mec = color, mfc = 1-(1-color)/3, **carb, transform = ax.transAxes)
ax.text(xl+dx/4, yl + dy/2, ref, size = 7, va = 'bottom', ha = 'left', transform = ax.transAxes, rotation = 45)

X = config.SLAP_adjust('NBS18-CO2', ref)
eX = 0.0 * 1.96
Y = config.SLAP_adjust('IAEA-603-CO2', ref)
eY = 0.005 * 1.96
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.Wostbrock2020.D17O_VSMOW_of_NBS18_CO2]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.Wostbrock2020.D17O_VSMOW_of_IAEA603_CO2]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
for _ in (-1, 1):
	ax.plot(
		[X, X+0.011*1.96*_], [Y, Y], '-',
		dashes = (2.5,1.5),
		lw = errb['elinewidth'],
		color = 1-(1-color)/3,
	)
	ax.text(X-0.011*2.4*_, Y, '?',
		color = 1-(1-color)/3,
		va = 'center',
		ha = 'center',
		weight = 'bold',
		size = 7,
	)

ax.plot(X, Y, mec = color, mfc = 1-(1-color)/3, **gas, label = f'$CO_2$, {ref}')
ax.plot(xl, yl-dy, mec = color, mfc = 1-(1-color)/3, **gas, transform = ax.transAxes)
xl += dx

ref = 'Barkan et al. (2019)'
color = np.array(config.COLORS['Barkan'])

X = config.SLAP_adjust('NBS18-CO2', ref)
Y = config.SLAP_adjust('IAEA-603-CO2', ref)
eX = 0.005 * 1.96
eY = 0.006 * 1.96
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = 1-(1-color)/3, mfc = 1-(1-color)/3, **gas, label = f'$CO_2$, {ref}')
ax.text(xl+dx/4, yl + dy/2, ref, size = 7, va = 'bottom', ha = 'left', transform = ax.transAxes, rotation = 45)

_X, _Y = X, Y
X -= CO2g_H2O_fractionations['Barkan & Luz (2012)']['DD17O']
X += CO2g_H2O_fractionations['Guo & Zhou (2019)']['DD17O']
Y -= CO2g_H2O_fractionations['Barkan & Luz (2012)']['DD17O']
Y += CO2g_H2O_fractionations['Guo & Zhou (2019)']['DD17O']
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.Barkan2019.D17O_VSMOW_of_NBS18_CO2]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.Barkan2019.D17O_VSMOW_of_IAEA603_CO2]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')

ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = color, mfc = 1-(1-color)/3, **gas)
ax.plot(xl, yl-dy, mec = color, mfc = 1-(1-color)/3, **gas, transform = ax.transAxes)
arr = FancyArrowPatch((_X, _Y), (X, Y), **arrstyle)
ax.add_patch(arr)
ax.annotate(
	'switch to\n${}_{}$different${}_{}$\n$Δ^{17}O^{eq}_{CO_2/H_2O}$',
	(-0.22, 1.0),
	xycoords = arr,
	ha = 'center',
	va = 'center',
	size = 8,
	linespacing = 1.3,
)

_X, _Y = X, Y
# X += 0.05 # Liu & Cao (2026)
# Y += 0.05 # Liu & Cao (2026)
X -= (0.5230 - 0.528) * np.log(1.01025) * 1e3 # Wostbrock et al. (2020, Chemical Geology)
Y -= (0.5230 - 0.528) * np.log(1.01025) * 1e3 # Wostbrock et al. (2020, Chemical Geology)
with open('output/labs-results.toml', 'a') as f:
	f.write(f'[reprocessed.Barkan2019.D17O_VSMOW_of_NBS18]\n')
	f.write(f'value = {X:6f}\n')
	f.write(f'c95 = {eX:6f}\n\n')
	f.write(f'[reprocessed.Barkan2019.D17O_VSMOW_of_IAEA603]\n')
	f.write(f'value = {Y:6f}\n')
	f.write(f'c95 = {eY:6f}\n\n')
ax.errorbar(X, Y, eY, eX, ecolor = 1-(1-color)/3, **errb)
ax.plot(X, Y, mec = 1-(1-color)/2, mfc = 1-(1-color)/3, **carb)
ax.plot(xl, yl, mec = 1-(1-color)/2, mfc = 1-(1-color)/3, **carb, transform = ax.transAxes)
xl += dx

arrstyle['arrowstyle'] = '<' + arrstyle['arrowstyle']
arrstyle['shrinkA'] = 13
arrstyle['shrinkB'] = 11
arr = FancyArrowPatch((_X, _Y), (X, Y), **arrstyle)
ax.add_patch(arr)
ax.annotate(
	'$Δ^{17}O^{acid}_{CO_2/CaCO_3}$',
	(-0.4, 0.95),
	xycoords = arr,
	ha = 'center',
	va = 'center',
	size = 8,
)

ax.text(xl-dx/3, yl, '         $CaCO_3$', size = 7, va = 'center', ha = 'right', transform = ax.transAxes)
ax.text(xl-dx/3, yl-dy, ' $CO_2$', size = 7, va = 'center', ha = 'right', transform = ax.transAxes)

fig.savefig('output/labs.pdf')
ppl.close(fig)












### ---------------------------------------------------------------
### ---------------------------------------------------------------
### ---------------------------------------------------------------




synthetic_marker = ('s', 4)
marine_cc_marker = ('d', 5)
marine_ar_marker = (config.hdiamond, 5)
biogenic_marker = ('D', 4)
speleo_marker = ('v', 5)
dvh_marker = ('o', 4)

for renormalize, alternative_D17O_603_of_NBS18 in (
	(False, None),
	(True, config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)') - config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')),
	(True, False),
):
		fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, config.ONE_COLUMN_WIDTH*1.2))
		fig.subplots_adjust(.22, .45, .99, .91)
		ax = ppl.subplot(111)
		# ax.axis([-0.200, -0.000, -0.230, -0.030])
		# ax.set_xlabel(r'                                                  $1\,/\,T^2$')
		ax.set_ylabel(r'$Δ^{17}O^\text{carbonate}_\text{IAEA-603} - Δ^{17}O^\text{water}_\text{VSMOW}$ (‰)')
		ax.text(1-0.02, 0-0.02, r'$1/\,T^2$', transform = ax.transAxes, va = 'top', ha = 'right', size = 9)
		ax.yaxis.set_major_locator(MultipleLocator(50e-3))

		tmin = -6
		tmax = 230
		twinax = ax.twiny()

		Tl = np.array([0, 20, 50, 100, 200])[::-1]
		twinax.xaxis.set_ticks((Tl+273.15)**-2)
		twinax.xaxis.set_ticklabels(
			[f'{t:.0f}$\\,$°C$\\hspace{{-0.5}}$' for t in Tl],
			rotation = -30,
			ha = 'right',
			va = 'bottom',
			size = 9,
			linespacing = 1,
		)
		ax.tick_params(axis = 'x', length = 0, labelbottom = False)
		ax.grid(False, axis = 'x')
		if alternative_D17O_603_of_NBS18:
			ax.text(
				0.33, 0.67,
				'ALTERNATIVE\n' r'$Δ^{17}O^\text{NBS18}_\text{IAEA-603}$ = ' f'{alternative_D17O_603_of_NBS18*1000:.0f} ppm',
				color = 'r',
				alpha = 0.3,
				size = 12,
				rotation = 30,
				transform = ax.transAxes,
				va = 'center',
				ha = 'center',
			)

		# --------------------------------------------------- Read our data
		data = pd.read_csv('output/calib_results.csv', index_col = 'Sample')
		data['1/T'] = (data['T'] + 273.15)**-1
		data['1/T2'] = (data['T'] + 273.15)**-2

		# ---------------------------------------------------- Fit our data
		T = data['T'].to_numpy() + 273.15
		D17O_603 = data['D17O_603'].to_numpy()
		if alternative_D17O_603_of_NBS18:
			# Renormalize carbonates to postulated NBS18 composition
			D17O_603 = config.D17O_adjust_carb(
				(1000 + data['d18O_VPDB']) * 1.03092 - 1000,
				D17O_603,
				D17O_IAEA603_old = 0.,
				D17O_NBS18_old = config.D17O_VSMOW_OF_NBS18_CO2_CHAILLOT - config.D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT,
				D17O_IAEA603_new = 0.,
				D17O_NBS18_new = alternative_D17O_603_of_NBS18,
			)

		# Recompute DD17O based on renormalized values
		DD17O_obs = (D17O_603 - data['D17O_VSMOW']).to_numpy()
		DD17O_obs_sigma = data['95CL_DD17O'].to_numpy()
		DD17O_theory = np.array([
			(guo_zhou_2019.calcite_water_eq_theta(_) - 0.528)
			* guo_zhou_2019.calcite_water_eq_klna18(_)
			for _ in T
		])

		if config.ONLY_CONSIDER_NATURAL_SAMPLES_FOR_REGRESSION:
			data['klna18_mismatch'] = 1e99
			data.loc[['LGB', 'DVH'], 'klna18_mismatch'] = 1

		klna18_mismatch = data['klna18_mismatch'].to_numpy()

		def _residuals(p):
			return (DD17O_obs - DD17O_theory - p['offset']) / DD17O_obs_sigma / klna18_mismatch

		_p = lmfit.Parameters()
		_p.add('offset', value = 0.07)
		fitresult = lmfit.minimize(
			_residuals,
			_p,
			method = 'least_squares',
			scale_covar = False,
		)

		D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY = -fitresult.uvars['offset'].n
		D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY = (
			D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			- config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
			+ config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')
		)
		DD17Oeq_afo_T = lambda t: (
			(guo_zhou_2019.calcite_water_eq_theta(t) - 0.528)
			* guo_zhou_2019.calcite_water_eq_klna18(t)
			- D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
		)

		reprocessed = {}
		mew = 0.7

		# ------------------------------------------------ Plot our data
		_data = data.loc[['DVH', 'LGB', 'CC170', 'CC200']]

		D17O_603 = _data['D17O_603'].to_numpy()
		if renormalize and alternative_D17O_603_of_NBS18:
			# Renormalize carbonates to postulated NBS18 composition
			D17O_603 = config.D17O_adjust_carb(
				(1000 + _data['d18O_VPDB']) * 1.03092 - 1000,
				D17O_603,
				D17O_IAEA603_old = 0.,
				D17O_NBS18_old = config.D17O_VSMOW_OF_NBS18_CO2_CHAILLOT - config.D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT,
				D17O_IAEA603_new = 0.,
				D17O_NBS18_new = alternative_D17O_603_of_NBS18,
			)

		# Recompute DD17O based on renormalized values
		_data['DD17O'] = (D17O_603 - _data['D17O_VSMOW']).to_numpy()

		ax.errorbar(
			_data['1/T2'],
			_data['DD17O'],
			_data['95CL_DD17O'],
			zorder = 1000,
			**(config.kw_errorbar | dict(marker = 'None', capsize = 2, ecolor = 'k')),
			)

		__data = _data.loc[['DVH', 'LGB']]
		ax.plot(
			__data['1/T2'],
			__data['DD17O'],
			'wo',
			mec = 'k',
			zorder = 2000,
			label = 'this study (natural calcites)',
			ms = 4,
			)

		__data = _data.loc[['CC170', 'CC200']]
		ax.plot(
			__data['1/T2'],
			__data['DD17O'],
			'ws',
			mec = 'k',
			zorder = 1000,
			label = 'this study (synthetic calcites near $^{18}O$ equilibrium)',
			ms = 4,
			)

		reprocessed['this study'] = dict(
			symbol = dvh_marker,
			mec = 'k',
			mfc = 'w',
			dw = data['d18O_VSMOW'],
			Dw = data['D17O_VSMOW'],
			dc = data['d18O_VPDB'],
			Dc = data['D17O_603'],
		)

		# ----------------------- Plot theoretical model wih best-fit offset
		Ti = np.linspace(273.15**-2, 1e-99, 101)**-0.5

		ax.plot(
			Ti**-2,
			DD17Oeq_afo_T(Ti),
			'k-',
			lw = 0.7,
			dashes = (6,2,2,2),
			zorder = 1,
			label = 'best-fit Guo & Zhou (2019) trend',
			)

		# ax.plot(
		# 	Ti**-2,
		# 	COEFS['a0'] + COEFS['a1']/Ti + COEFS['a2']/Ti**2,
		# 	'r+',
		# 	zorder = 10000,
		# )

		# ax.plot(
		# 	Ti**-2,
		# 	ALTCOEFS['a0'] + ALTCOEFS['a1']/Ti + ALTCOEFS['a2']/Ti**2,
		# 	'b+',
		# 	zorder = 10000,
		# )


		# ------------------------------------------ Plot renormalized Passey et al. (2014)
		col = np.array(config.COLORS['Passey'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				biogenic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				huth_2022.renormalized_passey_2014_data,
				'Passey et al. (2014) biogenic',
			),
		):
			# water data
			dw = np.array([r['w_d18O'] for r in _data])
			Dw = np.array([r['w_D17O'] for r in _data])
			# carb data
			dc = np.array([r['c_d18O'] for r in _data])
			Dc = np.array([r['c_D17O'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = huth_2022.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = huth_2022.D17O_VSMOW_OF_IAEA603,
					D17O_NBS18_old = huth_2022.TENTATIVE_D17O_VSMOW_OF_NBS18,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = huth_2022.D17O_VSMOW_OF_IAEA603

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ---------------------------------------------------------------- Plot Voarintsoa et al. (2020)
		col = np.array(config.COLORS['Barkan'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				synthetic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				voarintsoa_2020.data,
				'Voarintsoa et al. (2020) synthetic carbonates',
			),
		):
			# water data
			dw = np.array([r['d18O_h2o'] for r in _data])
			Dw = np.array([r['D17O_h2o'] for r in _data])
			# carb data
			dc = np.array([r['d18O_carb'] for r in _data])
			Dc = np.array([r['D17O_co2_carb'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = voarintsoa_2020.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = voarintsoa_2020.D17O_VSMOW_OF_IAEA603_CO2,
					D17O_NBS18_old =   voarintsoa_2020.D17O_VSMOW_OF_NBS18_CO2,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = voarintsoa_2020.D17O_VSMOW_OF_IAEA603_CO2

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ---------------------------------------------------------------- Plot Bergel et al. (2020)
		col = np.array(config.COLORS['Barkan'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				marine_ar_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				bergel_2020.data,
				'Bergel et al. (2020) biogenic aragonite',
			),
		):
			# water data
			dw = np.array([r['d18O_h2o'] for r in _data])
			Dw = np.array([r['D17O_h2o'] for r in _data])

			# carb data
			dc = np.array([r['d18O_carb'] for r in _data])
			Dc = np.array([r['D17O_co2_carb'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = bergel_2020.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = bergel_2020.D17O_VSMOW_OF_IAEA603_CO2,
					D17O_NBS18_old =   bergel_2020.D17O_VSMOW_OF_NBS18_CO2,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = bergel_2020.D17O_VSMOW_OF_IAEA603_CO2

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ---------------------------------------------------------------- Plot Wostbrock et al. (2020)
		col = np.array(config.COLORS['Wostbrock'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				dvh_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				wostbrock_2020.gca_devils_hole,
				'Wostbrock et al. (2020) Devils Hole',
			),
			(
				synthetic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				wostbrock_2020.gca_synthetic_with_ca,
				'Wostbrock et al. (2020) synthetic calcite',
			),
			(
				synthetic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				wostbrock_2020.gca_synthetic_without_ca,
				'',
			),
			(
				marine_cc_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				wostbrock_2020.gca_marine_calcite,
				'Wostbrock et al. (2020) marine calcite',
			),
			(
				marine_ar_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				wostbrock_2020.gca_marine_aragonite,
				'Wostbrock et al. (2020) marine aragonite',
			),
		):

			# water data
			dw = np.array([r['w_d18O'] for r in _data])
			Dw = np.array([r['w_D17O'] for r in _data])

			# carb data
			dc = np.array([r['c_d18O'] for r in _data])
			Dc = np.array([r['c_D17O'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = wostbrock_2020.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = wostbrock_2020.D17O_VSMOW_OF_IAEA603,
					D17O_NBS18_old =   wostbrock_2020.D17O_VSMOW_OF_NBS18,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)
			if 'Devils' in label and renormalize:
				print(f"Devils Hole (Wostbrock et al., 2020):\n  water Δ17O_VSMOW = {Dw[0]:.3f} ‰\n  calcite Δ17O_IAEA603 = {Dc[0]-D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY:.3f} ‰")
				with open(f'output/{"alt_" if alternative_D17O_603_of_NBS18 else ""}DVH_comparison.csv', 'w') as f:
					f.write('quantity,W20,SE_W20,this_study,SE_this_study,offset,SE_offset')
					f.write(f'\nD17O_VSMOW_of_water,{Dw[0]:.3f},0.003')
					f.write(f',{data.loc["DVH", "D17O_VSMOW"]:.3f},{data.loc["DVH", "95CL_D17O_VSMOW"]/1.96:.3f}')
					f.write(f',{data.loc["DVH", "D17O_VSMOW"]-Dw[0]:.3f},{((data.loc["DVH", "95CL_D17O_VSMOW"]/1.96)**2 + 0.003**2)**.5:.3f}')
					f.write(f'\nD17O_603_of_carb,{Dc[0] - D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY:.3f},0.007')
					f.write(f',{data.loc["DVH", "D17O_603"]:.3f},{data.loc["DVH", "95CL_D17O_603"]/1.96:.3f}')
					f.write(f',{data.loc["DVH", "D17O_603"] - Dc[0] + D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY:.3f},{((data.loc["DVH", "95CL_D17O_603"]/1.96)**2 + 0.007**2)**.5:.3f}')

				DVH_comparison = pd.read_csv(f'output/{"alt_" if alternative_D17O_603_of_NBS18 else ""}DVH_comparison.csv', index_col = 'quantity')
				config.df_to_toml(DVH_comparison, 'output/labs-results.toml', 'a', f'{"alt_" if alternative_D17O_603_of_NBS18 else ""}DVH_comparison')

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = wostbrock_2020.D17O_VSMOW_OF_IAEA603


			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ------------------------------------------ Plot Huth et al. (2022)
		col = np.array(config.COLORS['Huth'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				synthetic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				huth_2022.data,
				'Huth et al. (2022) synthetic',
			),
		):
			# water data
			dw = np.array([r['w_d18O'] for r in _data])
			Dw = np.array([r['w_D17O'] for r in _data])

			# carb data
			dc = np.array([r['c_d18O'] for r in _data])
			Dc = np.array([r['c_D17O'] for r in _data])

			if renormalize:

				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = huth_2022.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = huth_2022.D17O_VSMOW_OF_IAEA603,
					D17O_NBS18_old = huth_2022.TENTATIVE_D17O_VSMOW_OF_NBS18,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = huth_2022.D17O_VSMOW_OF_IAEA603

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ---------------------------------------------------------------- Plot Kelson et al. (2022)
		col = np.array(config.COLORS['Kelson'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				biogenic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				kelson_2022.data,
				'Kelson et al. (2022) bivalves',
			),
		):
			# water data
			dw = np.array([r['w_d18O'] for r in _data])
			Dw = np.array([r['w_D17O'] for r in _data])
			# carb data
			dc = np.array([r['c_d18O'] for r in _data])
			Dc = np.array([r['c_D17O'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = kelson_2022.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = kelson_2022.D17O_VSMOW_OF_IAEA603,
					D17O_NBS18_old = huth_2022.TENTATIVE_D17O_VSMOW_OF_NBS18,
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = kelson_2022.D17O_VSMOW_OF_IAEA603

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		# ---------------------------------------------------------------- Plot Roy et al. (2025)
		col = np.array(config.COLORS['Roy'])
		for symbol, mec, mfc, ec, _data, label in (
			(
				synthetic_marker,
				1-(1-col)*1,
				(1,1,1,0),
				1-(1-col)*0.33,
				roy_2025.data,
				'Roy et al. (2025) synthetic',
			),
		):
			# water data
			dw = np.array([
				(r['starting_water_d18O_VSMOW'] + r['final_water_d18O_VSMOW'])/2
				for r in _data
			])
			Dw = np.array([
				(r['starting_water_D17O_VSMOW'] + r['final_water_D17O_VSMOW'])/2
				for r in _data
			])

			# carb data
			dc = np.array([r['CO2acid_d18O_VSMOW'] for r in _data])
			dc = (1000 + dc) / 1.01025 - 1000
			Dc = np.array([r['CO2acid_D17O_VSMOW'] for r in _data])

			if renormalize:
				# adjust D17O
				Dw = config.D17O_adjust_water(
					d18O = dw,
					D17O_old = Dw,
					D17O_SLAP_old = roy_2025.D17O_VSMOW_OF_SLAP,
					D17O_SLAP_new = config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
				)

				# adjust D17O
				Dc = config.D17O_adjust_carb(
					d18O = dc,
					D17O_old = Dc,
					D17O_IAEA603_old = (
						roy_2025.D17O_VSMOW_OF_NBS19_CO2
						- config.SLAP_adjust('NBS19-CO2', 'Chaillot et al. (2025)')
						+ config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
					),
					D17O_NBS18_old =   (
						roy_2025.D17O_VSMOW_OF_NBS19_CO2
						- config.SLAP_adjust('NBS19-CO2', 'Chaillot et al. (2025)')
						+ config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')
					),
					D17O_IAEA603_new = (
						config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
					),
					D17O_NBS18_new = (
						config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)')
						if alternative_D17O_603_of_NBS18
						else D17O_VSMOW_OF_NBS18_FIT_FROM_THEORY
					),
				)

			if renormalize:
				if alternative_D17O_603_of_NBS18:
					_D603 = config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)')
				else:
					_D603 = D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
			else:
				_D603 = (
					roy_2025.D17O_VSMOW_OF_NBS19_CO2
					- config.SLAP_adjust('NBS19-CO2', 'Chaillot et al. (2025)')
					+ config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
				)

			T = np.array([r['T'] for r in _data]) + 273.15
			ax.plot(
				T**-2,
				Dc - Dw - _D603,
				ls = 'None',
				marker = symbol[0],
				ms = symbol[1],
				mec = mec,
				mew = mew,
				mfc = mfc,
				zorder = 100,
				label = label,
				)

			reprocessed[label] = dict(
				symbol = symbol,
				mec = mec,
				mew = mew,
				mfc = mfc,
				dw = dw,
				Dw = Dw,
				dc = (1000 + dc) / 1.03092 - 1000,
				Dc = Dc,
			)

		ax.legend(
			loc = 'lower center',
			bbox_to_anchor=(0.535, 0),
			bbox_transform = fig.transFigure,
			labelspacing = 0.18,
			handlelength = 2.75,
			fontsize = 7.5,
			frameon = False
		)

		xlimits = np.array([(tmax + 273.15)**-2, (tmin + 273.15)**-2])
		ylimits = np.array([-0.08, 0.120])
		ax.axis([*xlimits, *ylimits])
		twinax.axis([*xlimits, *ylimits])


		fig.savefig(f'output/interlab-calibrations{"-alt" if alternative_D17O_603_of_NBS18 else "" if renormalize else "-original"}.pdf')
		ppl.close(fig)


R = reprocessed

fig = ppl.figure(figsize = (config.TWO_COLUMN_WIDTH, config.TWO_COLUMN_WIDTH*0.3))
fig.subplots_adjust(.01, .22, .67, .94, .05, .5)

ax_w = ppl.subplot(121)
ax_c = ppl.subplot(122)

ax_w.set_xlabel(r'$δ^{18}O_\text{VSMOW}$ of waters (‰)')
ax_c.set_xlabel(r'$δ^{18}O_\text{VPDB}$ of carbonates (‰)')



y = 0
yt = []
for l,r in R.items():
	yt.append(l)
	y += 1
	ax_w.plot(r['dw'], r['dw']*0+y, ls = 'None', marker = r['symbol'][0], ms = r['symbol'][1], mec = r['mec'], mfc = r['mfc'])
ax_w.axis([None, None, y+1, 0])
ax_w.set_yticks(np.arange(len(yt))+1)
ax_w.set_yticklabels(['' for _ in yt])
ax_w.yaxis.tick_right()
ax_w.tick_params(axis = 'y', length = 0)

for x, name in (
	(0, 'VSMOW  '),
	(-55.5, '  SLAP'),
):
	ax_w.axvline(x, color = 'k', lw = 0.8)
	ax_w.text(x, -0.2, name, ha = 'center', size = 8)
ax_w.xaxis.set_major_locator(MultipleLocator(10))


y = 0
for l,r in R.items():
	print(l)
	y += 1
	ax_c.plot(r['dc'], r['dc']*0+y, ls = 'None', marker = r['symbol'][0], ms = r['symbol'][1], mec = r['mec'], mfc = r['mfc'])
ax_c.axis([None, None, y+1, 0])
ax_c.set_yticks(np.arange(len(yt))+1)
ax_c.set_yticklabels(yt, size = 8)
ax_c.yaxis.tick_right()
ax_c.tick_params(axis = 'y', length = 0)

for x, name in (
	(-2.20, 'IAEA-603'),
	(-23.01, 'NBS18'),
):
	ax_c.axvline(x, color = 'k', lw = 0.8)
	ax_c.text(x, -0.2, name, ha = 'center', size = 8)
ax_c.xaxis.set_major_locator(MultipleLocator(10))


fig.savefig('output/all_unknown_vs_anchors.pdf')


fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, config.ONE_COLUMN_WIDTH))
fig.subplots_adjust(.3, .3, .95, .95)

ax = ppl.subplot(111)

ax.set_xlabel(r'$Δ^{17}O_\text{VSMOW}$ of Devils Hole water (‰)')
ax.set_ylabel(r'$Δ^{17}O_\text{VSMOW}$ of Devils Hole calcite (‰)')
ax.xaxis.set_major_locator(MultipleLocator(0.01))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
ax.axis([5e-3, 55e-3, -55e-3, -05e-3])

X = DVH_comparison.loc['D17O_VSMOW_of_water']['W20']
eX = DVH_comparison.loc['D17O_VSMOW_of_water']['SE_W20']*1.96
Y = DVH_comparison.loc['D17O_603_of_carb']['W20'] + D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
eY = DVH_comparison.loc['D17O_603_of_carb']['SE_W20']*1.96
ax.errorbar(
	X,
	Y,
	eY,
	eX,
	**(config.kw_errorbar | dict(
		marker = 'None',
		ecolor = config.COLORS['Wostbrock'],
	)),
)
ax.plot(
	X,
	Y,
	'o',
	mfc = [1-(1-c)/3 for c in config.COLORS['Wostbrock']],
	mec = config.COLORS['Wostbrock'],
	mew = 1,
	label = 'Wostbrock et al. (2020)',
)

X = DVH_comparison.loc['D17O_VSMOW_of_water']['this_study']
eX = DVH_comparison.loc['D17O_VSMOW_of_water']['SE_this_study']*1.96
Y = DVH_comparison.loc['D17O_603_of_carb']['this_study'] + D17O_VSMOW_OF_IAEA603_FIT_FROM_THEORY
eY = DVH_comparison.loc['D17O_603_of_carb']['SE_this_study']*1.96
ax.errorbar(
	X,
	Y,
	eY,
	eX,
	**(config.kw_errorbar | dict(
		marker = 'None',
		ecolor = 'k',
	)),
)
ax.plot(
	X,
	Y,
	'o',
	mfc = '.75',
	mec = 'k',
	mew = 1,
	label = 'Wostbrock et al. (2020)',
)

fig.savefig(f'output/DVH_comparison{"-alt" if alternative_D17O_603_of_NBS18 else ""}.pdf')

vprint(f'Finished {Path(__file__).name}')
