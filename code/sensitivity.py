#! /usr/bin/env uv run

import lmfit
# import tomllib
import config
import hayles_2018
import guo_zhou_2019
import petroff_2021

import numpy as np
import pandas as pd
import uncertainties as uc

from pathlib import Path
from config import vprint
from transLegend import transLegend, get_transform_aspect
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import FancyBboxPatch
# from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple


def fit_theta(
	data,
	POSTULATED_D17O_VSMOW_OF_SLAP    =  config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
	POSTULATED_D17O_603_OF_NBS18     =  (
		config.D17O_VSMOW_OF_NBS18_CO2_CHAILLOT
		- config.D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT
	),
	POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.0874 # Ellis & Passey (2023) adjusted for D17O_VSMOW_OF_SLAP = 0
):

	data = data.copy()

	# Renormalize waters to postulated SLAP composition
	data['D17O_VSMOW'] = (
		data['D17O_VSMOW']
		+ (POSTULATED_D17O_VSMOW_OF_SLAP - config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE)
		* np.log(1 + data['d18O_VSMOW']/1e3)
		/ np.log(1 - 55.5/1e3)
	)

	# Renormalize carbonates to postulated NBS18 composition
	data['D17O_603'] = (
		data['D17O_603']
		+ (POSTULATED_D17O_603_OF_NBS18 - config.D17O_603_OF_NBS18_POSTULATED_HERE)
		* (np.log(1 + data['d18O_VPDB']/1e3) - np.log(1 - 2.37/1e3))
		/ (np.log(1 - 23.01/1e3) - np.log(1 - 2.37/1e3))
	)

	# Recompute DD17O based on renormalized values
	data['DD17O'] = data['D17O_603'] - data['D17O_VSMOW']
	data['95CL_DD17O'] = (data['95CL_D17O_603']**2 + data['95CL_D17O_VSMOW']**2)**.5

	out = {}

	# # Fit Δ17O_603_cc - Δ17O_VSMOW_water as a quadratic function of 1/T
	# X = np.array(data['1/T'])
	# Y = np.array(data['DD17O'])
	# eY = np.array(data['95CL_DD17O'] / data['weight'])
	# A = np.array([[1/e, x/e, (x**2)/e] for x,e in zip(X, eY)])
	# F, chisq = np.linalg.lstsq(A, Y/eY, rcond = None)[:2]
	# print(chisq, (chisq / (Y.size - 3))**.5)
	# CM = np.linalg.inv(A.T @ A) # * chisq / (Y.size - 3)
	# out = {
	# 	'a0': F[0],
	# 	'a1': F[1],
	# 	'a2': F[2],
	# 	'95CL_a0': CM[0,0]**0.5,
	# 	'95CL_a1': CM[1,1]**0.5,
	# 	'95CL_a2': CM[2,2]**0.5,
	# 	'correl_a0_a1': CM[0,1] / CM[0,0]**0.5 / CM[1,1]**0.5,
	# 	'correl_a0_a2': CM[0,2] / CM[0,0]**0.5 / CM[2,2]**0.5,
	# 	'correl_a1_a2': CM[1,2] / CM[1,1]**0.5 / CM[2,2]**0.5,
	# }

	# Prepare data to fit theta
	data['theta'] = (
		data['DD17O']
		+ POSTULATED_D17O_VSMOW_OF_IAEA603
	) / data['klna18'] + 0.528
	data['95CL_theta'] = data['95CL_DD17O'] / data['klna18']
	X = np.array(data['1/T'])
	Y = np.array(data['theta'])
	eY = np.array(
		data['95CL_theta']
		/ data['weight']
	)

	# Fit theta as an affine function of 1/T
	A = np.array([[1/e, x/e] for x,e in zip(X, eY)])
	F, chisq = np.linalg.lstsq(A, Y/eY, rcond = None)[:2]
	CM = np.linalg.inv(A.T @ A) * chisq / (Y.size - 2)
	out = out | {
		'b0': F[0],
		'b1': F[1],
		'95CL_b0': CM[0,0]**0.5,
		'95CL_b1': CM[1,1]**0.5,
		'correl_b0_b1': CM[0,1] / CM[0,0]**0.5 / CM[1,1]**0.5,
	}

	# Fit theta as an affine function of 1/T going through HIGH_T_LIMIT_OF_THETA
	A = np.array([[x/e] for x,e in zip(X, eY)])
	F, chisq = np.linalg.lstsq(A, (Y - config.HIGH_T_LIMIT_OF_THETA)/eY, rcond = None)[:2]
	CM = np.linalg.inv(A.T @ A) * chisq / (Y.size - 1)
	out = out | {
		'c0': config.HIGH_T_LIMIT_OF_THETA,
		'c1': F[0],
		'95CL_c1': CM[0,0]**0.5,
	}

	out['data'] = data.copy()
	return out

if __name__ == '__main__':

	vprint(f'Starting {Path(__file__).name}')

	# ---------------------------------------------------------------- Read data
	CULL_EQUILIBRIUM_SAMPLES = True

	data = pd.read_csv('output/calib_results.csv', index_col = 'Sample')
	data['1/T'] = (data['T'] + 273.15)**-1
	data['1/T2'] = (data['T'] + 273.15)**-2
	Xdiseq = [
		uc.ufloat(n, s / 1.96) - guo_zhou_2019.calcite_water_eq_klna18(t + 273.15)
		for n,s,t in zip(
			data['klna18'],
			data['95CL_klna18'],
			data['T'],
		)
	]
	Xn = np.array([x.n for x in Xdiseq])
	Xs = np.array([x.s for x in Xdiseq])
	data['weight'] = abs(Xs/Xn)
	if CULL_EQUILIBRIUM_SAMPLES:
		data = data.loc[['LGB', 'DVH', 'CC170', 'CC200']]
		data.weight = 1.
	print(data)


	# ---------------------------------------------------------------- Set up figure and subplots

	fig, _ = ppl.subplot_mosaic(
		'a\nb',
		figsize = (config.ONE_COLUMN_WIDTH, config.ONE_COLUMN_WIDTH * 1.8),
	)
	ax_theta, ax_DD17O = _['a'], _['b']
	fig.subplots_adjust(.225, .08, .985, .95, 0.4, 0.9)

	ax_theta.set_xlabel('$1\\,/\\,T$', labelpad = 5)
	ax_theta.set_ylabel('$^{17}θ^c_w$')
	ax_DD17O.set_xlabel(r'$10^3\cdot ln(^{18}α^\text{eq}_{c/w})$')
	ax_DD17O.set_ylabel(r'$Δ^{17}O^\text{carbonate}_\text{IAEA-603} - Δ^{17}O_\text{VSMOW}^\text{water}$ [‰]')

	Ti = np.array([0, 50, 100, 250][::-1])
	ax_theta.xaxis.set_ticks((Ti + 273.15)**-1)
	ax_theta.xaxis.set_ticklabels(
		[f'{t:.0f}$\\,$°C' for t in Ti],
		rotation = -20,
		ha = 'right',
		va = 'bottom',
		size = 9,
		linespacing = 1,
		)
	ax_theta.tick_params(
		axis = 'x',
		pad = 0,
		top = True,
		labeltop = True,
		bottom = False,
		labelbottom = False,
		)

	ax_DD17O_twin = ax_DD17O.twiny()
	Ti = np.array([0, 25, 50, 100, 250][::-1])
	ax_DD17O_twin.xaxis.set_ticks(guo_zhou_2019.calcite_water_eq_klna18(Ti+273.15))
	ax_DD17O_twin.xaxis.set_ticklabels(
		[f'{t:.0f}$\\,$°C' for t in Ti],
		rotation = -20,
		ha = 'right',
		va = 'bottom',
		size = 9,
		linespacing = 1,
		)
	ax_DD17O_twin.tick_params(
		axis = 'x',
		pad = 0,
		top = True,
		labeltop = True,
		bottom = False,
		labelbottom = False,
		)
	ax_DD17O.grid(False, axis = 'x')

	Tinv_DD17O  = np.linspace((300+273.15)**-1, 273.15**-1)
	Tinv_theta = np.linspace((500+273.15)**-1, 273.15**-1)

	ax_theta.axis([Tinv_theta[0], Tinv_theta[-1], 0.5195, 0.5356])
	ax_DD17O.axis([0, guo_zhou_2019.calcite_water_eq_klna18(273.15), -.045, 0.13])
	ax_DD17O_twin.axis(ax_DD17O.axis())

	ax_DD17O.yaxis.set_major_locator(MultipleLocator(4e-2))
	ax_theta.yaxis.set_major_locator(MultipleLocator(4e-3, offset = 2e-3))

	legfirst, = ax_DD17O.plot(
		[], [], ls = 'None', marker = 'None',
		label = (
			"$\\hspace{1.4}\\mathbf{Δ^{17}O_{VSMOW}^{SLAP} "
			+ "\\hspace{3.3} Δ^{17}O_{VSMOW}^{IAEA603} "
			+ "\\hspace{3.1} Δ^{17}O_{IAEA603}^{NBS18}}$"
		),
	)

	leg1, leg2 = {}, {}
	for inputs, color, dashes, id in [
		(
			dict(
				POSTULATED_D17O_VSMOW_OF_SLAP    =  0.000,
				POSTULATED_D17O_603_OF_NBS18     =  0.026,
				POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.068,
			), [0, 0, 0], (2.06, 2.06), 'base',
		),
		(
			dict(
				POSTULATED_D17O_VSMOW_OF_SLAP    = -0.015,
				POSTULATED_D17O_603_OF_NBS18     =  0.026,
				POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.068,
			), petroff_2021.six_color_sequence['blue'], (None, None), 'SLAP',
		),
		(
			dict(
			POSTULATED_D17O_VSMOW_OF_SLAP    =  0.000,
			POSTULATED_D17O_603_OF_NBS18     =  0.026,
			POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.100,
			), petroff_2021.six_color_sequence['red'], (None, None), '603',
		),
		(
			dict(
				POSTULATED_D17O_VSMOW_OF_SLAP    =  0.000,
				POSTULATED_D17O_603_OF_NBS18     =  0.050,
				POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.068,
			), petroff_2021.six_color_sequence['orange'], (None, None), 'NBS18',
		),
		# (
		# 	dict(
		# 		POSTULATED_D17O_VSMOW_OF_SLAP    = -0.015,
		# 		POSTULATED_D17O_603_OF_NBS18     =  0.050,
		# 		POSTULATED_D17O_VSMOW_OF_IAEA603 = -0.100,
		# 	), [0,1,0], (None, None), 'NBS18+603+SLAP',
		# ),
	]:
		reg = fit_theta(data, **inputs)
		_color = [1-(1-c)*0.7 for c in color]
		__color = [*color, .25]
		lw = 1

		eY = (
			reg['95CL_b0']**2
			+ (reg['95CL_b1'] * Tinv_theta)**2
			+ 2 * reg['95CL_b0'] * reg['95CL_b1'] * reg['correl_b0_b1'] * Tinv_theta
		)**.5

		ax_theta.fill_between(
			Tinv_theta,
			reg['b0'] + reg['b1'] * Tinv_theta + eY,
			reg['b0'] + reg['b1'] * Tinv_theta - eY,
			color = __color,
			lw = 0,
			zorder = 50,
			)

		ax_theta.errorbar(
			reg['data']['1/T'],
			reg['data']['theta'],
			reg['data']['95CL_theta'],
			ls = 'None',
			marker = 'None',
			ecolor = _color,
			elinewidth = lw,
			capthick = lw,
			capsize = 2,
			zorder = 10,
		)

		ax_theta.plot(
			reg['data']['1/T'],
			reg['data']['theta'],
			'wo',
			mec = _color,
			ms = 4,
			mew = lw,
			zorder = 20,
		)

		ax_DD17O.errorbar(
			reg['data']['klna18'],
			reg['data']['DD17O'],
			reg['data']['95CL_DD17O'],
			# reg['data']['95CL_klna18'],
			ls = 'None',
			marker = 'None',
			ecolor = _color,
			elinewidth = lw,
			capthick = lw,
			capsize = 2,
			zorder = 10,
		)

		ax_DD17O.plot(
			reg['data']['klna18'],
			reg['data']['DD17O'],
			'wo',
			mew = lw,
			mec = _color,
			ms = 4,
			zorder = 20,
		)

		# regression of DD as a quadratic function of T
		T = np.array(data['T']) + 273.15
		DD17O = np.array(reg['data']['DD17O'])
		eDD17O = np.array(data['95CL_DD17O']) / data['weight']

		def _residuals(p):
			return (
				DD17O
				- p['a2']/T**2
				# - p['a1']/T
				- p['a0']
			) / eDD17O

		_p = lmfit.Parameters()
		_p.add('a0', value = 0.080)
		# _p.add('a1', value = 0.)
		_p.add('a2', value = 0.)
		fitresult = lmfit.minimize(_residuals, _p, method = 'least_squares', scale_covar = True)
		a0 = fitresult.uvars['a0']
		# a1 = fitresult.uvars['a1']
		a2 = fitresult.uvars['a2']

		alpha = 0.2 if id == 'base' else 0.25

		leg1[id], = ax_DD17O.plot([], [], '-', color = color, alpha = alpha, lw = 6, solid_capstyle = 'butt')

		X = guo_zhou_2019.calcite_water_eq_klna18(1/Tinv_DD17O)
		Y = a0.n + a2.n * Tinv_DD17O**2

		leg2[id], = ax_DD17O.plot(
			X,
			Y,
			'-',
			color = color,
			lw = 1,
			solid_capstyle = 'butt',
			zorder = 1 if color == [0,0,0] else 0,
			dashes = dashes,
			label = (
				"$\\hspace{1.1}"
				+ f"\\hspace{{0}}{'\\hspace{1.1}' if inputs['POSTULATED_D17O_VSMOW_OF_SLAP'] == 0 else ''}"
				+ f"{inputs['POSTULATED_D17O_VSMOW_OF_SLAP']:.3f}\\ ‰"
				+ f"\\hspace{{4.3}} "
				+ f"{inputs['POSTULATED_D17O_VSMOW_OF_IAEA603']:.3f}\\ ‰"
				+ f"\\hspace{{4.5}} "
				+ f"{inputs['POSTULATED_D17O_603_OF_NBS18']:.3f}\\ ‰$"
			),
		)

		# ax_DD17O.plot(
		# 	guo_zhou_2019.calcite_water_eq_klna18(1/Tinv_theta),
		# 	(
		# 		guo_zhou_2019.calcite_water_eq_klna18(1/Tinv_theta)
		# 		* (reg['b0'] + reg['b1'] * Tinv_theta - 0.528)
		# 		- inputs['POSTULATED_D17O_VSMOW_OF_IAEA603']
		# 	),
		# 	':',
		# 	color = color,
		# 	lw = lw * 1.25,
		# 	zorder = 100 if color == [0,0,0] else 99,
		# )

		# Tinv_theta2 = np.array([0, 1/273.15])
		# ax_theta.plot(
		# 	Tinv_theta2,
		# 	config.HIGH_T_LIMIT_OF_THETA + reg['c1'] * Tinv_theta2,
		# 	'-',
		# 	lw = lw/2,
		# 	color = color,
		# 	# dashes = (6,1.75,2,1.75),
		# 	# label = f'\nAffine model constrained to\na zero-intercept of {config.HIGH_T_LIMIT_OF_THETA:.4f}',
		# 	label = f'Affine fit with high-T limit of {config.HIGH_T_LIMIT_OF_THETA:.4f}',
		# 	)

	ax_theta.plot(
		Tinv_theta[1/Tinv_theta < 573.15],
		guo_zhou_2019.calcite_water_eq_theta(1/Tinv_theta)[1/Tinv_theta < 573.15],
		'k',
		label = 'Theoretical predictions',
		lw = lw*0.7,
		ls = (0, (4,4)),
		)

	ax_theta.plot(
		Tinv_theta[1/Tinv_theta < 573.15],
		hayles_2018.theta_calcite_water_semiempirical(1/Tinv_theta)[1/Tinv_theta < 573.15],
		'k',
		lw = lw*0.7,
		ls = (4, (4,4)),
		)

	ax_theta.legend(
		loc = 'lower left',
		fontsize = 7.4,
		bbox_to_anchor = (0.03, 0.02),
		labelspacing = -0.2,
		frameon = False,
		borderpad = 0,
		handlelength = 2,
	)

	leg = ax_DD17O.legend(
		[legfirst] + [(leg1[_], leg2[_]) for _ in leg1],
		[legfirst.get_label()] + [leg2[_].get_label() for _ in leg1],
		loc = 'center',
		fontsize = 8,
		bbox_to_anchor = (0.5, 0.52),
		labelspacing = 0.4,
		bbox_transform = fig.transFigure,
		handlelength = 1.75,
		borderpad = 1,
		borderaxespad = 0.3,
		handler_map = {tuple: HandlerTuple(ndivide = 1)},
	)

	ax_theta.text(
		.5,
		.98,
		'Regression of $\\mathbf{^{17}θ^c_w}$ = a + b$\\,$/ T',
		size = 9,
		weight = 'bold',
		transform = ax_theta.transAxes,
		va = 'top',
		ha = 'center',
		zorder = 300,
	)

	ax_theta.text(
		.96,
		.04,
		'a',
		size = 12,
		weight = 'bold',
		transform = ax_theta.transAxes,
		va = 'bottom',
		ha = 'right',
		zorder = 300,
	)

	# ax_theta.text(
	# 	.96,
	# 	.97,
	# 	'b',
	# 	size = 11,
	# 	weight = 'bold',
	# 	transform = ax_theta.transAxes,
	# 	va = 'top',
	# 	ha = 'right',
	# 	zorder = 300,
	# )

	ax_DD17O.text(
		0.02,
		0.02,
		(
			'Regression of\n'
			r'$\mathbf{Δ^{17}O^\mathbf{carbonate}_\mathbf{IAEA603} - Δ^{17}O^\mathbf{water}_\mathbf{VSMOW}}$ = a + b$\,$/ T²'
		),
		size = 9,
		weight = 'bold',
		transform = ax_DD17O.transAxes,
		va = 'bottom',
		ha = 'left',
		zorder = 300,
		linespacing = 1.6,
	)

	ax_DD17O.text(
		0.96,
		0.96,
		'b',
		size = 12,
		weight = 'bold',
		transform = ax_DD17O.transAxes,
		va = 'top',
		ha = 'right',
		zorder = 300,
		linespacing = 1.6,
	)

	# ax_DD17O.text(
	# 	.96,
	# 	.97,
	# 	'a',
	# 	size = 11,
	# 	weight = 'bold',
	# 	transform = ax_DD17O.transAxes,
	# 	va = 'top',
	# 	ha = 'right',
	# 	zorder = 300,
	# )

	leg.get_frame().set_facecolor('None')
	_tl = transLegend(leg)
	aspect = get_transform_aspect(_tl)
	for x,y,w,h,c in (
		(.165, .569, .79, .14, 'gray'),
		(.165, .407, .19, .15, 'blue'),
		(.475, .254, .19, .15, 'red'),
		(.765, .103, .19, .15, 'orange'),
	):
		c = petroff_2021.six_color_sequence[c]
		ax_DD17O.add_patch(
			FancyBboxPatch(
				(x, y), w, h,
				boxstyle = "round,pad=0.0,rounding_size=0.01",
				transform = _tl,
				fc = (*c, .5),
				lw = 0,
				zorder = 0,
				clip_on = False,
				mutation_aspect = aspect,
			)
		)

	# ---------------------------------------------------------------- Save results
	fig.savefig('output/sensitivity.pdf')
	ppl.close(fig)

	vprint(f'Finished {Path(__file__).name}')
