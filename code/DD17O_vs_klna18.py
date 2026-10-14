#! /usr/bin/env uv run

import lmfit
import config
import guo_zhou_2019

import numpy as np
import pandas as pd
import uncertainties as uc

from pathlib import Path
from config import vprint
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple

def add_inset(ax, xmin, xmax, ymin, ymax, **kwargs):
	fig = ppl.gcf()

	disp_coords = ax.transAxes.transform([[xmin, ymin], [xmax, ymax]])
	(fx0, fy0), (fx1, fy1) = fig_coords = fig.transFigure.inverted().transform(disp_coords)

	x0, y0 = min(fx0, fx1), min(fy0, fy1)
	width, height = abs(fx1 - fx0), abs(fy1 - fy0)

	return fig.add_axes([x0, y0, width, height], **kwargs)

def clever_rounding(
	a, b, c,
	precision = 3,
	T = np.linspace(573.15**-2, 273.15**-2, 1000)**-0.5,
):
	_c = np.round(c, precision)

	for digits in np.arange(10)-4:
		_b = np.round(b, digits)
		maxdiff = (np.abs(_b-b)/T).max()
		if maxdiff < 10**-precision:
			break

	for digits in np.arange(10)-4:
		_a = np.round(a, digits)
		maxdiff = (np.abs(_a-a)/T**2).max()
		if maxdiff < 10**-precision:
			break

	maxdiff = np.abs(
		_a/T**2 + _b/T + _c
		- a/T**2 - b/T - c
	).max()

	return ((_a, _b, _c), maxdiff)



if __name__ == '__main__':

	vprint(f'Starting {Path(__file__).name}')
	alternative_D17O_603_of_NBS18_ref = (
		'Wostbrock et al. (2020)'
		# 'Ellis & Passey (2023)'
	)

	with open('output/DD17O-results.toml', 'w') as f:
		pass

	for alternative_D17O_603_of_NBS18 in (True, False):
		# ---------------------------------------------------------------- Set up figure
		fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, config.ONE_COLUMN_WIDTH * 2))
		fig.subplots_adjust(.21, .07, .99, .94, .11, .11)

		sqTinv = np.linspace(573.15**-2, 263.15**-2, 101)
		Tinv = sqTinv**0.5
		Tk = sqTinv**-0.5

		for CARBREF in ('NBS18', 'IAEA-603'):
			vprint(f'CARBREF = {CARBREF}{", alternative_D17O_603_of_NBS18" if alternative_D17O_603_of_NBS18 else ""}')

			nbs18ref = (CARBREF == 'NBS18')
			iaea603ref = (CARBREF == 'IAEA-603')

			# ---------------------------------------------------------------- Read data
			data = pd.read_csv('output/calib_results.csv', index_col = 'Sample')
			data['1/T'] = (data['T'] + 273.15)**-1
			data['1/T2'] = (data['T'] + 273.15)**-2
			data['weight'] = 1e-99
			data.loc[['LGB', 'DVH'], 'weight'] = 1

			if alternative_D17O_603_of_NBS18:
				data['D17O_603']	 = config.D17O_adjust_carb(
					d18O = (1000 + data['d18O_VPDB'].to_numpy()) * 1.03092 - 1000,
					D17O_old = data['D17O_603'],
					D17O_IAEA603_old = 0.,
					D17O_NBS18_old = config.D17O_VSMOW_OF_NBS18_CO2_CHAILLOT - config.D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT,
					D17O_IAEA603_new = 0.,
					D17O_NBS18_new = config.SLAP_adjust('NBS18', alternative_D17O_603_of_NBS18_ref) - config.SLAP_adjust('IAEA-603', alternative_D17O_603_of_NBS18_ref),
				)
				data['DD17O'] = data['D17O_603'] - data['D17O_VSMOW']

			if nbs18ref:
				if alternative_D17O_603_of_NBS18:
					alternative_D17O_603_of_NBS18 = (
						config.SLAP_adjust('IAEA-603', alternative_D17O_603_of_NBS18_ref)
						- config.SLAP_adjust('NBS18', alternative_D17O_603_of_NBS18_ref)
					)
					data['DD17O'] += alternative_D17O_603_of_NBS18
				else:
					data['DD17O'] += config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
					data['DD17O'] -= config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')
			print(alternative_D17O_603_of_NBS18)

			# ---------------------------------------------------------------- Regression --> D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY
			T = data['T'].to_numpy() + 273.15
			DD17O_obs = data['DD17O'].to_numpy()
			DD17O_obs_sigma = data['95CL_DD17O'].to_numpy() / 1.96
			DD17O_theory = np.array([
				(guo_zhou_2019.calcite_water_eq_theta(_) - 0.528)
				* guo_zhou_2019.calcite_water_eq_klna18(_)
				for _ in T
			])
			klna18_obs = data['klna18'].to_numpy()
			klna18_obs_sigma = data['95CL_klna18'].to_numpy() / 1.96

			weights = data['weight'].to_numpy()

			def _residuals(p):
				return (DD17O_obs - DD17O_theory - p['offset']) / DD17O_obs_sigma * weights

			_p = lmfit.Parameters()
			_p.add('offset', value = 0.07)
			fitresult = lmfit.minimize(_residuals, _p, method = 'least_squares', scale_covar = False)
			with open('output/DD17O_fit_report.txt', 'w') as f:
				f.write(lmfit.fit_report(fitresult))
			D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY = -fitresult.uvars['offset']
			DD17Oeq_afo_T = lambda t: (
				(guo_zhou_2019.calcite_water_eq_theta(t) - 0.528)
				* guo_zhou_2019.calcite_water_eq_klna18(t)
				- D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY
			)

			# ---------------------------------------------------------------- Quadratic approximation of DD17Oeq_afo_T()
			_T = np.linspace((250+273.15)**-4, (273.15)**-4, 100)**-0.25 # higher density at low T
			X = _T**-1
			Y = np.array([_.n for _ in DD17Oeq_afo_T(_T)])
			quadfit = np.polyfit(X, Y, 2)
			quadfit, maxdiff = clever_rounding(*quadfit)


			# ---------------------------------------------------------------- Save results
			altstr = 'alt' if alternative_D17O_603_of_NBS18 else ''
			with open('output/DD17O-results.toml', 'a') as f:
				f.write(f'\n[{altstr}bestfit.D17O_VSMOW_OF_{CARBREF.replace("-", "")}]\n')
				f.write(f'value = {D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.n:.6f}\n')
				f.write(f'se = {D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.s:.6f}\n')
				f.write(f'c95 = {D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.s*1.96:.6f}\n')
				f.write(f'\n[{altstr}bestfit.D17O_VSMOW_OF_{CARBREF.replace("-", "")}.coefs]\n')
				f.write(f'a2 = {quadfit[0]}\n')
				f.write(f'a1 = {quadfit[1]}\n')
				f.write(f'a0 = {quadfit[2]}\n')
				f.write(f'maxdiff = {maxdiff}\n\n')

			# ---------------------------------------------------------------- Set up axes
			ax2 = ppl.subplot(211 if iaea603ref else 212)
			ppl.ylabel(
				(
					# f'$Δ^{{17}}O\\,^\\text{{carb/{CARBREF}}}_\\text{{water/VSMOW}}$'
					# f'$Δ^{{17}}O\\,^\\text{{carb/{CARBREF}}}_\\text{{water/VSMOW}} = Δ^{{17}}O\\,^\\text{{c}}_\\text{{w}} - Δ^{{17}}O\\,_{{VSMOW}}^\\text{{{CARBREF}}}$'
					f'$Δ^{{17}}O\\,^\\text{{c}}_\\text{{w}} - Δ^{{17}}O\\,_{{VSMOW}}^\\text{{{CARBREF}}}$'
					# f' = $Δ^{{17}}O\\,^\\text{{carb}}_\\text{{{CARBREF}}} - Δ^{{17}}O\\,^\\text{{water}}_{{VSMOW}}$'
					'   [‰]'
				),
				size = 8,
			)
			if nbs18ref:
				ax2.set_xlabel(
					'$1000\\cdot$ln($^{18}α^{c}_{w}$)',
					# labelpad = 20,
				)
			ax = ax2.twiny()

			ax_limits = np.array([0, guo_zhou_2019.calcite_water_eq_klna18(273.15), -0.065, 0.145])
			if nbs18ref:
				ax_limits[-2:] += config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
				ax_limits[-2:] -= config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')
			# if alternative_D17O_603_of_NBS18:
			# 	ax_limits[-2:] += -0.02 if nbs18ref else 0.02

			ax.axis(ax_limits)
			Tl = np.array([0, 20, 50, 100, 250])[::-1]
			ax2.axis(ax_limits)
			ax.xaxis.set_ticks(guo_zhou_2019.calcite_water_eq_klna18(Tl+273.15))
			if iaea603ref:
				ax.xaxis.set_ticklabels(
					[f'{t:.0f}$\\,$°C' for t in Tl],
					rotation = -30,
					ha = 'right',
					va = 'bottom',
					size = 9,
					linespacing = 1,
				)
			else:
				ax.xaxis.set_ticklabels(
					['' for t in Tl],
				)
				ax.tick_params(axis = 'x', length = 0)

			ax2.grid(False, axis = 'x')

			ax.yaxis.set_major_locator(MultipleLocator(0.04))
			ax2.xaxis.set_major_locator(MultipleLocator(10))

			# ---------------------------------------------------------------- Plot D17O_VSMOW_OF_CARBREF on yaxis

			blank, = ax.plot([], [], ls = 'None', marker = 'None')
			new_estimate_color = [1-(1-_)/2 for _ in config.COLORS['diseq']]
			bfD, = ax.plot(0, -D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.n,
				'D',
				ms = 5,
				mew = 1,
				mec = 'k',
				mfc = 'w',
				clip_on = False,
				zorder = 1000,
				# label = 'This study (best fit)'
				label = f'$-Δ^{{17}}O^\\text{{{CARBREF}}}_\\text{{VSMOW}}$ (this study)',
			)

			_t = Tk[0]
			ax.add_artist(
				FancyArrowPatch(
					posA = [0, -D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.n],
					posB = [guo_zhou_2019.calcite_water_eq_klna18(_t), DD17Oeq_afo_T(_t).n],
					arrowstyle = '<-',
					shrinkA = 4,
					shrinkB = 0,
					connectionstyle = "arc3,rad=-0.08",
					mutation_scale = 9,
					lw = 1,
					color = 'k',
					clip_on = False,
				)
			)


			epcolor = config.COLORS['Ellis']
			bfD_ep, = ax.plot(0, -config.SLAP_adjust(CARBREF, 'Ellis & Passey (2023)'),
				'D',
				ms = 5,
				mew = 1,
				mec = 'k',
				mfc = 1-(1-epcolor)*2/3,
				clip_on = False,
				zorder = 1000,
				# label = 'Ellis & Passey (2023)',
				label = f'$-Δ^{{17}}O^\\text{{{CARBREF}}}_\\text{{VSMOW}}$ (Ellis & Passey, 2023)',
			)

			wcolor = config.COLORS['Wostbrock']
			bfD_wb, = ax.plot(0, -config.SLAP_adjust(CARBREF, 'Wostbrock et al. (2020)'),
				'D',
				ms = 5,
				mew = 1,
				mec = 'k',
				mfc = 1-(1-wcolor)*2/3,
				clip_on = False,
				zorder = 1000,
				# label = 'Wostbrock et al. (2020a)'
				label = f'$-Δ^{{17}}O^\\text{{{CARBREF}}}_\\text{{VSMOW}}$ (Wostbrock et al., 2020a)',
			)

			# ---------------------------------------------------------------- Plot theoretical model wih best-fit offset
			DD = DD17Oeq_afo_T(Tk)
			DDn = np.array([_.n for _ in DD])
			DDe = np.array([_.s for _ in DD]) * 1.96
			# gze = ax.fill_between(
			# 	guo_zhou_2019.calcite_water_eq_klna18(Tk),
			# 	DDn - DDe,
			# 	DDn + DDe,
			# 	lw = 0,
			# 	color = config.COLORS['diseq'],
			# 	alpha = 0.2,
			# 	zorder = 0,
			# )

			gzt, = ax.plot(
				guo_zhou_2019.calcite_water_eq_klna18(Tk),
				DDn,
				'k-',
				lw = 1,
				label = 'Best fit Guo & Zhou (2019) trend',
				solid_capstyle = 'round',
				)

			# a,b,c = quadfit
			# Tkk = np.linspace(1e-20, 273.15**-2)**-0.5
			# ax.plot(
			# 	guo_zhou_2019.calcite_water_eq_klna18(Tkk),
			# 	a/Tkk**2 + b/Tkk + c,
			# 	'r:',
			# 	lw = 1,
			# 	zorder = 10000,
			# 	label = 'Quadratic approximation of Guo & Zhou (2019) trend',
			# 	)

			ax.scatter(
				guo_zhou_2019.calcite_water_eq_klna18(data['T'] + 273.15),
				[_.n for _ in DD17Oeq_afo_T(data['T'] + 273.15)],
				marker = config.kw_errorbar['marker'],
				# cmap = cmap,
				# edgecolors = 'k',
				# c = Xn,
				c = 'k',
				s = 3**2,
				linewidths = 0,
				zorder = 201,
			)

			# ---------------------------------------------------------------- Plot disequilibrium arrows
			for t, x, y in zip(
				data['T'],
				data['klna18'],
				data['DD17O'],
			):
				if abs(y - DD17Oeq_afo_T(t + 273.15).n) > 0.01:
					ax.add_artist(
						FancyArrowPatch(
							posA = [
									guo_zhou_2019.calcite_water_eq_klna18(t + 273.15),
									DD17Oeq_afo_T(t + 273.15).n
								],
							posB = (x, y),
							arrowstyle = '->',
							shrinkA = 0,
							shrinkB = 4,
							mutation_scale = 8,
							lw = 0.8,
							color = 1-(1-config.COLORS['diseq'])*2/3,
						)
					)
					# ax.annotate("",
					# 	xytext = [
					# 		guo_zhou_2019.calcite_water_eq_klna18(t + 273.15),
					# 		DD17Oeq_afo_T(t + 273.15).n
					# 	],
					# 	xy = [x, y],
					# 	arrowprops = {"arrowstyle": "->"},
					# )

			# ---------------------------------------------------------------- Quantify isotopic disequilibrium
			Xdiseq = [
				uc.ufloat(n, s / 1.96) - guo_zhou_2019.calcite_water_eq_klna18(t + 273.15)
				for n,s,t in zip(
					data['klna18'],
					data['95CL_klna18'],
					data['T'],
				)
			]
			Ydiseq = [
				uc.ufloat(n, s / 1.96) - DD17Oeq_afo_T(t + 273.15)
				for n, s, t in zip(
					data['DD17O'],
					data['95CL_DD17O'],
					data['T'],
				)
			]
			Xn = np.array([x.n for x in Xdiseq])
			Xs = np.array([x.s for x in Xdiseq])
			Yn = np.array([y.n for y in Ydiseq])
			Ys = np.array([y.s for y in Ydiseq])
			# Xnrank = Xn.argsort().argsort()
			# data['Xnrank'] = Xnrank
			data['klna18_diseq_nominal'] = Xn
			data['klna18_diseq_sigma'] = Xs
			data['D17O_diseq_nominal'] = Yn
			data['D17O_diseq_sigma'] = Ys

			# ---------------------------------------------------------------- Plot our data
			ax.errorbar(
				data['klna18'],
				data['DD17O'],
				data['95CL_DD17O'],
				# data['95CL_klna18'],
				ecolor = 'k',
				zorder = 100,
				**(config.kw_errorbar | dict(marker = 'None')),
				)

			diseq_color = np.array(config.COLORS['diseq'])
			cmap = LinearSegmentedColormap.from_list(
				'ax3ccmap',
				['w', diseq_color]
			)

			thisstudy = {}
			for _data, _marker, _label, _ms in (
				(data[data['T'] < 50], 'o', 'natural calcites', config.kw_errorbar['ms']*1.1),
				(data[data['T'] > 50], 's', 'synthetic calcites', config.kw_errorbar['ms']),
			):
				thisstudy[_marker], = ax.plot(
					[], [], 'w'+_marker,
					ms = _ms,
					mec=  'k',
					label = _label,
				)
				ax.scatter(
					_data['klna18'],
					_data['DD17O'],
					marker = _marker,
					cmap = cmap,
					vmin = data['klna18_diseq_nominal'].min(),
					vmax = data['klna18_diseq_nominal'].max(),
					edgecolors = 'k',
					c = _data['klna18_diseq_nominal'],
					s = _ms**2,
					zorder = 201,
				)

			ax.legend(
				[
					thisstudy['o'],
					thisstudy['s'],
					blank,
					gzt,
					blank,
					bfD,
					bfD_ep,
					bfD_wb,
				],
				[
					thisstudy['o'].get_label(),
					thisstudy['s'].get_label(),
					# f'\n$\\hspace{{-2.5}}\\text{{Estimates of }}-Δ^{{17}}O^\\text{{{CARBREF}}}_\\text{{VSMOW}}$',
					'',
					gzt.get_label(),
					'',
					bfD.get_label(),
					bfD_ep.get_label(),
					bfD_wb.get_label(),
				],
				loc = 'lower left',
				fontsize = 7.5,
				handler_map = {tuple: HandlerTuple(ndivide = 1)},
				frameon = False,
				handlelength = 1.3,
				labelspacing = 0.25,
				borderpad = 0,
				bbox_to_anchor = (0.02,0),
			)

			if not alternative_D17O_603_of_NBS18:
				# ---------------------------------------------------------------- Inset plot
				ax3color = diseq_color * 0.75
				ax3 = add_inset(ax, 0.66, 0.975, 0.68, 0.975)
				for spine in ax3.spines.values():
					spine.set_visible(False)
					# spine.set_color(diseq_color)
				ax3.add_patch(
					FancyBboxPatch(
						(0, 0), 1, 1,
						transform = ax3.transAxes,
						boxstyle = "round,pad=0,rounding_size=0.05",
						linewidth = 0.8,
						edgecolor = ax3color,
						facecolor = 'w',
						clip_on = False,
						zorder = 10,
					)
				)

				ax3.xaxis.set_major_locator(MultipleLocator(1))
				ax3.yaxis.set_major_locator(MultipleLocator(20e-3))
				ax3.set_xlabel('$δ^{18}O$ diseq. [‰]', size = 7.5, color = ax3color)
				ax3.set_ylabel('$Δ^{17}O$ diseq. [‰]', size = 7.5, color = ax3color)
				ax3.tick_params(axis = 'both', which = 'major', labelsize = 8, length = 2, width = 0.8, colors = ax3color, zorder = 10)
				ax3.axhline(0, color = ax3color, lw = 0.8, zorder = 15, dashes = (3,2))
				ax3.axvline(0, color = ax3color, lw = 0.8, zorder = 15, dashes = (3,2))

				ax3.errorbar(
					Xn,
					Yn,
					1.96*Ys,
					1.96*Xs,
					ls = 'None',
					marker = 'None',
					ecolor = 1-(1-ax3color)/3,
					capthick = 0.8,
					elinewidth = 0.8,
					capsize = 2,
					zorder = 20,
				)

				for _data, _marker in (
					(data[data['T'] < 50], 'o'),
					(data[data['T'] > 50], 's'),
				):
					ax3.scatter(
						_data['klna18_diseq_nominal'],
						_data['D17O_diseq_nominal'],
						marker = _marker,
						cmap = cmap,
						vmin = data['klna18_diseq_nominal'].min(),
						vmax = data['klna18_diseq_nominal'].max(),
						edgecolors = ax3color,
						linewidths = 0.8,
						c = _data['klna18_diseq_nominal'],
						s = 4**2,
						zorder = 21,
					)

			ax.text(
				0.03, 0.975,
				{'IAEA-603': 'a', 'NBS18': 'b'}[CARBREF],
				weight = 'bold',
				size = 12,
				va = 'top',
				ha = 'left',
				transform = ax.transAxes,
			)

			if alternative_D17O_603_of_NBS18:
				ax.text(
					0.71, 0.97,
					 f' alternative assumption: \n $Δ^{{17}}O^\\text{{NBS18}}_\\text{{IAEA-603}} = {alternative_D17O_603_of_NBS18*1000:.0f}$ ppm ',
					va = 'top',
					ha = 'center',
					size = 9,
					transform = ax.transAxes,
					linespacing = 1.6,
					bbox=dict(
						boxstyle = "round, pad = 0.1, rounding_size = 0.3",
						fc = '0.95',
						ec = 'k',
						lw = 0.7,
					)
				)

		# ---------------------------------------------------------------- Save results
		fig.savefig(
			'output/DD17O_vs_klna18_with_W20_offset.pdf'
			if alternative_D17O_603_of_NBS18 else
			'output/DD17O_vs_klna18.pdf'
		)
		ppl.close(fig)

		altstr = 'alt' if alternative_D17O_603_of_NBS18 else ''
		if iaea603ref:
			config.df_to_toml(data, 'output/DD17O-results.toml', 'a', f'{altstr}data')

	vprint(f'Finished {Path(__file__).name}')
