#! /usr/bin/env uv run

import lmfit
# import tomllib
import config
import guo_zhou_2019

import numpy as np
import pandas as pd

from scipy import stats
from pathlib import Path
from config import vprint
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator

def carb_D17O_adjust(
	d18O,
	D17O,
	D17O_IAEA603_old,
	D17O_NBS18_old,
	d18O_IAEA603 = (1000-2.37) * 1.03092 - 1000,
	D17O_IAEA603_new = config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)'),
	d18O_NBS18 = (1000-23.01) * 1.03092 - 1000,
	D17O_NBS18_new = config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)'),
):
	x = (
		np.log(1+d18O/1e3)
		- np.log(1+d18O_IAEA603/1e3)
	) / (
		np.log(1+d18O_NBS18/1e3)
		- np.log(1+d18O_IAEA603/1e3)
	)
	return D17O + x * (D17O_NBS18_new - D17O_NBS18_old) + (1-x) * (D17O_IAEA603_new - D17O_IAEA603_old)

if __name__ == '__main__':

	vprint(f'Starting {Path(__file__).name}')

	# ---------------------------------------------------------------- Set up figure
	fig = ppl.figure(figsize = (config.TWO_COLUMN_WIDTH, config.TWO_COLUMN_WIDTH * 1))
	fig.subplots_adjust(.11, .25, .98, .92, .1, .1)

	reprocessed = {}

	for CARBREF in ('IAEA-603', 'NBS18'):

		reprocessed[CARBREF] = {}

		nbs18ref = (CARBREF == 'NBS18')
		iaea603ref = (CARBREF == 'IAEA-603')

		ax_T = ppl.subplot(221 if iaea603ref else 223)
		ax_klna18 = ppl.subplot(222 if iaea603ref else 224, sharey = ax_T)

		ax_T.set_ylabel(
			(
				f'$Δ^{{17}}O\\,^\\text{{c}}_\\text{{w}} - Δ^{{17}}O\\,_{{VSMOW}}^\\text{{{CARBREF}}}$'
				# f' = $Δ^{{17}}O\\,^\\text{{carb}}_\\text{{{CARBREF}}} - Δ^{{17}}O\\,^\\text{{water}}_{{VSMOW}}$'
				'   [‰]'
			),
			size = 8,
		)
		if nbs18ref:
			ax_klna18.set_xlabel(
				'$1000\\cdot$ln($^{18}α^{c}_{w}$)',
				# labelpad = 20,
			)
			ax_T.set_xlabel(
				r'$1\,/\,T^2$',
				# labelpad = 20,
			)

		tmin = -6
		tlimits = np.array([0, (tmin + 273.15)**-2])
		ylimits = np.array([-0.06, 0.15])
		if nbs18ref:
			ylimits -= 0.04
		twinax_T = ax_T.twiny()
		twinax_T.axis([*tlimits, *ylimits])
		ax_T.axis([*tlimits, *ylimits])
		Tl = np.array([0, 20, 50, 100, 250])[::-1]
		twinax_T.xaxis.set_ticks((Tl+273.15)**-2)
		if iaea603ref:
			twinax_T.xaxis.set_ticklabels(
				[f'{t:.0f}$\\,$°C' for t in Tl],
				rotation = -30,
				ha = 'right',
				va = 'bottom',
				size = 9,
				linespacing = 1,
			)
			ax_T.tick_params(axis = 'x', length = 0, labelbottom = False)
		else:
			twinax_T.tick_params(axis = 'x', length = 0, labeltop = False)
		ax_T.grid(False, axis = 'x')

		klnalimits = np.array([0, guo_zhou_2019.calcite_water_eq_klna18(273.15)+2])
		twinax_klna18 = ax_klna18.twiny()
		twinax_klna18.axis([*klnalimits, *ylimits])
		ax_klna18.axis([*klnalimits, *ylimits])
		Tl = np.array([0, 20, 50, 100, 250])[::-1]
		twinax_klna18.xaxis.set_ticks(guo_zhou_2019.calcite_water_eq_klna18(Tl+273.15))
		if iaea603ref:
			twinax_klna18.xaxis.set_ticklabels(
				[f'{t:.0f}$\\,$°C' for t in Tl],
				rotation = -30,
				ha = 'right',
				va = 'bottom',
				size = 9,
				linespacing = 1,
			)
			ax_klna18.tick_params(axis = 'x', length = 0, labelbottom = False)
		else:
			twinax_klna18.tick_params(axis = 'x', length = 0, labeltop = False)
		ax_klna18.tick_params(axis = 'y', length = 0, labelleft = False)

		# ---------------------------------------------------------------- Read our data
		data = pd.read_csv('output/calib_results.csv', index_col = 'Sample')
		data['1/T'] = (data['T'] + 273.15)**-1
		data['1/T2'] = (data['T'] + 273.15)**-2

		if nbs18ref:
			data['DD17O'] += config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
			data['DD17O'] -= config.SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')

		# ---------------------------------------------------------------- Regression --> D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY
		T = data['T'].to_numpy() + 273.15
		DD17O_obs = data['DD17O'].to_numpy()
		DD17O_obs_sigma = data['95CL_DD17O'].to_numpy()
		DD17O_theory = np.array([
			(guo_zhou_2019.calcite_water_eq_theta(_) - 0.528)
			* guo_zhou_2019.calcite_water_eq_klna18(_)
			for _ in T
		])
		klna18_obs = data['klna18_observed'].to_numpy()
		klna18_obs_sigma = data['95CL_klna18_observed'].to_numpy() / 1.96

		if config.ONLY_CONSIDER_NATURAL_SAMPLES_FOR_REGRESSION:
			data['klna18_mismatch'] = 1e99
			data.loc[['LGB', 'DVH'], 'klna18_mismatch'] = 1

		klna18_mismatch = data['klna18_mismatch'].to_numpy()

		def _residuals(p):
			return (DD17O_obs - DD17O_theory - p['offset']) / DD17O_obs_sigma / klna18_mismatch

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

		for ax, twinax in ((ax_T, twinax_T), (ax_klna18, twinax_klna18)):
			isT = (ax == ax_T)

			# ---------------------------------------------------------------- Plot Huth et al. (2022)
			import huth_2022
			col = np.array(config.COLORS['Wostbrock'])
			for symbol, mec, mfc, ec, _data, label in (
				(
					config.ldiamond,
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					huth_2022.data,
					'Huth et al. (2022) synthetic',
				),
				(
					'D',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					huth_2022.renormalized_passey_2014_data,
					'Passey et al. (2014) biogenic',
				),
			):
				dw = np.array([r['w_d18O'] for r in _data])
				Dw = np.array([r['w_D17O'] for r in _data])
				eDw = np.array([np.nan for r in _data])
				Dw += (
					config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
					- huth_2022.D17O_VSMOW_OF_SLAP
				) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

				dc = np.array([r['c_d18O'] for r in _data])
				Dc = np.array([r['c_D17O'] for r in _data])
				eDc = np.array([np.nan for r in _data])
				# Dc = carb_D17O_adjust(
				# 	dc,
				# 	Dc,
				# 	config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)'),
				# 	config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)'),
				# )

				T = np.array([r['T'] for r in _data]) + 273.15
				ax.errorbar(
					T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
					Dc - Dw - config.SLAP_adjust(CARBREF, 'Wostbrock et al. (2020)'),
					# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
					(eDw**2 + eDc**2)**0.5,
					ecolor = ec,
					mec = mec,
					zorder = 100,
					label = label,
					**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
					)

				reprocessed[CARBREF][label] = dict(
					symbol = symbol,
					mec = mec,
					mfc = mfc,
					dw = dw,
					Dw = Dw,
					dc = (1000 + dc) / 1.03092 - 1000,
					Dc = Dc,
				)

			# ---------------------------------------------------------------- Plot Kelson et al. (2022)
			import kelson_2022
			col = np.array(config.COLORS['Wostbrock'])
			for symbol, mec, mfc, ec, _data, label in (
				(
					'v',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					kelson_2022.data,
					'Kelson et al. (2022) bivalves',
				),
			):
				dw = np.array([r['w_d18O'] for r in _data])
				Dw = np.array([r['w_D17O'] for r in _data])
				eDw = np.array([np.nan for r in _data])
				Dw += (
					config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
					- kelson_2022.D17O_VSMOW_OF_SLAP
				) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

				dc = np.array([r['c_d18O'] for r in _data])
				Dc = np.array([r['c_D17O'] for r in _data])
				eDc = np.array([np.nan for r in _data])
				# Dc = carb_D17O_adjust(
				# 	dc,
				# 	Dc,
				# 	config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)'),
				# 	config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)'),
				# )

				T = np.array([r['T'] for r in _data]) + 273.15
				ax.errorbar(
					T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
					Dc - Dw - config.SLAP_adjust(CARBREF, 'Wostbrock et al. (2020)'),
					# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
					(eDw**2 + eDc**2)**0.5,
					ecolor = ec,
					mec = mec,
					zorder = 100,
					label = label,
					**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
					)

				reprocessed[CARBREF][label] = dict(
					symbol = symbol,
					mec = mec,
					mfc = mfc,
					dw = dw,
					Dw = Dw,
					dc = (1000 + dc) / 1.03092 - 1000,
					Dc = Dc,
				)

			# ---------------------------------------------------------------- Plot Wostbrock et al. (2020)
			import wostbrock_2020
			col = np.array(config.COLORS['Wostbrock'])
			for symbol, mec, mfc, ec, _data, label in (
				(
					'o',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					wostbrock_2020.gca_devils_hole,
					'Wostbrock et al. (2020) Devils Hole',
				),
				(
					config.rdiamond,
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					wostbrock_2020.gca_synthetic_with_ca,
					'Wostbrock et al. (2020) synthetic calcite',
				),
				(
					'd',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					wostbrock_2020.gca_marine_calcite,
					'Wostbrock et al. (2020) marine calcite',
				),
				(
					config.hdiamond,
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					wostbrock_2020.gca_marine_aragonite,
					'Wostbrock et al. (2020) marine aragonite',
				),
			):
				dw = np.array([r['w_d18O'] for r in _data])
				Dw = np.array([r['w_D17O'] for r in _data])
				eDw = np.array([
					# (
					# 	r['w_D17O_SD']/r['Nw']**0.5 if r['Nw'] else 0.005
					# ) * (
					# 	stats.t.ppf(1 - 0.05/2, r['Nw'] - 1) if r['Nw'] > 1 else 1.96
					# )
					# 0.005/r['Nw']**0.5 if r['Nw'] else 0.005
					np.nan
					for r in _data
				])
				Dw += (
					config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
					- wostbrock_2020.D17O_VSMOW_OF_SLAP
				) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

				dc = np.array([r['c_d18O'] for r in _data])
				Dc = np.array([r['c_D17O'] for r in _data])
				eDc = np.array([
					# (
					# 	r['c_D17O_SD']/r['Nc']**0.5 if r['Nc'] else 0.005
					# ) * (
					# 	stats.t.ppf(1 - 0.05/2, r['Nc'] - 1) if r['Nc'] > 1 else 1.96
					# )
					# 0.005/r['Nc']**0.5
					np.nan
					for r in _data
				])
				# Dc = carb_D17O_adjust(
				# 	dc,
				# 	Dc,
				# 	config.SLAP_adjust('IAEA-603', 'Wostbrock et al. (2020)'),
				# 	config.SLAP_adjust('NBS18', 'Wostbrock et al. (2020)'),
				# )

				T = np.array([r['T'] for r in _data]) + 273.15
				ax.errorbar(
					T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
					Dc - Dw - config.SLAP_adjust(CARBREF, 'Wostbrock et al. (2020)'),
					# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
					(eDw**2 + eDc**2)**0.5,
					ecolor = ec,
					mec = mec,
					zorder = 100,
					label = label,
					**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
					)

				reprocessed[CARBREF][label] = dict(
					symbol = symbol,
					mec = mec,
					mfc = mfc,
					dw = dw,
					Dw = Dw,
					dc = (1000 + dc) / 1.03092 - 1000,
					Dc = Dc,
				)

			# ---------------------------------------------------------------- Plot Voarintsoa et al. (2020)
			import voarintsoa_2020
			col = np.array(config.COLORS['Barkan'])
			for symbol, mec, mfc, ec, _data, label in (
				(
					's',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					voarintsoa_2020.data,
					'Voarintsoa et al. (2020) synthetic carbonates',
				),
			):
				dw = np.array([r['d18O_h2o'] for r in _data])
				Dw = np.array([r['D17O_h2o'] for r in _data])
				eDw = np.array([
					np.nan
					for r in _data
				])
				Dw += (
					config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
					- voarintsoa_2020.D17O_VSMOW_OF_SLAP
				) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

				dc = np.array([r['d18O_carb'] for r in _data])
				Dc = np.array([r['D17O_co2_carb'] for r in _data])
				eDc = np.array([
					np.nan
					for r in _data
				])
				# Dc = carb_D17O_adjust(
				# 	(1000 + dc) * 1.03092 - 1000,
				# 	Dc,
				# 	config.SLAP_adjust('IAEA-603-CO2', 'Barkan et al. (2019)'),
				# 	config.SLAP_adjust('NBS18-CO2', 'Barkan et al. (2019)'),
				# )

				T = np.array([r['T'] for r in _data]) + 273.15
				ax.errorbar(
					T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
					Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Barkan et al. (2019)'),
					# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
					(eDw**2 + eDc**2)**0.5,
					ecolor = ec,
					mec = mec,
					zorder = 100,
					label = label,
					**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
					)

				reprocessed[CARBREF][label] = dict(
					symbol = symbol,
					mec = mec,
					mfc = mfc,
					dw = dw,
					Dw = Dw,
					dc = (1000 + dc) / 1.03092 - 1000,
					Dc = Dc,
				)

			if False:
				# ---------------------------------------------------------------- Plot Sha et al. (2020)
				import sha_2020
				col = np.array(config.COLORS['Barkan'])
				for symbol, mec, mfc, ec, _data, label in (
					(
						'^',
						1-(1-col)*1,
						1-(1-col)*0,
						1-(1-col)*0.33,
						sha_2020.data,
						'Sha et al. (2020) speleothems',
					),
				):
					dw = np.array([r['w_d18O'] for r in _data])
					Dw = np.array([r['w_D17O'] for r in _data])
					eDw = np.array([
						np.nan
						for r in _data
					])
					Dw += (
						config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
						- sha_2020.D17O_VSMOW_OF_SLAP
					) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

					dc = np.array([r['c_d18O'] for r in _data])
					Dc = np.array([r['c_D17O'] for r in _data])
					eDc = np.array([
						np.nan
						for r in _data
					])
					# Dc = carb_D17O_adjust(
					# 	(1000 + dc) * 1.03092 - 1000,
					# 	Dc,
					# 	config.SLAP_adjust('IAEA-603-CO2', 'Barkan et al. (2019)'),
					# 	config.SLAP_adjust('NBS18-CO2', 'Barkan et al. (2019)'),
					# )

					T = np.array([r['T'] for r in _data]) + 273.15
					ax.errorbar(
						T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
						Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Barkan et al. (2019)'),
						# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
						(eDw**2 + eDc**2)**0.5,
						ecolor = ec,
						mec = mec,
						zorder = 100,
						label = label,
						**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
						)

					reprocessed[CARBREF][label] = dict(
						symbol = symbol,
						mec = mec,
						mfc = mfc,
						dw = dw,
						Dw = Dw,
						dc = (1000 + dc) / 1.03092 - 1000,
						Dc = Dc,
					)

			# ---------------------------------------------------------------- Plot Roy et al. (2025)
			import roy_2025
			col = np.array(config.COLORS['Roy'])
			for symbol, mec, mfc, ec, _data, label in (
				(
					's',
					1-(1-col)*1,
					1-(1-col)*0,
					1-(1-col)*0.33,
					roy_2025.data,
					'Roy et al. (2025) synthetic',
				),
			):
				dw = np.array([
					(r['starting_water_d18O_VSMOW'] + r['final_water_d18O_VSMOW'])/2
					for r in _data
				])
				Dw = np.array([
					(r['starting_water_D17O_VSMOW'] + r['final_water_D17O_VSMOW'])/2
					for r in _data
				])
				eDw = np.array([
					np.nan
					for r in _data
				])
				Dw += (
					config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
					- roy_2025.D17O_VSMOW_OF_SLAP
				) * np.log(1 + dw/1e3) / np.log(1-55.5e-3)

				dc = (np.array([r['CO2acid_d18O_VSMOW'] for r in _data]) + 1000) / 1.01025 - 1000 # calcite vs VSMOW
				Dc_vs_nbs19 = np.array([r['CO2acid_D17O_VSMOW'] for r in _data]) - roy_2025.D17O_VSMOW_OF_NBS19_CO2
				Dc_vs_iaea603 = (
					Dc_vs_nbs19
					+ config.SLAP_adjust('NBS19-CO2', 'Chaillot et al. (2025)')
					- config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
				)
				eDc = np.array([
					np.nan
					for r in _data
				])
				# Dc = carb_D17O_adjust(
				# 	(1000 + dc) * 1.03092 - 1000,
				# 	Dc,
				# 	config.SLAP_adjust('IAEA-603-CO2', 'Barkan et al. (2019)'),
				# 	config.SLAP_adjust('NBS18-CO2', 'Barkan et al. (2019)'),
				# )

				T = np.array([r['T'] for r in _data]) + 273.15
				if iaea603ref:
					ax.errorbar(
						T**-2 if isT else np.log((1000+dc)/(1000+dw))*1e3,
						Dc_vs_iaea603 - Dw,
						# Dc - Dw - config.SLAP_adjust(f'{CARBREF}-CO2', 'Chaillot et al. (2025)'),
						(eDw**2 + eDc**2)**0.5,
						ecolor = ec,
						mec = mec,
						zorder = 100,
						label = label,
						**(config.kw_errorbar | {'mfc': mfc, 'marker': symbol}),
						)

				reprocessed[CARBREF][label] = dict(
					symbol = symbol,
					mec = mec,
					mfc = mfc,
					dw = dw,
					Dw = Dw,
					dc = (1000 + dc) / 1.03092 - 1000,
					Dc = Dc,
				)


			# ---------------------------------------------------------------- Plot our data
			ax.errorbar(
				data['1/T2'] if isT else data['klna18_observed'],
				data['DD17O'],
				data['95CL_DD17O'],
				ecolor = 'k',
				mec = 'k',
				zorder = 1000,
				label = 'this study',
				**config.kw_errorbar,
				)

			reprocessed[CARBREF]['this study'] = dict(
				symbol = config.kw_errorbar['marker'],
				mec = 'k',
				mfc = 'w',
				dw = data['d18O_VSMOW'],
				Dw = data['D17O_VSMOW'],
				dc = data['d18O_VPDB'],
				Dc = data['D17O_603'] if iaea603ref else (
					data['D17O_603']
					+ config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
					- config.SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
				),
			)

			# ---------------------------------------------------------------- Plot D17O_VSMOW_OF_CARBREF on yaxis

			new_estimate_color = [1-(1-_)/2 for _ in config.COLORS['diseq']]
			bfD, = twinax.plot(0, -D17O_VSMOW_OF_CARBREF_FIT_FROM_THEORY.n,
				'D',
				ms = 4,
				mew = 1,
				mec = 'k',
				mfc = 'w',
				clip_on = False,
				zorder = 1000,
				label = 'This study (best fit)'
			)

			epcolor = np.array(config.COLORS['Ellis'])
			bfD_ep, = twinax.plot(0, -config.SLAP_adjust(CARBREF, 'Ellis & Passey (2023)'),
				'D',
				ms = 4,
				mew = 1,
				mec = 'k',
				mfc = 1-(1-epcolor)/2,
				clip_on = False,
				zorder = 1000,
				label = 'Ellis & Passey (2023)'
			)

			wcolor = np.array(config.COLORS['Wostbrock'])
			bfD_wb, = twinax.plot(0, -config.SLAP_adjust(CARBREF, 'Wostbrock et al. (2020)'),
				'D',
				ms = 4,
				mew = 1,
				mec = 'k',
				mfc = 1-(1-wcolor)/2,
				clip_on = False,
				zorder = 1000,
				label = 'Wostbrock et al. (2020a)'
			)

			# ---------------------------------------------------------------- Plot theoretical model wih best-fit offset
			sqTinv = np.linspace(1e-20, 263.15**-2, 101)
			Tinv = sqTinv**0.5
			Tk = sqTinv**-0.5

			DD = DD17Oeq_afo_T(Tk)
			DDn = np.array([_.n for _ in DD])
			DDe = np.array([_.s for _ in DD]) * 1.96
			# gze = ax.fill_between(
			# 	(
			# 		sqTinv if isT
			# 		else guo_zhou_2019.calcite_water_eq_klna18(Tk)),
			# 	DDn - DDe,
			# 	DDn + DDe,
			# 	lw = 0,
			# 	color = config.COLORS['diseq'],
			# 	alpha = 0.2,
			# 	zorder = 0,
			# )

			gzt, = ax.plot(
				(
					sqTinv if isT
					else guo_zhou_2019.calcite_water_eq_klna18(Tk)),
				DDn,
				'k-',
				lw = 0.8,
				label = 'Guo & Zhou (2019) trend',
				)

		if iaea603ref:
			ax_T.legend(
				loc = 'lower center',
				bbox_to_anchor=(0.5, 0),
				bbox_transform = fig.transFigure,
				labelspacing = 0.2,
				ncols = 2
			)

	fig.savefig('output/compilation.pdf')
	ppl.close(fig)

	R = reprocessed['IAEA-603']

	fig = ppl.figure(figsize = (config.TWO_COLUMN_WIDTH, config.TWO_COLUMN_WIDTH*0.7))
	fig.subplots_adjust(.02, .1, .6, .97, .2, .5)

	ax_w = ppl.subplot(211)
	ax_c = ppl.subplot(212)

	ax_w.set_xlabel(r'$δ^{18}O_\text{VSMOW}$ of waters (‰)')
	ax_c.set_xlabel(r'$δ^{18}O_\text{VPDB}$ of carbonates (‰)')



	y = 0
	yt = []
	for l,r in R.items():
		yt.append(l)
		y += 1
		ax_w.plot(r['dw'], r['dw']*0+y, ls = 'None', marker = r['symbol'], mec = r['mec'], mfc = r['mfc'])
	ax_w.axis([None, None, 0, y+1])
	ax_w.set_yticks(np.arange(len(yt))+1)
	ax_w.set_yticklabels(yt, size = 8)
	ax_w.yaxis.tick_right()
	ax_w.tick_params(axis = 'y', length = 0)

	for x, name in (
		(0, 'VSMOW'),
		(-55.5, 'SLAP'),
	):
		ax_w.axvline(x, color = 'k', lw = 0.8)
		ax_w.text(x, y+1.2, name, ha = 'center', size = 8)
	ax_w.xaxis.set_major_locator(MultipleLocator(10))


	y = 0
	for l,r in R.items():
		print(l)
		y += 1
		ax_c.plot(r['dc'], r['dc']*0+y, ls = 'None', marker = r['symbol'], mec = r['mec'], mfc = r['mfc'])
	ax_c.axis([None, None, 0, y+1])
	ax_c.set_yticks(np.arange(len(yt))+1)
	ax_c.set_yticklabels(yt, size = 8)
	ax_c.yaxis.tick_right()
	ax_c.tick_params(axis = 'y', length = 0)

	for x, name in (
		(-2.20, 'IAEA-603'),
		(-23.01, 'NBS18'),
	):
		ax_c.axvline(x, color = 'k', lw = 0.8)
		ax_c.text(x, y+1.2, name, ha = 'center', size = 8)
	ax_c.xaxis.set_major_locator(MultipleLocator(10))


	fig.savefig('output/all_unknown_vs_anchors.pdf')

	vprint(f'Finished {Path(__file__).name}')
