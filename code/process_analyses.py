#! /usr/bin/env uv run

import stdz
import config

import numpy as np
import pandas as pd

from pathlib import Path
from config import vprint
from matplotlib import transforms
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator


def check_d626b_outliers(
	data: pd.DataFrame,
	p_grubbs: float = 0.95,
	flag: bool = True,
	log: bool = True,
	rename: bool = False,
) -> None:

	outliers = grubbs_two_sided_test(data['d626b'].to_numpy(), alpha = 1 - p_grubbs)

	if flag:
		data['d626b_outlier'] = outliers

	if log:
		for status, (uid, r) in zip(outliers, data.iterrows()):
			if status:
				vprint(f"Analysis {uid} ({r['Sample']}) is a d626b outlier (d626b = {r['d626b']:.1f} ‰)")

	if rename:
		data.loc[outliers, 'Sample'] = data.loc[outliers, 'Sample'] + '_outlier_' + data.loc[outliers].index


def grubbs_two_sided_test(
	x: np.typing.ArrayLike,
	alpha: float = 0.05,
) -> np.array:
	"""
	Return array of booleans indicating which elements of `x` are
	Grubbs outliers for a p-value of (1-alpha).
	"""

	from scipy.stats import t

	X = np.asarray(x)
	mask = np.array([True for _ in X])

	while True:
		Xmean = X[mask].mean()
		Xsd = X[mask].std(ddof = 1)
		N = X[mask].size
		Z = abs(X - Xmean) / Xsd
		Z[~mask] = 0
		G_calculated = Z.max()
		t_critical = t.isf(alpha / 2 / (2 * N), N - 2)
		G_critical = ((N - 1) / np.sqrt(N)) * np.sqrt(t_critical**2 / (N - 2 + t_critical**2))

		if G_calculated > G_critical:
			mask[Z.argmax()] = False
		else:
			break

	return ~mask


def d626b_plot(
	data: stdz.Dataset | pd.DataFrame,
	p_grubbs: float = 0.95,
	filename: str = 'output/d626b_outliers.pdf',
) -> None:

	if isinstance(data, pd.DataFrame):
		X = data.loc[:, ['d626b', 'D17O_VSMOW_residual', 'd626b_outlier']]
	elif isinstance(data, stdz.Dataset):
		X = data.data.loc[:, ['d626b', 'D17O_VSMOW_residual', 'd626b_outlier']]

	fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, 2.8))
	ppl.subplots_adjust(.2, .2, .97, .8)

	_all, = ppl.plot(
		X['d626b'],
		X['D17O_VSMOW_residual'] * 1000,
		'ko',
		mfc = 'w',
		ms = 5,
		mew = 1,
	)

	_outliers, = ppl.plot(
		X.loc[X.d626b_outlier, 'd626b'],
		X.loc[X.d626b_outlier, 'D17O_VSMOW_residual'] * 1000,
		'x',
		mew = 1,
		ms = 8,
		mfc = 'none',
		mec = 'k',
	)

	if True in X.d626b_outlier:
		ppl.legend(
			[(_all, _out)],
			[f'Grubbs outliers (p = {p_grubbs})\nbased only on δ$_{{626b}}$ values'],
			fontsize = 9,
			handlelength = 1.5,
			loc = 'lower center',
			bbox_to_anchor = (.5, 1.01),
			)

	ppl.axhline(0, color = 'k', lw = 0.75, dashes = (4,2))
	ppl.axvline(0, color = 'k', lw = 0.75, dashes = (4,2))
	ppl.xlabel('δ$_{626b}$ (‰)')
	ppl.ylabel('Δ$^{17}$O residuals (ppm)')
	ppl.margins(0.2)

	fig.savefig(filename)
	ppl.close(fig)


def plot_DVHW(
	df_old: pd.DataFrame,
	df_new: pd.DataFrame,
	filename: str = 'output/DVHW.pdf',
) -> None:
	fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, 3.7))
	ppl.subplots_adjust(.22, .12, .84, .93)
	ax1 = ppl.subplot(211)
	ppl.ylabel('$δ^{18}O_{VSMOW}$ (‰)')
	ax2 = ppl.subplot(212, sharex = ax1)
	ppl.ylabel('$Δ^{17}O_{VSMOW}$ (ppm)')

	df_old.loc[df_old.index.str.startswith('DVHW'), 'Depth'] = df_old.loc[df_old.index.str.startswith('DVHW')].index.str[4:].astype(float) * 0.3048
	avg = df_new.loc['DVHW', :].to_dict()
	ax1.xaxis.set_major_locator(MultipleLocator(10))

	ax1.text(
		1,
		1.07,
		'Devils Hole Cave #2',
		transform = ax1.transAxes,
		weight = 'bold',
		size = 9,
		ha = 'right',
	)
	ax1.axhspan(
		avg['water_d18O_VSMOW'] - avg['95CL_water_d18O_VSMOW'] / 2,
		avg['water_d18O_VSMOW'] + avg['95CL_water_d18O_VSMOW'] / 2,
		color = 'r',
		alpha = .15,
		lw = 0,
	)
	ax1.text(
		1.02,
		avg['water_d18O_VSMOW'],
		f"{avg['water_d18O_VSMOW']:.2f}\n $±\\,${avg['95CL_water_d18O_VSMOW']:.2f}".replace('-', '$-$'),
		color = 'r',
		size = 8.5,
		va = 'center',
		transform = transforms.blended_transform_factory(
			ax1.transAxes,
			ax1.transData,
		),
	)
	ax1.errorbar(
		df_old['Depth'],
		df_old['water_d18O_VSMOW'],
		df_old['95CL_water_d18O_VSMOW'],
		ecolor = 'r',
		mec = 'r',
		**config.kw_errorbar,
	)
	df_old = df_old[df_old.Depth>0]

	ax2.axhspan(
		(avg['water_D17O_VSMOW'] - avg['95CL_water_D17O_VSMOW'] / 2) * 1000,
		(avg['water_D17O_VSMOW'] + avg['95CL_water_D17O_VSMOW'] / 2) * 1000,
		color = 'm',
		alpha = .25,
		lw = 0,
	)
	ax2.text(
		1.02,
		avg['water_D17O_VSMOW']*1000,
		f"{1000*avg['water_D17O_VSMOW']:.0f}$\\,±\\,${1000*avg['95CL_water_D17O_VSMOW']:.0f}\nN = {avg['N']:.0f}",
		color = 'm',
		size = 8.5,
		va = 'center',
		transform = transforms.blended_transform_factory(
			ax2.transAxes,
			ax2.transData,
		),
	)
	ax2.errorbar(
		df_old['Depth'],
		df_old['water_D17O_VSMOW'] * 1000,
		df_old['95CL_water_D17O_VSMOW'] * 1000,
		ecolor = 'm',
		mec = 'm',
		**config.kw_errorbar,
	)

	for x,y,ey,n in zip(
		df_old['Depth'],
		df_old['water_D17O_VSMOW']*1e3,
		df_old['95CL_water_D17O_VSMOW']*1e3,
		df_old['N'],
	):
		ax2.text(
			x-0.3,
			y-ey-3,
			f'N$\\,=\\,${n}',
			va = 'top',
			ha = 'center',
			size = 8,
			color = 'k',
		)

	ax1.margins(0.22)
	ax2.margins(0.22)
	ax1.text(1-0.02,  0.04,'a',weight='bold', size = 9, ha = 'right', transform = ax1.transAxes)
	ax2.text(1-0.02,1-0.04,'b',weight='bold', size = 9, ha = 'right', transform = ax2.transAxes, va = 'top')
	ppl.xlabel('Water depth (m)')
	fig.savefig(filename)
	ppl.close(fig)


if __name__ == '__main__':

	vprint(f'Starting {Path(__file__).name}')

	# READ ANLYTICAL DATA
	df = pd.read_csv('inputs/rawdata.csv', index_col = 'UID')

	# Define sessions
	df.loc[df.Run.str.startswith('2023-06'), 'Session'] = '2023-06'
	df.loc[df.Run.str.startswith('2023-08'), 'Session'] = '2023-08'
	df.loc[df.Run == '2023-07-31', 'Session'] = '2023-08'
	df.loc[df.Run.str.startswith('2024-07'), 'Session'] = '2024-07'

	# Exclude two analytical outliers
	df = df.drop(df[df.index == '5D727D28'].index)
	df = df.drop(df[df.index == 'CE4051F7'].index)

	# --- QUALITY CONTROL: CHECK FOR d626b OUTLIERS --- #
	vprint('Quality control')

	# Define water anchors for a rough, first-order standardziation
	water_anchors = {
		'VSMOW2': {'d18O_VSMOW':   0.0, 'D17O_VSMOW': 0.000},
		'SLAP2':  {'d18O_VSMOW': -55.5, 'D17O_VSMOW': 0.000},
	} # this is just for screening so which normalization scale is used doesn't really matter

	# Adjust D17O in the case of non-zero SLAP value:
	if config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE:
		for w in water_anchors:
			water_anchors[w]['D17O_VSMOW'] += (
				np.log(1+water_anchors[w]['d18O_VSMOW']/1e3)
				/ np.log(1-55.5e-3)
				* config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE
			)

	# Flag for outliers within each session
	for session in df.Session.unique():
		sdf = df.loc[df.Session == session, ['Sample', 'd626b']]
		check_d626b_outliers(sdf, p_grubbs = config.P_GRUBBS)
		df.loc[df.Session == session, 'd626b_outlier'] = sdf['d626b_outlier']

	# Perform rough standardization (needed for plotting)
	data = stdz.Dataset(df)
	data.standardize_D17O(
		'd628', 'd627', 'd18O_VSMOW', 'd17O_VSMOW',
		anchors = water_anchors,
		relative_sigma_precision = 1e-2,
	)
	d626b_plot(data, p_grubbs = config.P_GRUBBS)

	# Delete outliers
	df = df[df.d626b_outlier == False]

	# --- STANDARDIZE CARBONATES --- #
	vprint('Standardize carbonates')

	carb_anchors = {
		'IAEA603': {
			'd18O_VPDB': -2.37,
			'D17O_603': 0.,
		},
		'NBS18': {
			'd18O_VPDB': -23.01,
			'D17O_603': config.D17O_VSMOW_OF_NBS18_CO2_CHAILLOT - config.D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT
		},
	}

	carb_df = df[df.Type == 'carbonate'].copy()
	carb_ds = stdz.Dataset(carb_df)

	carb_ds.standardize_D17O(
		anchors = carb_anchors,
		d18O_key_in = 'd628',
		d17O_key_in = 'd627',
		d18O_key_out = 'd18O_VPDB',
		d17O_key_out = 'd17O_603',
		relative_sigma_precision = config.RELATIVE_SIGMA_PRECISION,
	)
	carb_ds.samples.to_csv('output/carb_samples.csv', float_format = '%.6f')
	carb_ds.data.to_csv('output/carb_analyses.csv', float_format = '%.6f')
	# config.df_to_toml(carb_ds.samples, 'output/carb_samples.toml')

	# Plot residuals as function of time
	fig = ppl.figure(figsize = (11, 8))
	carb_ds.plot_residuals('D17O_603', ylabel = '$Δ^{17}O_{603}$ residuals   [‰]')
	fig.savefig('output/carb_D17O_residuals.pdf')

	fig = ppl.figure(figsize = (11, 8))
	carb_ds.plot_residuals('d18O_VPDB', ylabel = '$δ^{18}O_{VPDB}$ residuals   [‰]')
	fig.savefig('output/carb_d18O_residuals.pdf')

	# Plot Δ17O residuals vs δ18O residuals
	fig = ppl.figure(figsize = (5, 5))
	carb_ds.plot_D17O_residuals(
		d18O_key_out = 'd18O_VPDB',
		D17O_key_out = 'D17O_603',
	)
	fig.savefig('output/carb_kde_D17O_residuals.pdf')

	fr = carb_ds.standardization['D17O_603']['fitresult']
	cdiff = fr.uvars['D17O_DVH2']-fr.uvars['D17O_LGB2']
	with open('output/processing-results.toml', 'w') as f:
		f.write('\n[standardization.carbonates]\n')
		f.write(f'N = {len(carb_ds.data)}\n')
		Nf = len(carb_ds.data) - len(carb_ds.samples)
		f.write(f'Nf = {Nf}\n')
		r = carb_ds.data.D17O_603_residual.to_numpy()
		rmse = ((r**2).sum()/(len(carb_ds.data) - len(carb_ds.samples)))**0.5
		f.write(f'D17O_603_sigma = {rmse:.6f}\n')
		r = carb_ds.data.d18O_VPDB_residual.to_numpy()
		rmse = ((r**2).sum()/(len(carb_ds.data) - len(carb_ds.samples)))**0.5
		f.write(f'd18O_VPDB_sigma = {rmse:.6f}\n')


	# --- STANDARDIZE WATERS --- #
	vprint('Standardize waters (first run)')

	water_df = df[df.Type == 'water'].copy()

	# Specify initial composition of equilibrated CO2
	match config.CO2g_H2O_FRACTIONATION:
		case 'Guo & Zhou (2019)':
			starting_CO2 = {'d18O_VSMOW': 24.71, 'D17O_VSMOW': -0.0835}
		case 'Barkan & Luz (2012)':
			starting_CO2 = {'d18O_VSMOW': 24.25, 'D17O_VSMOW': -0.150}

	water_df['eqCO2_d18O_VSMOW'] = starting_CO2['d18O_VSMOW']
	water_df['eqCO2_D17O_VSMOW'] = starting_CO2['D17O_VSMOW']
	water_df['eqCO2_d13C_VPDB'] = -3.7

	# Standardize waters, treating DVHW samples as separate samples
	water_ds = stdz.Dataset(water_df)
	water_ds.standardize_water_eq_CO2_D17O(
		water_anchors,
		alpha18_CO2g_H2O = config.alpha18_CO2g_H2O,
		theta17_CO2g_H2O = config.theta17_CO2g_H2O,
		relative_sigma_precision = config.RELATIVE_SIGMA_PRECISION,
	)
	water_samples_old = water_ds.samples.copy()
	water_ds.samples.to_csv('output/water_samples_before_merging.csv', float_format = '%.6f')
	water_ds.data.to_csv('output/water_analyses_before_merging.csv', float_format = '%.6f')

	def rmswd(ufloats):

		from scipy import stats
		from uncertainties import covariance_matrix

		_x = np.array([u for u in ufloats])
		x = np.array([u.n for u in ufloats])
		C = np.array(covariance_matrix(ufloats))
		Cinv = np.linalg.inv(C)

		# compute GLS weighted mean
		ones = np.ones_like(x)
		xbar = (ones @ Cinv @ x) / (ones @ Cinv @ ones)
		_xbar = (ones @ Cinv @ _x) / (ones @ Cinv @ ones)

		resid = x - xbar
		chi2 = resid @ Cinv @ resid
		dof = len(x) - 1
		rmswd_val = np.sqrt(chi2 / dof)
		p_value = stats.chi2.sf(chi2, dof)

		return rmswd_val, p_value

	R = water_ds.standardization['water_D17O_VSMOW']['fitresult']
	X = [R.uvars[_] for _ in R.uvars if 'DVHW' in _ and 'd18O' in _]
	rmswdX, pvX = rmswd(X)
	Y = [R.uvars[_] for _ in R.uvars if 'DVHW' in _ and 'D17O' in _]
	rmswdY, pvY = rmswd(Y)

	vprint('Standardize waters (second run)')

	# Treat all DVHW samples as a single sample
	water_df.loc[water_df.Sample.str.startswith('DVHW'), 'Sample'] = 'DVHW'

	# Standardize waters again
	water_ds = stdz.Dataset(water_df)
	water_ds.standardize_water_eq_CO2_D17O(
		water_anchors,
		alpha18_CO2g_H2O = config.alpha18_CO2g_H2O,
		theta17_CO2g_H2O = config.theta17_CO2g_H2O,
		relative_sigma_precision = config.RELATIVE_SIGMA_PRECISION,
	)
	water_ds.samples.to_csv('output/water_samples.csv', float_format = '%.6f')
	water_ds.data.to_csv('output/water_analyses.csv', float_format = '%.6f')
	# config.df_to_toml(water_ds.samples, 'output/water_samples.toml')

	with open('output/processing-results.toml', 'a') as f:
		f.write('\n[standardization.waters]\n')
		f.write(f'N = {len(water_ds.data)}\n')
		Nf = len(water_ds.data) - len(water_ds.samples)
		f.write(f'Nf = {Nf}\n')
		r = water_ds.data.CO2_D17O_VSMOW_residual.to_numpy()
		rmse = ((r**2).sum()/(len(water_ds.data) - len(water_ds.samples)))**0.5
		f.write(f'D17O_VSMOW_sigma = {rmse:.6f}\n')
		r = water_ds.data.CO2_d18O_VSMOW_residual.to_numpy()
		rmse = ((r**2).sum()/(len(water_ds.data) - len(water_ds.samples)))**0.5
		f.write(f'd18O_VSMOW_sigma = {rmse:.6f}\n')
		f.write('\n[DVHW.d18O_RMSWD]\n')
		f.write(f'value = {rmswdX:.6f}\n')
		f.write(f'p = {pvX:.6f}\n')
		f.write('\n[DVHW.D17O_RMSWD]\n')
		f.write(f'value = {rmswdY:.3e}\n')
		f.write(f'p = {pvY:.3e}\n')

	fig = ppl.figure(figsize = (11, 8))
	water_ds.plot_residuals('CO2_D17O_VSMOW', ylabel = '$Δ^{17}O_{VSMOW}$ residuals   [‰]')
	fig.savefig('output/water_D17O_residuals.pdf')

	fig = ppl.figure(figsize = (11, 8))
	water_ds.plot_residuals('CO2_d18O_VSMOW', ylabel = '$δ^{18}O_{VSMOW}$ residuals   [‰]')
	fig.savefig('output/water_d18O_residuals.pdf')

	fig = ppl.figure(figsize = (5, 5))
	ppl.subplots_adjust(.15, .15, .95, .95)
	water_ds.plot_D17O_residuals(
		d18O_key_out = 'CO2_d18O_VSMOW',
		D17O_key_out = 'CO2_D17O_VSMOW',
	)
	fig.savefig('output/water_kde_D17O_residuals.pdf')

	# Plot DVHW variability
	vprint('Plot DVHW variability')
	plot_DVHW(water_samples_old, water_ds.samples)



	# rmswdX = ((((X-X.mean())/eX)**2).sum()/(X.size-1))**0.5 / 1.96
	# print(f'{rmswdX =:.2}')
	# rmswdY = ((((Y-Y.mean())/eY)**2).sum()/(Y.size-1))**0.5 / 1.96
	# print(f'{rmswdY =:.2}')


	fr = water_ds.standardization['water_D17O_VSMOW']['fitresult']
	wdiff = fr.uvars['D17O_DVHW']-fr.uvars['D17O_LGBW']
	cwdiff = cdiff - wdiff
	with open('output/processing-results.toml', 'a') as f:
		f.write('\n[D17Odiff_DVH_LGB]\n')
		f.write(f'value = {cwdiff.n:.6f}\n')
		f.write(f'se = {cwdiff.s:.6f}\n')
		f.write(f'c95 = {cwdiff.s*1.96:.6f}\n')

	# --- COMBINE RESULTS --- #
	vprint('Combine water and carbonate results')

	# Prepare calibration data
	results = pd.DataFrame(dict(
		Sample =      ['LGB',  'DVH',  'CC75',  'CC120',  'CC170',  'CC200',  'CC250'],
		T =           [  7.9,   33.7,     75.,     120.,     170.,     200.,     250.],
		CarbSample =  ['LGB2', 'DVH2', 'CC75',  'CC120',  'CC170',  'CC200',  'CC250'],
		WaterSample = ['LGBW', 'DVHW', 'CC75W', 'CC120W', 'CC170W', 'CC200W', 'CC250W'],
	))

	results.set_index('WaterSample', inplace = True)
	results['N_water'] = water_ds.samples['N']
	results['d18O_VSMOW'] = water_ds.samples['water_d18O_VSMOW']
	results['95CL_d18O_VSMOW'] = water_ds.samples['95CL_water_d18O_VSMOW']
	results['D17O_VSMOW'] = water_ds.samples['water_D17O_VSMOW']
	results['95CL_D17O_VSMOW'] = water_ds.samples['95CL_water_D17O_VSMOW']

	results.set_index('CarbSample', inplace = True)
	results['N_carb'] = carb_ds.samples['N']
	results['d18O_VPDB'] = carb_ds.samples['d18O_VPDB']
	results['95CL_d18O_VPDB'] = carb_ds.samples['95CL_d18O_VPDB']
	results['D17O_603'] = carb_ds.samples['D17O_603']
	results['95CL_D17O_603'] = carb_ds.samples['95CL_D17O_603']

	results.set_index('Sample', inplace = True)

	if config.FORCE_LGB_WATER_d18O_TO_LONG_TERM_AVERAGE:
		# δ18O values from Daëron et al. (2019)
		results.loc['LGB', 'd18O_VSMOW'] = -7.39
		results.loc['LGB', '95CL_d18O_VSMOW'] = 1.96 * 0.09

	results['alpha18'] = (1000 + results['d18O_VPDB']) / (1000 + results['d18O_VSMOW']) * 1.03092
	results['klna18'] = 1000 * np.log(results['alpha18'])
	results['95CL_klna18'] = 1000 * (
		(results['95CL_d18O_VPDB'] / (1000 + results['d18O_VPDB']))**2
		+ (results['95CL_d18O_VSMOW'] / (1000 + results['d18O_VSMOW']))**2
	)**0.5

	results['DD17O'] = results['D17O_603'] - results['D17O_VSMOW']
	results['95CL_DD17O'] = (results['95CL_D17O_603']**2 + results['95CL_D17O_VSMOW']**2)**0.5

	### Compute offset from Guo & Zhou (2019), expressed in weighed deviation (number of sigmas)
	# klna18_obs = results['klna18'].to_numpy()
	# klna18_obs_sigma = results['95CL_klna18'].to_numpy() / 1.96
	# klna18_theory = (results['T'] + 273.15).map(guo_zhou_2019.calcite_water_eq_klna18)
	# results['klna18_mismatch'] = abs(klna18_obs - klna18_theory) / klna18_obs_sigma
	# results['klna18_diseq'] = klna18_obs - klna18_theory

	results['T'] = results['T'].map(lambda x: f'{x:.1f}')
	results.to_csv('output/calib_results.csv', float_format = '%.6f')

	fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, 5))
	ppl.subplots_adjust(.18, .13, .96, .98, .5, .3)
	ax1 = ppl.subplot(211)
	ax1.set_xlabel('$δ^{18}O_\\text{VSMOW}$ [‰]')
	ax1.set_ylabel('$Δ^{17}O_\\text{VSMOW}$ [‰]')
	color = np.array((0,0,0.9))
	for marker, label, X in (
		('o', 'unknowns', water_ds.samples[water_ds.samples['95CL_water_D17O_VSMOW'] > 0]),
		('s', 'anchors', water_ds.samples[water_ds.samples['95CL_water_D17O_VSMOW'].isna()]),
	):
		match label:
			case 'unknowns':
				ax1.errorbar(
					X['water_d18O_VSMOW'],
					X['water_D17O_VSMOW'],
					X['95CL_water_D17O_VSMOW'],
					ecolor = color,
					mec = color,
					label = label,
					**(config.kw_errorbar | {'marker': marker}),
				)
			case 'anchors':
				ax1.plot(
					X['water_d18O_VSMOW'],
					X['water_D17O_VSMOW'],
					marker,
					mec = color,
					mfc = 1-(1-color)/3,
					ms = config.kw_errorbar['ms'],
					mew = config.kw_errorbar['mew'],
					label = label,
				)
	ax1.yaxis.set_major_locator(MultipleLocator(0.02))
	kw_legend = dict(
		fontsize = 9,
	)
	ax1.legend(**kw_legend)

	ax2 = ppl.subplot(212)
	ax2.set_xlabel('$δ^{18}O_\\text{VPDB}$ [‰]')
	ax2.set_ylabel('$Δ^{17}O_\\text{IAEA-603}$ [‰]')
	color = np.array((0.9,0,0))
	for marker, label, X in (
		('o', 'unknowns', carb_ds.samples[carb_ds.samples['95CL_D17O_603'] > 0]),
		('s', 'anchors', carb_ds.samples[carb_ds.samples['95CL_D17O_603'].isna()]),
	):
		match label:
			case 'unknowns':
				ax2.errorbar(
					X[X.index != 'CAMORA']['d18O_VPDB'],
					X[X.index != 'CAMORA']['D17O_603'],
					X[X.index != 'CAMORA']['95CL_D17O_603'],
					ecolor = color,
					mec = color,
					label = label,
					**(config.kw_errorbar | {'marker': marker}),
				)
			case 'anchors':
				ax2.plot(
					X['d18O_VPDB'],
					X['D17O_603'],
					marker,
					mec = color,
					mfc = 1-(1-color)/3,
					ms = config.kw_errorbar['ms'],
					mew = config.kw_errorbar['mew'],
					label = label,
				)
	ax2.yaxis.set_major_locator(MultipleLocator(0.04))
	ax2.legend(**kw_legend)

	fig.savefig('output/unknowns_vs_anchors.pdf')

	vprint(f'Finished {Path(__file__).name}')
