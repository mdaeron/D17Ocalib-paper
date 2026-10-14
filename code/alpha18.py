#! /usr/bin/env uv run

import config
import oneil_1969
import gabitov_2012
import levitt_2018
import daeron_2019
import hayles_2018
import guo_zhou_2019

import numpy as np
import pandas as pd

from pathlib import Path
from config import vprint, hdiamond
from matplotlib import pyplot as ppl
from matplotlib.patches import FancyArrowPatch

# Hales et al. (2018) original klna18 prediction
def klna18_hayles(T):
	return 1000 * hayles_2018.ln_alpha18_calcite_water(T)

# Guo & Zhou (2019) original klna18 prediction
def klna18_guozhou(T):
	return guo_zhou_2019.calcite_water_eq_klna18(T)

if __name__ == '__main__':

	vprint(f'Starting {Path(__file__).name}')

	# Load our data
	data = pd.read_csv('output/calib_results.csv')
	data['1/T'] = (data['T'] + 273.15)**-1
	data['1/T2'] = (data['T'] + 273.15)**-2

	# O'Neil et al. (1969) data
	ondata = pd.DataFrame(oneil_1969.calcite_water_samples)
	ondata = ondata[ondata['T'] > 100]
	ondata['1/T2'] = (ondata['T'] + 273.15)**-2
	ondata['klna18'] = ondata['klna18']
	ondata['95CL_klna18'] = ondata['klna18_cl']

	# Gabitov et al. (2012) data
	gdata = pd.DataFrame(gabitov_2012.calcite_water_samples)
	gdata['1/T2'] = (gdata['T'] + 273.15)**-2

	# Levitt et al. (2018) data
	ldata = pd.DataFrame(levitt_2018.calcite_water_samples)
	ldata = ldata[ldata.logR < -7.8]
	ldata['T'] = ldata['T'].round(0)
	ldata['1/T2'] = (ldata['T'] + 273.15)**-2
	ldata = ldata.groupby('T')[['klna18', '1/T2']].mean()

	# Create figure
	fig = ppl.figure(figsize = (config.TWO_COLUMN_WIDTH, config.TWO_COLUMN_WIDTH/3))
	fig.subplots_adjust(.07, .1, .74, .85)

	# Set up subplot (a)
	ax1 = ppl.subplot(121)
	ppl.xlabel(r'$1\,/\,T^2$', labelpad = 6)
	ppl.ylabel(r'$1000\cdot$ln($^{18}α^{c}_{w}$)')
	Ti = np.array([0, 50, 100, 250, 700][::-1])
	ax1.xaxis.set_ticks((Ti + 273.15)**-2)
	ax1.xaxis.set_ticklabels(
		[f'{t:.0f}$\\,$°C' for t in Ti],
		rotation = -30,
		ha = 'right',
		va = 'bottom',
		size = 9,
		linespacing = 1,
		)
	ax1.tick_params(
		axis = 'x',
		pad = 0,
		top = True,
		labeltop = True,
		bottom = False,
		labelbottom = False,
		)

	# Set up subplot (b)
	ax2 = ppl.subplot(122)
	ppl.xlabel(r'$1\,/\,T^2$', labelpad = 6)
	Ti = np.array([5,15,25,35][::-1])
	ax2.xaxis.set_ticks((Ti + 273.15)**-2)
	ax2.xaxis.set_ticklabels(
		[f'{t:.0f}$\\,$°C' for t in Ti],
		rotation = -30,
		ha = 'right',
		va = 'bottom',
		size = 9,
		linespacing = 1,
		)
	ax2.tick_params(
		axis = 'x',
		pad = 0,
		top = True,
		labeltop = True,
		bottom = False,
		labelbottom = False,
		)

	# Label subplots
	for ax, letter, x, y, ha, va in (
		(ax1, 'a', 1-0.02, 0.02, 'right', 'bottom'),
		(ax2, 'b', 1-0.02, 0.02, 'right', 'bottom'),
	):
		ax.text(x, y, letter, va = va, ha = ha, size = 12, weight = 'bold', transform = ax.transAxes)

	for ax in [ax1, ax2]:
		_lw = 0.8

		# Plot our observations
		if ax == ax2:
			ax.errorbar(
				data['1/T2'],
				data['klna18'],
				data['95CL_klna18'],
				ecolor = 'k',
				zorder = 99,
				**(config.kw_errorbar | {'marker': 'None'}),
			)

		for _data, _marker, _label, _ms in (
			(data[data['T'] < 50], 'o', 'this study (natural calcites)', 5),
			(data[data['T'] > 50], 's', 'this study (synthetic calcites)', 5),
		):
			ax.plot(
				_data['1/T2'],
				_data['klna18'],
				ls = 'None',
				marker = _marker,
				mec = 'k',
				mfc = 'w',
				mew = _lw,
				ms = _ms,
				zorder = 100,
				label = _label,
			)

		# Legend spacer
		if ax == ax2:
			ax.plot([], [], marker = 'None', ls = 'None', label = ' ')

		ax.plot(
			ondata['1/T2'],
			ondata['klna18'],
			ls = 'None',
			marker = 's',
			mec = config.COLORS['ONeil'],
			mfc = 'None',
			mew = _lw,
			ms = 5,
			zorder = 50,
			label = "O'Neil et al. (1969)",
		)

	# Plot Coplen (2007) data
	ax2.plot(
		[(33.7+273.15)**-2],
		[np.log((1000+14.56)/(1000-13.54)) * 1000],
		label = 'Coplen (2007)',
		ls = 'None',
		marker = hdiamond,
		ms = 6,
		mec = config.COLORS['Coplen'],
		mew = _lw,
		zorder = 2,
		mfc = 'None',
		)

	# Plot Daëron et al. (2019) data
	ax2.errorbar(
		[(daeron_2019.samples[_]['T']+273.15)**-2 for _ in daeron_2019.samples],
		[daeron_2019.samples[_]['klna18'] for _ in daeron_2019.samples],
		[daeron_2019.samples[_]['SE_klna18']*1.96 for _ in daeron_2019.samples],
		ecolor = config.COLORS['Daëron'],
		zorder = 0,
		**(config.kw_errorbar | {'marker': 'None', 'elinewidth': _lw, 'capthick': _lw}),
	)

	ax2.plot(
		[(daeron_2019.samples[_]['T']+273.15)**-2 for _ in daeron_2019.samples],
		[daeron_2019.samples[_]['klna18'] for _ in daeron_2019.samples],
		label = 'Daëron et al. (2019)',
		ls = 'None',
		marker = hdiamond,
		ms = 6,
		mec = config.COLORS['Daëron'],
		mew = _lw,
		zorder = 1,
		mfc = 'w',
		)

	# Plot Wostbrock et al. (2020) data
	ax2.plot(
		[(33.7+273.15)**-2],
		[27.994],
		label = 'Wostbrock et al. (2020)',
		ls = 'None',
		marker = hdiamond,
		ms = 6,
		mec = config.COLORS['Wostbrock'],
		mew = _lw,
		zorder = 1,
		mfc = 'None',
	)

	# Legend spacer
	ax2.plot([], [], marker = 'None', ls = 'None', label = ' ')

	# Plot Gabitov et al. (2012) data
	ax2.plot(
		gdata['1/T2'],
		gdata['klna18_slow'],
		ls = 'None',
		marker = config.rdiamond,
		mec = config.COLORS['Gabitov'],
		mfc = 'None',
		mew = _lw,
		ms = 6,
		zorder = 50,
		label = 'Gabitov et al. (2012)',
	)

	# Plot Levitt et al. (2018) data
	ax2.plot(
		ldata['1/T2'],
		ldata['klna18'],
		ls = 'None',
		marker = config.ldiamond,
		mec = config.COLORS['Levitt'],
		mfc = 'None',
		mew = _lw,
		ms = 6,
		zorder = 50,
		label = 'Levitt et al. (2018)',
	)

	# Legend spacer
	ax2.plot([], [], marker = 'None', ls = 'None', label = ' ')

	# Plot theoretical predictions
	_sqTinv = np.linspace(2e-6, (273.15-5)**-2)
	_Tk = _sqTinv**-0.5

	for ax in [ax1, ax2]:

		# Plot Hayles et al. (2018) prediction
		ax.plot(
			_sqTinv,
			klna18_hayles(_Tk),
			'k-',
			lw = 0.7,
			dashes = (2,2),
			label = 'Hayles et al. (2018)',
			zorder = 0,
		)

		# Plot Guo & Zhou (2019) prediction
		ax.plot(
			_sqTinv,
			klna18_guozhou(_Tk),
			'k-',
			lw = 0.7,
			dashes = (5,2),
			label = 'Guo & Zhou (2019)',
			zorder = 0,
		)

	ax2.legend(
		loc = 'center left',
		bbox_to_anchor = (1.0, 0.5),
		fontsize = 8,
		labelspacing = 0.25,
		frameon = False,
		markerfirst = True,
		handlelength = 1.65,
	)

	# add inset rectangle
	ax2r = np.array([39., 4., 27.2, 34.2])
	ax2r[:2] = (ax2r[:2] + 273.15)**-2
	ax2.axis(ax2r)

	inset_color = '0.75'

	ax1r = [0, (273.15-5)**-2, -2, 36]
	ax1.axis(ax1r)
	ax1.plot(
		[ax2r[0], ax2r[0], ax2r[1], ax2r[1], ax2r[0]],
		[ax2r[3], ax2r[2], ax2r[2], ax2r[3], ax2r[3]],
		'-',
		color = inset_color,
		lw = 1,
		zorder = -1,
	)
	# ax1.text(ax2r[1]-2e-7, ax2r[2]+.1, 'b', va = 'bottom', ha = 'right', size = 9, weight = 'bold', color = inset_color)
	ax1.add_artist(
		FancyArrowPatch(
			posA = [ax2r[1], (ax2r[2]+ax2r[3])/2],
			posB = [ax2r[1]+35e-7, (ax2r[2]+ax2r[3])/2],
			arrowstyle = '->',
			shrinkA = 2,
			shrinkB = 0,
			mutation_scale = 12,
			lw = 1,
			color = inset_color,
			clip_on = False,
		)
	)

	# Save figure
	fig.savefig('output/alpha18.pdf')
	# ppl.show()
	ppl.close(fig)

	vprint(f'Finished {Path(__file__).name}')
