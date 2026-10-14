#! /usr/bin/env uv run
"""
1. Generate random triple oxygen compositions postulating a "true" Δ17O value of SLAP
2. Compute how these compositions would change if they were standardized assuming a different Δ17O value for SLAP
3. Verify that this effect very closely amounts to a change in Δ17O that varies linearly with ln(1+δ18O):
     newΔ17O - oldΔ17O = (newΔ17O_SLAP - oldΔ17O_SLAP) * ln(1+δ18O) / ln(1+δ18O_SLAP)
"""

import config
import numpy as np
from matplotlib import pyplot as ppl
from matplotlib.ticker import MultipleLocator

np.random.seed(1234)

R18_OF_VSMOW = 0.0020052
R17_OF_VSMOW = 0.00038475
LAMBDA = 0.528

# >> ALL QUANTITIES ARE UNSCALED <<
d18O_OF_SLAP = -55.5e-3
TRUE_D17O_OF_SLAP = 0.
POSTULATED_D17O_OF_SLAP = -11e-6

TRUE_d17_OF_SLAP = (
	np.exp(TRUE_D17O_OF_SLAP)
	* (1 + d18O_OF_SLAP) ** LAMBDA
	- 1
)

MU_d18O = d18O_OF_SLAP / 2
SIGMA_d18O = 10e-3
MU_D17O = 0.
SIGMA_D17O = 100e-6


N = 200

draws = np.random.normal(size = (N,2))

# in order to test that the residuals scale with POSTULATED_D17O_OF_SLAP,
# multiply POSTULATED_D17O_OF_SLAP by an arbitrary factor and save the corresponding plots
for factor in (1, 10):

	_POSTULATED_D17O_OF_SLAP = POSTULATED_D17O_OF_SLAP * factor

	POSTULATED_d17O_OF_SLAP = (
		np.exp(_POSTULATED_D17O_OF_SLAP)
		* (1 + d18O_OF_SLAP) ** LAMBDA
		- 1
	)

	d18O = MU_d18O + draws[:,0] * SIGMA_d18O
	D17O_true = MU_D17O + draws[:,1] * SIGMA_D17O
	d17O_true = np.exp(D17O_true) * (1 + d18O) ** LAMBDA - 1

	d18Oi = np.linspace(0, d18O_OF_SLAP, 101)
	d17Oi_true = (1 + d18Oi) ** LAMBDA - 1

	# 2pt-standardization in d17O space (equivalent to R17 space)
	d17O_reprocessed = d17O_true / TRUE_d17_OF_SLAP * POSTULATED_d17O_OF_SLAP
	D17O_reprocessed = (
		np.log(1 + d17O_reprocessed)
		- np.log(1 + d18O) * LAMBDA
	)
	d17Oi_reprocessed = d17Oi_true / TRUE_d17_OF_SLAP * POSTULATED_d17O_OF_SLAP
	D17Oi_reprocessed = (
		np.log(1 + d17Oi_reprocessed)
		- np.log(1 + d18Oi) * LAMBDA
	)

	Xi = np.array([0, np.log(1 + d18O_OF_SLAP)])
	X = np.log(1 + d18O)
	Y = D17O_reprocessed - D17O_true
	Y_approx = (
		(_POSTULATED_D17O_OF_SLAP - TRUE_D17O_OF_SLAP)
		* X / np.log(1 + d18O_OF_SLAP)
	)

	red = (.9, 0, 0)
	kw_plot_rm = dict(
		ls = 'None',
		marker = 's',
		mec = red,
		mfc = 'w',
		mew = 1,
		ms = 4,
		zorder = 100,
	)
	kw_plot_random = dict(
		ls = 'None',
		marker = '+',
		mec = 'k',
		mew = 0.6,
		ms = 6,
		zorder = 0,
		alpha = 0.5,
	)
	kw_leg = dict(
		fontsize = 8,
		labelspacing = 0.3,
	)
	fig = ppl.figure(figsize = (config.ONE_COLUMN_WIDTH, 4))
	ppl.subplots_adjust(0.21, 0.12, 0.99, 0.99, 0.1, 0.15)

	ax = ppl.subplot(211)

	ax.plot(0, 0, label = 'ref. materials', **kw_plot_rm)
	ax.text(0, 0, 'VSMOW\n\n', color = red, va = 'center', ha = 'center', zorder = 100, linespacing = 1, size = 9)

	ax.plot(
		1e3 * np.log(1+d18O_OF_SLAP),
		1e6 * (_POSTULATED_D17O_OF_SLAP - TRUE_D17O_OF_SLAP),
		**kw_plot_rm)
	ax.text(
		1e3 * np.log(1+d18O_OF_SLAP),
		1e6 * (_POSTULATED_D17O_OF_SLAP - TRUE_D17O_OF_SLAP),
		'\n\nSLAP', color = red, va = 'center', ha = 'center', size = 9,
	)

	ax.plot(
		Xi * 1e3,
		1e6 * (
			(_POSTULATED_D17O_OF_SLAP - TRUE_D17O_OF_SLAP)
			* Xi / np.log(1 + d18O_OF_SLAP)
		), '-', lw = 1, color = red, label = 'linear approximation',
	)

	ax.plot(
		X * 1e3,
		Y * 1e6,
		label = 'random measurements',
		**kw_plot_random,
	)

	ppl.ylabel('Change in $Δ^{17}O$ [ppm]')
	x1, x2, y1, y2 = ppl.axis()
	dx, dy = x2-x1, y2-y1
	e = 0.1
	ppl.axis([
		x1 - dx * e,
		x2 + dx * e,
		y1 - dy * e,
		y2 + dy * e,
	])
	ax.yaxis.set_major_locator(MultipleLocator(3*factor))
	ax.legend(handlelength = 1.2, **kw_leg)
	ax.text(0.97, 0.03, 'a', va  ='bottom', ha = 'right', weight = 'bold', size = 12, transform = ax.transAxes)


	ax = ppl.subplot(212, sharex = ax)
	ax.plot(
		X * 1e3,
		(Y - Y_approx) * 1e6,
		**kw_plot_random,
	)
	ax.plot(
		np.log(1+d18Oi) * 1e3,
		(
			D17Oi_reprocessed
			- (
				(_POSTULATED_D17O_OF_SLAP - TRUE_D17O_OF_SLAP)
				* np.log(1+d18Oi) / np.log(1 + d18O_OF_SLAP)
			)
		) * 1e6,
		'k-',
		dashes = (8,2,2,2),
		lw = 1,
		label = 'residuals for analytes\nwith $Δ^{17}O_{VSMOW} = 0$'
	)
	ppl.xlabel('$10^3$ $\\cdot$ ln(1 + $δ^{18}O_{VSMOW}$)')
	ppl.ylabel('Linear approx. residuals [ppm]')

	y1, y2 = ppl.axis()[-2:]
	dy = y2-y1
	ppl.axis([
		None,
		None,
		y1 - dy * e,
		y2 + dy * e,
	])
	ax.legend(handlelength = 3, **kw_leg)
	ax.text(0.97, 0.97, 'b', va  ='top', ha = 'right', weight = 'bold', size = 12, transform = ax.transAxes)

	if factor != 1:
		fig.savefig(f'output/linear-conversion-x{factor:.0f}.pdf')
	else:
		fig.savefig('output/linear-conversion.pdf')
		with open('output/linear-conversion.toml', 'w') as f:
			f.write(f'N = {N}')
			f.write(f'\nRMSE_ppm = {(Y - Y_approx).std()*1e6:.3f}')
			f.write(f'\n\n[d18O_random]')
			f.write(f'\nmu = {MU_d18O*1e3}')
			f.write(f'\nsigma = {SIGMA_d18O*1e3}')
			f.write(f'\n\n[D17O_random]')
			f.write(f'\nmu = {MU_D17O*1e3}')
			f.write(f'\nsigma = {SIGMA_D17O*1e3}')
