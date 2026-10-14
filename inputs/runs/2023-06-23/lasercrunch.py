#! /usr/bin/env uv run

from pathlib import Path
from csv import DictReader
from math import isnan
import numpy as np
from matplotlib import pyplot as ppl
from matplotlib.transforms import blended_transform_factory


INPUT_DIR       = 'input'
INPUT_FMT       = ['csv', 'dat']
WG_NAME         = 'Ref Bottle'
T_MULTIPLIER    = 24*60
T_START         = -1e9
T_STOP          = +1e9
P_MULTIPLIER    = 1.333
GAP_WIDTH       = 4.0
IGNORE_FIRST    = 0
SIGMA_FILTER    = 0.1
BLOCK_COLOR     = [.95, .95, .95]
BLOCK_COLOR_WG  = [ .9,  .95,   1]
ISO_FIELDS      = (
	'626',
	'636',
	'628',
	'627',
	'626b',
	)
ISO_COLORS = {
	'626' : [  0,  0,  0],
	'636' : [  0, .7,  0],
	'628' : [.95,  0,  0],
	'627' : [  1,  0,  1],
	'626b': [.35, .7, .7],
	}
FLOAT_FIELDS = ('t', 'P', *ISO_FIELDS)

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


def read_input(input_dir = INPUT_DIR, input_fmt = INPUT_FMT):
	"""
	Reads input files and returns a list of dicts
	"""

	data = []

	for fmt in input_fmt:
		for file in sorted(Path(input_dir).glob(f'*.{fmt}')):
			with open(file) as fid:
				data += list(DictReader(fid, delimiter = '\t'))

	run = Path('.').absolute().stem

	return data, run


def cleanup(data):
	"""
	Cleanup by removing lines with missing fields
	"""
	out = [r.copy() for r in data
		if all([is_number(r[f]) and not isnan(float(r[f])) for f in FLOAT_FIELDS])
		]

	t0 = float(out[0]['t'])
	for r in out:
		if r['SampName'] == WG_NAME:
			r['SampName'] = 'WG'
		for f in FLOAT_FIELDS:
			r[f] = float(r[f])
		r['P'] *= P_MULTIPLIER
		r['t'] -= t0
		r['t'] *= T_MULTIPLIER

	return [r for r in out if r['t'] > T_START and r['t'] < T_STOP]


def chop(data):

	out = []
	currentsample = ''
	t_prev = data[0]['t']
	for r in data:
		if r['SampName'] != currentsample or (r['t'] - t_prev) > GAP_WIDTH:
			currentsample = r['SampName']
			out.append({})

		t_prev = r['t']

		if not out[-1]:
			out[-1]['Sample'] = r['SampName']
			for k in FLOAT_FIELDS:
				out[-1][k] = [r[k]]
		else:
			for k in FLOAT_FIELDS:
				out[-1][k].append(r[k])

	for r in out:
		for k in FLOAT_FIELDS:
			r[k] = np.array(r[k])
		r['mask'] = np.array([k >= IGNORE_FIRST for k,t in enumerate(r['t'])])

	return out

def crunch(data):

	out = [r.copy() for r in data]

	for r in out:
		for k in ISO_FIELDS[1:]:
			r['_R'+k] = r[k] / r[ISO_FIELDS[0]]
			m = np.median(r['_R'+k])
			r['mask'] = r['mask'] & (abs((r['_R'+k] / m - 1) * 1e3) < SIGMA_FILTER)

	wg = {
		'R'+k: np.array([_ for r in out for _ in r['_R'+k][r['mask']] if r['Sample'] == 'WG']).mean()
		for k in ISO_FIELDS[1:]
		}

	for r in out:
		for k in ISO_FIELDS[1:]:
			r['_d'+k+'_drifting'] = (r['_R'+k] / wg['R'+k] - 1) * 1e3
		r['t_start'] = r['t'][r['mask']].min()
		r['t_stop'] = r['t'][r['mask']].max()

	return out


def dedrift(data):

	out = [r.copy() for r in data]

	for k,r in enumerate(out):
		if r['Sample'] == 'WG':
			for i in ISO_FIELDS[1:]:
				r['_d'+i] = r['_d'+i+'_drifting'] - r['_d'+i+'_drifting'][r['mask']].mean()
		else:
			if k > 0 and k < (len(out)-1) and out[k-1]['Sample'] == 'WG' and out[k+1]['Sample'] == 'WG':
				for i in ISO_FIELDS[1:]:
					x = np.hstack((
						out[k-1]['t'][out[k-1]['mask']],
						out[k+1]['t'][out[k+1]['mask']],
						))
					y = np.hstack((
						out[k-1]['_d'+i+'_drifting'][out[k-1]['mask']],
						out[k+1]['_d'+i+'_drifting'][out[k+1]['mask']],
						))
					a,b = np.polyfit(x, y, 1)
					r['_d'+i] = (
						(1000 + r['_d'+i+'_drifting']) / (1000 + a * r['t'] + b) - 1
						) * 1e3
			else:
				raise ValueError(f'Block {k+1} is not bracketed by WG blocks.')

		for i in ISO_FIELDS[1:]:
			r['d'+i] = r['_d'+i][r['mask']].mean()

	return out


def plot_all(data, figsize = (16,8.5)):

	kwargs_valid = dict(mew = 0, ms = 3, alpha = 0.2, marker = 'o', ls = 'None')
	kwargs_rejected = dict(mew = 1, ms = 3, alpha = 1, marker = 'x', ls = 'None', mec = 'k')
	kwargs_blockmean = dict(mew = 1.5, mfc = 'w', ms = 4, alpha = 1, marker = 'D', ls = 'None')

	fig = ppl.figure(figsize = figsize)
	ax_p = ppl.subplot(511)
	ax_c = ppl.subplot(512, sharex = ax_p)
	ax_r = ppl.subplot(513, sharex = ax_p)
	ax_d = ppl.subplot(514, sharex = ax_p)
	ax_residuals = ppl.subplot(515, sharex = ax_p)
	ppl.subplots_adjust(top = 0.85)

	ppl.sca(ax_p)
	ppl.ylabel('P / mbar')

	btrans = blended_transform_factory(ppl.gca().transData, ppl.gca().transAxes)
	for r in data:
		ppl.axvspan(r['t_start'], r['t_stop'], color = BLOCK_COLOR_WG if r['Sample'] == 'WG' else BLOCK_COLOR, lw = 1)
		ppl.text((2*r['t_start'] + r['t_stop'])/3, 1.03, r['Sample'], va = 'bottom', ha = 'left', rotation = 60, size = 7, transform = btrans)

	for r in data:
		ppl.plot(r['t'][r['mask']], r['P'][r['mask']], 'k', **kwargs_valid)
		ppl.plot(r['t'][~r['mask']], r['P'][~r['mask']], **kwargs_rejected)

	ppl.sca(ax_c)
	ppl.ylabel('concentration')

	for r in data:
		ppl.axvspan(r['t_start'], r['t_stop'], color = BLOCK_COLOR_WG if r['Sample'] == 'WG' else BLOCK_COLOR, lw = 1)

	for r in data:
		for i in ISO_FIELDS:
			ppl.plot(r['t'][r['mask']], r[i][r['mask']], mfc = ISO_COLORS[i], **kwargs_valid)
			ppl.plot(r['t'][~r['mask']], r[i][~r['mask']], **kwargs_rejected)

	ppl.sca(ax_r)
	ppl.ylabel('drifting delta')

	for r in data:
		ppl.axvspan(r['t_start'], r['t_stop'], color = BLOCK_COLOR_WG if r['Sample'] == 'WG' else BLOCK_COLOR, lw = 1)

	for r in data:
		for i in ISO_FIELDS[1:]:
			ppl.plot(r['t'][r['mask']], r['_d'+i+'_drifting'][r['mask']], mfc = ISO_COLORS[i], **kwargs_valid)
			ppl.plot(r['t'][~r['mask']], r['_d'+i+'_drifting'][~r['mask']], **kwargs_rejected)

	ppl.sca(ax_d)
	ppl.ylabel('delta')

	for r in data:
		ppl.axvspan(r['t_start'], r['t_stop'], color = BLOCK_COLOR_WG if r['Sample'] == 'WG' else BLOCK_COLOR, lw = 1)

	for r in data:
		for i in ISO_FIELDS[1:]:
			ppl.plot(r['t'][r['mask']], r['_d'+i][r['mask']], mfc = ISO_COLORS[i], **kwargs_valid)
			ppl.plot(r['t'][~r['mask']], r['_d'+i][~r['mask']], **kwargs_rejected)
			ppl.plot((r['t_start']+r['t_stop'])/2, r['d'+i], mec = ISO_COLORS[i], **kwargs_blockmean)

	ppl.sca(ax_residuals)
	ppl.ylabel('delta residuals')
	ppl.xlabel('t')

	for r in data:
		ppl.axvspan(r['t_start'], r['t_stop'], color = BLOCK_COLOR_WG if r['Sample'] == 'WG' else BLOCK_COLOR, lw = 1)

	for r in data:
		for i in ISO_FIELDS[1:]:
			ppl.plot(r['t'][r['mask']], r['_d'+i][r['mask']] - r['d'+i], mfc = ISO_COLORS[i], **kwargs_valid)
			ppl.plot(r['t'][~r['mask']], r['_d'+i][~r['mask']] - r['d'+i], **kwargs_rejected)

	ppl.axis([data[0]['t_start']-4, data[-1]['t_stop']+4, None, None])
	ppl.text(0.5, 0.97, run, transform = fig.transFigure, size = 18, weight = 'bold', va = 'top', ha = 'center')

	ppl.savefig('crunched.pdf')
	ppl.show()
	ppl.close(fig)


def save_to_csv(data, run, filename = 'crunched.csv'):
	with open(filename, 'w') as f:
		f.write('Run,Sample,P,' + ','.join([f'd{i}' for i in ISO_FIELDS[1:]]))
		for r in data:
			if r['Sample'] != 'WG':
				f.write(f'\n{run},{r["Sample"]},{r["P"][r["mask"]].mean():.6f},' + ','.join([f"{r['d'+i]:.6f}" for i in ISO_FIELDS[1:]]))


def save_config(filename = 'config.csv'):
	with open(filename, 'w') as f:
		f.write('Option,Value')
		f.write(f'\nWG_NAME,{WG_NAME}')
		f.write(f'\nGAP_WIDTH,{GAP_WIDTH}')
		f.write(f'\nIGNORE_FIRST,{IGNORE_FIRST}')
		f.write(f'\nSIGMA_FILTER,{SIGMA_FILTER}')

if __name__ == '__main__':

	data, run = read_input()
	data = cleanup(data)
	data = chop(data)
	data = crunch(data)
	data = dedrift(data)
	plot_all(data)
	save_to_csv(data, run)
	save_config()
