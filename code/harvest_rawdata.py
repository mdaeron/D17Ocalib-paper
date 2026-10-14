#! /usr/bin/env uv run

import json
import hashlib

from pathlib import Path
from config import vprint

def dict_hash(x):
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # We need to sort arguments so {'a': 1, 'b': 2} is
    # the same as {'b': 2, 'a': 1}
    encoded = json.dumps(x, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()

vprint(f'Starting {Path(__file__).name}')

out = []
for file in sorted(Path('inputs').glob('runs/2*/crunched.csv')):
	vprint(f'Adding run {file.parent.stem}')
	with open(file) as f:
		if out:
			out += f.readlines()[1:]
		else:
			out += f.readlines()

out = [l.strip().split(',') for l in out]
out = [{k:v for k,v in zip(out[0], l)} for l in out[1:]]

for r in out:
	try:
		t, r['Sample'], r['UID'], *extra = r['Sample'].split(':')
		r['Type'] = {'W': 'water', 'C': 'carbonate', 'F': 'fluorinated carbonate', 'G': 'carbon dioxide'}[t]
		if t == 'F':
			r['Sample'] = r['Sample'] + '_fluorinated'
		if extra:
			if t == 'C':
				if extra[0].endswith('C'):
					r['Tacid'] = extra[0][:-1]
			elif t == 'W':
				for _ in extra:
					if _.endswith('uL'):
						r['WaterVolume'] = _
					if _.endswith('umol'):
						r['AmountCO2'] = _
				v = float(r['WaterVolume'][:-2]) / 1e6
				n = float(r['AmountCO2'][:-4]) / 1e6
				r['H2O/CO2'] = f'{v / 18.015e-3 / n:.0f}'
				if r['Sample'].startswith('C13D'):
					r['eqCO2_d18O_VSMOW'] = '10.6'
					r['eqCO2_D17O_VSMOW'] = '-0.15'
				else:
					r['eqCO2_d18O_VSMOW'] = '24.24'
					r['eqCO2_D17O_VSMOW'] = '-0.15'
	except ValueError:
		pass

out = [r for r in out if 'PESNIN' not in r['Sample'] and 'GG' not in r['Sample']]
# out = [r for r in out if '2024-07' not in r['Run'] or r['Type'] == 'water']
out = [r for r in out if not (r['Run'] == '2024-07-10' and 'DRIFTING' in r['Sample'])]

for r in out:
	r['newUID'] = dict_hash(r)[:8].upper()
	r['Sample'] = r['Sample'].replace('JAUTZY_', '')
	r['Sample'] = r['Sample'].replace('-', '')
	if r['Run'] in ['2024-07-09', '2024-07-10'] and r['Type'] == 'carbonate':
		r['Run'] = (
			r['Run']
			+ {
				'240412': 'a',
				'240416': '',
				'240417': '',
			}[r['UID'][-7:-1]]
		)
	r['UID'] = r['newUID']

out = [r for r in out if not (r['Run'] == '2024-07-09a')]

keys = [
	'Run',
	'UID',
	'Sample',
	'Type',
	'P',
	'd636',
	'd628',
	'd627',
# 	'd638',
	'd626b',
	'Tacid',
	'WaterVolume',
	'AmountCO2',
	'H2O/CO2',
	'eqCO2_d18O_VSMOW',
	'eqCO2_D17O_VSMOW',
]
with open('inputs/rawdata.csv', 'w') as f:
	f.write(
		'\n'.join(
			[','.join(keys)]
			+ [','.join([r[k] if k in r else '' for k in keys]) for r in out]
			)
		)

vprint(f'Finished {Path(__file__).name}')
