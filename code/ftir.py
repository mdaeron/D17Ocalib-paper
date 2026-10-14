#! /usr/bin/env uv run

import config
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as ppl

data = {}
for p in sorted(Path('inputs/ftir').glob('*.tsv')):
	df = pd.read_csv(p, sep = '\t', decimal=',').rename(columns = {'X': 'wn', 'Y': 'A'})
	data['calcite' if 'ISP' in p.stem else p.stem.replace('075', '75')] = df

kw = {
	'CC75': dict(lw = 1),
	'CC120': dict(lw = 1),
	'CC170': dict(lw = 1),
	'CC200': dict(lw = 1),
	'CC250': dict(lw = 1),
	'calcite': dict(dashes = (12,2), lw = 0.5, color = 'k'),
	'aragonite': dict(dashes = (2,2), lw = 0.5, color = 'k'),
}

fig = ppl.figure(figsize = (6,3))
ppl.subplots_adjust(0.025, 0.15, 0.975, 0.95)
ax = ppl.subplot(111)
for k,v in data.items():
	ax.plot(v.wn, v.A / v.A.to_numpy().max(), '-', label = k, **kw[k])
ax.legend()
ax.axis([1900, 500, 0, 1.05])
ax.set_yticks([])
ax.set_xlabel('Wavenumber (1/cm)')
fig.savefig('output/ftir.pdf')
