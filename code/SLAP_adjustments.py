#! /usr/bin/env uv run

import config
import tomli_w

refs = (
	'Barkan et al. (2019)',
	'Wostbrock et al. (2020)',
	'Ellis & Passey (2023)',
	'Chaillot et al. (2025)',
)

standards = (
	'IAEA-603',
	'NBS18',
	'IAEA-603-CO2',
	'NBS18-CO2',
	'DVH-C',
	'DVH-W',
)

D17O_VSMOW_OF_SLAP = (
	-0.015,
	-0.011,
	0.0,
	0.034,
)

toml_dict = {}
with open('output/SLAP_adjustments.csv', 'w') as f:
	f.write('D17O_VSMOW_of_SLAP_new,,,,,' + ','.join([f'{_:.3f}' for _ in D17O_VSMOW_OF_SLAP]))
	f.write('\n,Ref,D17O_VSMOW_of_SLAP_original,RM,D17O_VSMOW_original,')
	for ref in refs:
		key = ref
		toml_dict[key] = {}
		for standard in standards:
			try:
				old = config.SLAP_adjust(
					standard = standard,
					ref = ref,
					adjust = False,
				)
				Do = {
					'Barkan et al. (2019)': 0.000,
					'Wostbrock et al. (2020)': -0.015,
					'Ellis & Passey (2023)': -0.011,
					'Chaillot et al. (2025)': 0.000,
				}
				f.write(f'\n,{ref},{Do[ref]:.3f},{standard},{old:.3f}')
				for D in D17O_VSMOW_OF_SLAP:
					new = config.SLAP_adjust(
						standard = standard,
						ref = ref,
						new_D17O_VSMOW_OF_SLAP = D,
					)
					f.write(f',{new:.4f}')
					if D == config.D17O_VSMOW_OF_SLAP_POSTULATED_HERE:
						toml_dict[key][standard] = float(new)

			except UnboundLocalError:
				pass

with open('output/SLAP-results.toml', 'wb') as f:
	tomli_w.dump({'SLAPcorrected': toml_dict}, f)
