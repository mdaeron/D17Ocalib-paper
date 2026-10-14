import toml
import numpy as np
import petroff_2021
from matplotlib import style
from matplotlib  import markers

style.use('code/mydefaults.mplstyle')

VERBOSE = True

if VERBOSE:
	vprint = print
else:
	def vprint(*args, **kwargs):
		return None

# ONE_COLUMN_WIDTH = 3.35
# TWO_COLUMN_WIDTH = 7
ONE_COLUMN_WIDTH = 85/25.4
TWO_COLUMN_WIDTH = 175/25.4

HIGH_T_LIMIT_OF_THETA = 0.5305

D17O_VSMOW_OF_SLAP_POSTULATED_HERE    = 0.000
P_GRUBBS                              = 0.95
RELATIVE_SIGMA_PRECISION              = 1e-1
FORCE_LGB_WATER_d18O_TO_LONG_TERM_AVERAGE = True
ONLY_CONSIDER_NATURAL_SAMPLES_FOR_REGRESSION = True

def SLAP_adjust(
	standard,
	ref,
	new_D17O_VSMOW_OF_SLAP = D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
	adjust = True,
):
	match ref:
		case 'Barkan et al. (2019)': # Barkan et al. (2019)
			old_D17O_VSMOW_OF_SLAP = 0.0
		case 'Wostbrock et al. (2020)': # Wostbrock et al. (2020), https://doi.org/10.1016/j.chemgeo.2019.119432
			old_D17O_VSMOW_OF_SLAP = -0.015
		case 'Ellis & Passey (2023)': # Ellis & Passey (2023), https://doi.org/10.1016/j.chemgeo.2023.121616
			old_D17O_VSMOW_OF_SLAP = -0.011
		case 'Chaillot et al. (2025)': # Chaillot et al. (2025), https://doi.org/10.31223/X53Q43
			old_D17O_VSMOW_OF_SLAP = 0.0

	match standard:
		case 'IAEA-603':
			d18O_VPDB = -2.37
			d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000
			match ref:
				case 'Wostbrock et al. (2020)':
					D17O_VSMOW_old = -0.100
				case 'Ellis & Passey (2023)':
					D17O_VSMOW_old = -0.082
		case 'IAEA-603-CO2':
			match ref:
				case 'Chaillot et al. (2025)':
					d18O_VPDB = -2.37
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.00813 - 1000 # CO2 from 90 °C acid
					D17O_VSMOW_old = -0.1273
				case 'Barkan et al. (2019)':
					d18O_VPDB = -2.37
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.01025 - 1000 # CO2 from 25 °C acid
					D17O_VSMOW_old = -0.194
				case 'Wostbrock et al. (2020)':
					d18O_VPDB = -2.37
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.01025 - 1000 # CO2 from 25 °C acid
					D17O_VSMOW_old = -0.147
		case 'NBS18':
			d18O_VPDB = -23.01
			d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000
			match ref:
				case 'Wostbrock et al. (2020)':
					D17O_VSMOW_old = -0.048
				case 'Ellis & Passey (2023)':
					D17O_VSMOW_old = -0.033
		case 'NBS18-CO2':
			match ref:
				case 'Chaillot et al. (2025)':
					d18O_VPDB = -23.01
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.00813 - 1000 # CO2 from 90 °C acid
					D17O_VSMOW_old = -0.1013
				case 'Barkan et al. (2019)':
					d18O_VPDB = -23.01
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.01025 - 1000 # CO2 from 25 °C acid
					D17O_VSMOW_old = -0.163
				case 'Wostbrock et al. (2020)':
					d18O_VPDB = -23.01
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.01025 - 1000 # CO2 from 25 °C acid
					D17O_VSMOW_old = -0.100
		case 'NBS19-CO2':
			match ref:
				case 'Chaillot et al. (2025)':
					d18O_VPDB = -2.20
					d18O_VSMOW = (1000 + d18O_VPDB) * 1.03092 - 1000  # carbonate
					d18O_VSMOW = (1000 + d18O_VSMOW) * 1.00813 - 1000 # CO2 from 90 °C acid
					D17O_VSMOW_old = -0.1304
		case 'DVH-C':
			d18O_VSMOW = 14.621
			match ref:
				case 'Wostbrock et al. (2020)':
					D17O_VSMOW_old = -0.043
		case 'DVH-W':
			d18O_VSMOW = -13.388
			match ref:
				case 'Wostbrock et al. (2020)':
					D17O_VSMOW_old = 0.025

	if adjust:
		return D17O_VSMOW_old + (
			(new_D17O_VSMOW_OF_SLAP - old_D17O_VSMOW_OF_SLAP)
			* np.log(1 + d18O_VSMOW/1e3) / np.log(1-55.5e-3)
		)
	else:
		return D17O_VSMOW_old

D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT = SLAP_adjust('IAEA-603-CO2', 'Chaillot et al. (2025)')
assert D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT == -0.1273

D17O_VSMOW_OF_NBS18_CO2_CHAILLOT = SLAP_adjust('NBS18-CO2', 'Chaillot et al. (2025)')
assert D17O_VSMOW_OF_NBS18_CO2_CHAILLOT == -0.1013

D17O_603_OF_NBS18_POSTULATED_HERE = D17O_VSMOW_OF_NBS18_CO2_CHAILLOT - D17O_VSMOW_OF_IAEA603_CO2_CHAILLOT

d18O_VSMOW_OF_DHW_WOSTBROCK = -13.388 # Wostbrock et al. (2020), https://doi.org/10.1016/j.gca.2020.07.045
d18O_VSMOW_OF_DHC_WOSTBROCK =  14.621 # Wostbrock et al. (2020), https://doi.org/10.1016/j.gca.2020.07.045
D17O_VSMOW_OF_DHW_WOSTBROCK =   0.025 # Wostbrock et al. (2020), https://doi.org/10.1016/j.gca.2020.07.045
D17O_VSMOW_OF_DHC_WOSTBROCK =  -0.043 # Wostbrock et al. (2020), https://doi.org/10.1016/j.gca.2020.07.045
klna18_OF_DVH_WOSTBROCK     =  27.994 # Wostbrock et al. (2020), https://doi.org/10.1016/j.gca.2020.07.045

d18O_VSMOW_OF_DHW_WAFFEK = -13.46 # Affek et al. (2025), https://doi.org/10.1016/j.epsl.2025.119528
D17O_VSMOW_OF_DHW_WAFFEK =  0.045  # Affek et al. (2025), https://doi.org/10.1016/j.epsl.2025.119528

CO2g_H2O_FRACTIONATION = (
	'Guo & Zhou (2019)'
# 	'Barkan & Luz (2012)'
)

_CO2g_H2O_fractionation = {
	'Guo & Zhou (2019)': dict(
		alpha18_CO2g_H2O = 1.041461382,
		theta17_CO2g_H2O = 0.524606996,
		),
	'Barkan & Luz (2012)': dict(
		alpha18_CO2g_H2O = 1.041036,
		theta17_CO2g_H2O = 0.5229,
		),
	}[CO2g_H2O_FRACTIONATION]

alpha18_CO2g_H2O = _CO2g_H2O_fractionation['alpha18_CO2g_H2O']
theta17_CO2g_H2O = _CO2g_H2O_fractionation['theta17_CO2g_H2O']
alpha17_CO2g_H2O =  np.exp(theta17_CO2g_H2O * np.log(alpha18_CO2g_H2O))

kw_errorbar = dict(
	ls = 'None',
	marker = 'o',
	ms = 5,
	capsize = 3,
	capthick = 1,
	elinewidth = 1,
	mfc = 'w',
	mew = 1,
)

# COLORS = {}
# for col, names in (
# 	('blue'      , ['Barkan', 'diseq', 'Daëron']),
# 	('brown'     , ['Roy']),
# 	('lightgray' , []),
# 	('gold'    , ['Coplen']),
# 	('orange'       , ['Wostbrock']),
# 	('purple'    , ['Ellis', 'Kelson', 'Huth']),
# 	('lightblue' , []),
# 	('red'      , ['ONeil']),
# 	('darkgray'  , ['Levitt', 'Gabitov']),
# 	('tan'       , ['Passey']),
# ):
# 	for name in names:
# 		COLORS[name] = petroff_2021.ten_color_sequence[col]

COLORS = {}
for col, names in (
	('lightblue' , ['Roy', 'Levitt', 'Gabitov']),
	('blue'      , ['Coplen', 'diseq', 'Barkan']),
	('darkblue'  , ['ONeil']),
	('purple'    , ['Ellis', 'Kelson', 'Huth']),
	('red'       , ['Daëron']),
	('orange'    , ['Wostbrock']),
	('pistacchio', ['Passey']),
	('gray'      , []),
):
	for name in names:
		COLORS[name] = np.array(petroff_2021.eight_color_sequence[col])

# COLORS = {}
# for col, names in (
# 	('blue'      , ['Roy', 'Coplen', 'diseq', 'Barkan']),
# 	('plum'      , ['ONeil']),
# 	('purple'    , ['Ellis', 'Kelson', 'Huth']),
# 	('red'       , ['Daëron']),
# 	('orange'    , ['Wostbrock']),
# 	('gray'      , ['Passey', 'Levitt', 'Gabitov']),
# ):
# 	for name in names:
# 		COLORS[name] = np.array(petroff_2021.six_color_sequence[col])

hdiamond = markers.MarkerStyle("d")
hdiamond._transform.rotate_deg(90)

ldiamond = markers.MarkerStyle("d")
ldiamond._transform.rotate_deg(59)

rdiamond = markers.MarkerStyle("d")
rdiamond._transform.rotate_deg(-59)

def D17O_adjust(
	d18O,
	D17O_old,
	d18O_RM1,
	d18O_RM2,
	D17O_RM1_old,
	D17O_RM2_old,
	D17O_RM1_new,
	D17O_RM2_new,
):
	'''
	General 2-pt renormalization of Δ17O measurements
	'''
	x = (
		np.log(1+d18O/1e3)
		- np.log(1+d18O_RM2/1e3)
	) / (
		np.log(1+d18O_RM1/1e3)
		- np.log(1+d18O_RM2/1e3)
	)
	change_in_D17O_RM1 = D17O_RM1_new - D17O_RM1_old
	change_in_D17O_RM2 = D17O_RM2_new - D17O_RM2_old

	return D17O_old + x * change_in_D17O_RM1 + (1-x) * change_in_D17O_RM2

def D17O_adjust_water(
	d18O,
	D17O_old,
	D17O_SLAP_old,
	D17O_SLAP_new = D17O_VSMOW_OF_SLAP_POSTULATED_HERE,
):
	'''
	2-pt renormalization of water Δ17O measurements
	'''
	return D17O_adjust(
		d18O         = d18O,
		D17O_old     = D17O_old,
		d18O_RM1     = 0.,
		d18O_RM2     = -55.5,
		D17O_RM1_old = 0.,
		D17O_RM2_old = D17O_SLAP_old,
		D17O_RM1_new = 0.,
		D17O_RM2_new = D17O_SLAP_new,
	)

def D17O_adjust_carb(
	d18O,
	D17O_old,
	D17O_IAEA603_old,
	D17O_NBS18_old,
	D17O_IAEA603_new,
	D17O_NBS18_new,
):
	'''
	2-pt renormalization of carbonate Δ17O measurements.
	!! d18O values relative to VSMOW
	'''
	return D17O_adjust(
		d18O         = d18O,
		D17O_old     = D17O_old,
		d18O_RM1     = (1000-2.37) * 1.03092 - 1000,
		d18O_RM2     = (1000-23.01) * 1.03092 - 1000,
		D17O_RM1_old = D17O_IAEA603_old,
		D17O_RM2_old = D17O_NBS18_old,
		D17O_RM1_new = D17O_IAEA603_new,
		D17O_RM2_new = D17O_NBS18_new,
	)

def df_to_toml(df, dest, mode = 'w', prefix = ''):
	with open(dest, mode) as f:
		if prefix:
			toml.dump({prefix: df.to_dict(orient = 'index')}, f)
		else:
			toml.dump(df.to_dict(orient = 'index'), f)
