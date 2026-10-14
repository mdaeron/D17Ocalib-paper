"""
Standardize stable isotope measurements by comparison with reference materials of known composition
"""

__author__    = 'Mathieu Daëron'
__contact__   = 'daeron@lsce.ipsl.fr'
__copyright__ = 'Copyright (c) 2024 Mathieu Daëron'
__license__   = 'MIT License - https://opensource.org/licenses/MIT'
__date__      = '2024-07-09'
__version__   = '0.0.1'


import numpy
import pandas
import lmfit
import scipy.stats
import scipy.linalg
from sklearn.covariance import MinCovDet
from matplotlib import pyplot
from matplotlib.patches import Ellipse
from scipy.optimize import fsolve
from collections import defaultdict
from scipy.stats import gaussian_kde


def estimate_covariance(X, robust = True, support_fraction = None, assume_centered = True):
	# X.shape must be (N,2)
	if robust:
		return MinCovDet(
			support_fraction = support_fraction,
			assume_centered = assume_centered
		).fit(X).covariance_
	else:
		return numpy.cov(X.T)

def sanitize(x):
	return x.replace('-', '').replace('.', '')


def cov_ellipse(CM, p = None, r = None):
	"""
	Parameters
	----------
	CM : (2, 2) array
		Covariance matrix.
	p : float
		Confidence level, should be in (0, 1)

	Returns
	-------
	width, height, rotation :
		 The lengths of two axes and the rotation angle in degree
	for the ellipse.
	"""

	assert r is None or p is None

	if r is None and p is None:
		p = 0.95
		r2 = scipy.stats.chi2.ppf(p, 2)
	if r is not None and p is None:
		r2 = r**2
	if p is not None and r is None:
		r2 = scipy.stats.chi2.ppf(p, 2)

	val, vec = scipy.linalg.eigh(CM)
	width, height = 2 * (val[:, None] * r2)**.5
	rotation = numpy.degrees(numpy.arctan2(*vec[::-1, 0]))

	return width, height, rotation


class Dataset:

	def __init__(self, data = [], index_field = 'UID'):
		"""
		Hold data to be standardized
		"""
		if isinstance(data, pandas.DataFrame):
			self.data = data.copy()
		elif isinstance(data, list):
			self.data = pandas.DataFrame(data)
		else:
			raise ValueError("Invalid data type. Allowed data types are pandas.DataFrame and list of dictionnaries.")
		if index_field:
			try:
				self.data.set_index(index_field, inplace = True)
			except KeyError:
				pass


	def read_csv(self, filepath_or_buffer, *args, **kwargs):
		"""
		Read additional data from csv file
		"""
		df = pandas.read_csv(filepath_or_buffer, *args, **kwargs)
		self.data = pandas.concat([self.data, df], ignore_index = True)


	def to_csv(self, filepath_or_buffer, *args, **kwargs):
		"""
		Write data to csv file
		"""
		self.data.to_csv(filepath_or_buffer, *args, **kwargs)


	def standardize(
		self,
		k_in,
		k_out,
		anchors,
		method = 'correct_observations',
		constraints = {},
	):
		"""
		Standardize data based on a set of anchors, i.e. a dict of the form:
		anchors = {
			'RefMaterial_1': RM1_value,
			'RefMaterial_2': RM2_value, ...
		}
		"""

		fitparams = lmfit.Parameters()

		for s in self.data.Session.unique():
			fitparams.add('delta_scaling_' + sanitize(s), value = 1.0)
			fitparams.add(  'delta_of_wg_' + sanitize(s), value = 0.0)
		for s in self.data.Sample.unique():
			if s not in anchors:
				fitparams.add('x_' + sanitize(s), value = 0.0)

		for p in fitparams:
			if p in constraints:
				fitparams[p].expr = constraints[p]

		def observations():
			return self.data[k_in]

		def truevalues(p):
			return self.data.Sample.map({
				s: anchors[s] if s in anchors else float(p[f"x_{sanitize(s)}"])
				for s in self.data.Sample.unique()
			})

		def correct_observations(p):
			delta_scalinp_map = {
				s: float(p[f"delta_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			delta_of_wg_map = {
				s: float(p[f"delta_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			return 1e3 * (
				(1 + self.data[k_in] / 1e3 / self.data.Session.map(delta_scalinp_map))
				* (1 + self.data.Session.map(delta_of_wg_map) / 1e3)
				- 1
			)

		def predict_observations(p):
			delta_scalinp_map = {
				s: float(p[f"delta_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			delta_of_wg_map = {
				s: float(p[f"delta_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			return (
				(
					(1 + truevalues(p) / 1e3)
					/ (1 + self.data.Session.map(delta_of_wg_map) / 1e3)
					- 1
				)
				* 1e3
				* self.data.Session.map(delta_scalinp_map)
			)

		def prediction_residuals(p):
			return (observations() - predict_observations(p)).array

		def correction_residuals(p):
			return (correct_observations(p) - truevalues(p)).array

		residuals = {
			'correct_observations': correction_residuals,
			'predict_observations': prediction_residuals,
		}[method]

		fitresult = lmfit.minimize(
			residuals, fitparams, method = 'least_squares', scale_covar = True
		)

		self.data[k_out+'_anchor'] = self.data['Sample'].isin(anchors)
		self.data[k_in+'_predicted'] = predict_observations(fitresult.params)
		self.data[k_in+'_residual'] = self.data[k_in] - self.data[k_in+'_predicted']
		self.data[k_out+'_corrected'] = correct_observations(fitresult.params)
		self.data[k_out+'_true'] = truevalues(fitresult.params)
		self.data[k_out+'_residual'] = self.data[k_out+'_corrected'] - self.data[k_out+'_true']

		if not hasattr(self, 'standardization'):
			self.standardization = dict()

		self.standardization[k_out] = dict(
			fitresult = fitresult,
			report = lmfit.fit_report(fitresult),
			Nf = fitresult.nfree,
			t95 = scipy.stats.t.ppf(1 - 0.05 / 2, fitresult.nfree)
		)

		if not hasattr(self, 'stdz_lookup'):
			self.stdz_lookup = dict()

		self.stdz_lookup[k_out] = k_out

		if not hasattr(self, 'sessions'):
			self.sessions = dict()

		self.sessions[k_out] = pandas.DataFrame()
		sessions = self.sessions[k_out]

		sessions['N'] = self.data.value_counts(subset='Session')
		NaNu = self.data.groupby(['Session', k_out+'_anchor']).size().to_dict()
		sessions['Na'] = sessions.index.map({
			s: NaNu[(s, True)] for s in sessions.index
		})
		sessions['Nu'] = sessions.index.map({
			s: NaNu[(s,  False)] for s in sessions.index
		})
		sessions[f'{k_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"delta_scaling_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{k_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"delta_scaling_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{k_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"delta_of_wg_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{k_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"delta_of_wg_{sanitize(s)}"].stderr
			for s in sessions.index
		})


		samples = pandas.DataFrame()
		samples['N'] = self.data.value_counts(subset='Sample')
		samples[k_in] = (
			self.data
			.groupby('Sample')
			.agg({k_in: 'mean'})
		)
		samples[k_out] = samples.index.map({
			s: anchors[s] if s in anchors else fitresult.params[f"x_{sanitize(s)}"].value
			for s in samples.index
		})
		samples['SE_'+k_out] = samples.index.map({
			s: None if s in anchors else fitresult.params[f"x_{sanitize(s)}"].stderr
			for s in samples.index
		})
		samples['95CL_'+k_out] = samples['SE_'+k_out] * self.standardization[k_out]['t95']
		samples['SD_'+k_out] = (
			self.data
			.groupby('Sample')
			.agg({k_out+'_corrected': 'std'})
		)

		if hasattr(self, 'samples'):
			columns = self.samples.columns.tolist()
			for c in samples.columns:
				if c not in columns:
					columns.append(c)
			self.samples = self.samples.combine_first(samples)[columns]
		else:
			self.samples = samples.copy()

		df_anchors = pandas.DataFrame({k_out: [anchors[s] for s in anchors]}, index = anchors.keys())

		if hasattr(self, 'anchors'):
			self.anchors = self.anchors.combine_first(df_anchors)
		else:
			self.anchors = df_anchors.copy()
			self.anchors.index.name = 'Sample'

		# TODO: compute anchor mismatch ie mean and SE of residuals


	def standardize_D17O(
		self,
		d18O_key_in,
		d17O_key_in,
		d18O_key_out,
		d17O_key_out,
		anchors,
		method = 'correct_observations', # 'correct_observations' | 'predict_observations'
		constraints = {},
		residuals_17 = 'D17', # 'D17' | 'd17'
		lambda_17 = 0.528,
		gamma_17 = 0., # permil
		robust_cov_estimator = True,
		relative_sigma_precision = 1e-5,
		max_iter = 10,
	):
		"""
		Standardize triple oxygen isotopes based on a set of anchors, i.e. a dict of the form:
		anchors = {
			'RefMaterial_1': {d18O_key_out: ..., D17O_key_out: ...},
			'RefMaterial_2': {d18O_key_out: ..., D17O_key_out: ...},
		}
		"""

		D17O_key_out = d17O_key_out.replace('d17', 'D17')

		anchors = anchors.copy()
		for s in anchors:
			if d17O_key_out not in anchors[s]:
				anchors[s][d17O_key_out] = 1e3 * (
					numpy.exp((anchors[s][D17O_key_out] + gamma_17) / 1e3)
					* (1 + anchors[s][d18O_key_out] / 1e3)**lambda_17
					- 1
				)
			if D17O_key_out not in anchors[s]:
				anchors[s][D17O_key_out] = 1e3 * (
					numpy.log(1+anchors[s][d17O_key_out] / 1e3)
					- lambda_17 * numpy.log(1+anchors[s][d18O_key_out] / 1e3)
				) - gamma_17

		# TODO: check consistency between d17O and D17O in anchors

		fitparams = lmfit.Parameters()

		for s in self.data.Session.unique():
			fitparams.add('d18O_scaling_' + sanitize(s), value = 1.0)
			fitparams.add(  'd18O_of_wg_' + sanitize(s), value = 0.0)
			fitparams.add('d17O_scaling_' + sanitize(s), value = 1.0)
			fitparams.add(  'd17O_of_wg_' + sanitize(s), value = 0.0)
		for s in self.data.Sample.unique():
			if s not in anchors:
				fitparams.add('d18O_' + sanitize(s), value = 0.0)
				fitparams.add('D17O_' + sanitize(s), value = 0.0)

		for p in fitparams:
			if p in constraints:
				fitparams[p].expr = constraints[p]

		def observations():
			return (
				self.data[d18O_key_in],
				self.data[d17O_key_in],
				1e3 * (
					numpy.log(1+self.data[d17O_key_in]/1e3)
					- lambda_17 * numpy.log(1+self.data[d18O_key_in]/1e3)
				) - gamma_17,
			)

		def truevalues(p):

			d18true = self.data.Sample.map({
				s: anchors[s][d18O_key_out]
				if s in anchors
				else float(p[f"d18O_{sanitize(s)}"])
				for s in self.data.Sample.unique()
			})

			D17true = self.data.Sample.map({
				s: anchors[s][D17O_key_out]
				if s in anchors
				else float(p[f"D17O_{sanitize(s)}"])
				for s in self.data.Sample.unique()
			})

			d17true = 1e3 * (
				numpy.exp((D17true + gamma_17)/1e3)
				* (1+d18true/1e3)**lambda_17
				- 1
			)

			return (d18true, d17true, D17true)

		def correct_observations(p):
			d18O_scalinp_map = {
				s: float(p[f"d18O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18O_of_wg_map = {
				s: float(p[f"d18O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18corrected = 1e3 * (
				(1 + self.data[d18O_key_in] / 1e3 / self.data.Session.map(d18O_scalinp_map))
				* (1 + self.data.Session.map(d18O_of_wg_map) / 1e3)
				- 1
			)

			d17O_scalinp_map = {
				s: float(p[f"d17O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17O_of_wg_map = {
				s: float(p[f"d17O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17corrected = 1e3 * (
				(1 + self.data[d17O_key_in] / 1e3 / self.data.Session.map(d17O_scalinp_map))
				* (1 + self.data.Session.map(d17O_of_wg_map) / 1e3)
				- 1
			)

			D17corrected = 1e3 * (
				numpy.log(1+d17corrected/1e3)
				- lambda_17 * numpy.log(1+d18corrected/1e3)
			) - gamma_17

			return (d18corrected, d17corrected, D17corrected)

		def predict_observations(p):
			d18true, d17true, D17true = truevalues(p)

			d18O_scalinp_map = {
				s: float(p[f"d18O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18O_of_wg_map = {
				s: float(p[f"d18O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18predicted = (
				(
					(1 + d18true / 1e3)
					/ (1 + self.data.Session.map(d18O_of_wg_map) / 1e3)
					- 1
				)
				* 1e3
				* self.data.Session.map(d18O_scalinp_map)
			)

			d17O_scalinp_map = {
				s: float(p[f"d17O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17O_of_wg_map = {
				s: float(p[f"d17O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17predicted = (
				(
					(1 + d17true / 1e3)
					/ (1 + self.data.Session.map(d17O_of_wg_map) / 1e3)
					- 1
				)
				* 1e3
				* self.data.Session.map(d17O_scalinp_map)
			)

			D17predicted = 1e3 * (
				numpy.log(1+d17predicted/1e3)
				- lambda_17 * numpy.log(1+d18predicted/1e3)
			) - gamma_17

			return (d18predicted, d17predicted, D17predicted)

		def prediction_residuals(p, residuals_17 = residuals_17, ChM = None):

			d18predicted, d17predicted, D17predicted = predict_observations(p)
			d18obs, d17obs, D17obs = observations()

			if residuals_17 == 'D17':
				R = numpy.array([
					d18predicted - d18obs,
					D17predicted - D17obs,
				])

			elif residuals_17 == 'd17':
				R = numpy.array([
					d18predicted - d18obs,
					d17predicted - d17obs,
				])

			if ChM is None:
				return R.flatten()
			else:
				return (ChM @ R).flatten()

		def correction_residuals(p, residuals_17 = residuals_17, ChM = None):

			d18true, d17true, D17true = truevalues(p)
			d18corrected, d17corrected, D17corrected = correct_observations(p)

			if residuals_17 == 'D17':
				R = numpy.array([
					d18corrected - d18true,
					D17corrected - D17true,
				])

			elif residuals_17 == 'd17':
				R = numpy.array([
					d18corrected - d18true,
					d17corrected - d17true,
				])
#
			if ChM is None:
				return R.flatten()
			else:
				return (ChM @ R).flatten()

		residuals = {
			'correct_observations': correction_residuals,
			'predict_observations': prediction_residuals,
		}[method]

		fitresult = lmfit.minimize(
			residuals, fitparams, method = 'least_squares', scale_covar = True
		)

		_s1, _s2, _r = 0, 0, 10

		for _ in range(max_iter):
			CM = estimate_covariance(
				residuals(fitresult.params).reshape((2, self.data.shape[0])).T,
				robust = robust_cov_estimator,
				)
			if (
				abs(_s1/CM[0,0]**.5 - 1) < relative_sigma_precision
				and abs(_s2/CM[1,1]**.5 - 1) < relative_sigma_precision
				and abs(_r - CM[0,1]/CM[0,0]**.5/CM[1,1]**.5) < relative_sigma_precision
			):
				break

			_s1, _s2, _r = CM[0,0]**.5, CM[1,1]**.5, CM[0,1]/CM[0,0]**.5/CM[1,1]**.5

			ChM = scipy.linalg.cholesky(scipy.linalg.inv(CM))
			fitresult = lmfit.minimize(
				residuals,
				fitparams,
				method = 'least_squares',
				scale_covar = True, # ensures that we are insensitive to CM scaling
				kws = {'ChM': ChM},
				)
		else:
			raise ValueError(f'D17O Standardization failed to converge after {_+1} iterations.')

		self.data[D17O_key_out+'_anchor'] = self.data['Sample'].isin(anchors)

		self.data[d18O_key_out+'_corrected'] = correct_observations(fitresult.params)[0]
		self.data[d18O_key_out+'_true'] = truevalues(fitresult.params)[0]
		self.data[d18O_key_out+'_residual'] = self.data[d18O_key_out+'_corrected'] - self.data[d18O_key_out+'_true']

		self.data[D17O_key_out+'_corrected'] = correct_observations(fitresult.params)[2]
		self.data[D17O_key_out+'_true'] = truevalues(fitresult.params)[2]
		self.data[D17O_key_out+'_residual'] = self.data[D17O_key_out+'_corrected'] - self.data[D17O_key_out+'_true']

		if not hasattr(self, 'standardization'):
			self.standardization = dict()

		self.standardization[D17O_key_out] = dict(
			fitresult = fitresult,
			report = lmfit.fit_report(fitresult),
			rchisq = fitresult.redchi,
			Nf = fitresult.nfree // 2,
			t95 = scipy.stats.t.ppf(1 - 0.05 / 2, fitresult.nfree // 2)
		)

		if not hasattr(self, 'stdz_lookup'):
			self.stdz_lookup = dict()

		self.stdz_lookup[d18O_key_out] = D17O_key_out
		self.stdz_lookup[D17O_key_out] = D17O_key_out

		if not hasattr(self, 'sessions'):
			self.sessions = dict()

		self.sessions[D17O_key_out] = pandas.DataFrame()
		sessions = self.sessions[D17O_key_out]

		sessions['N'] = self.data.value_counts(subset='Session')
		NaNu = self.data.groupby(['Session', D17O_key_out+'_anchor']).size().to_dict()
		sessions['Na'] = sessions.index.map({
			s: NaNu[(s, True)] for s in sessions.index
		})
		sessions['Nu'] = sessions.index.map({
			s: NaNu[(s,  False)] for s in sessions.index
		})
		sessions[f'{d18O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d18O_scaling_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d18O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d18O_scaling_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d18O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d18O_of_wg_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d18O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d18O_of_wg_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d17O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d17O_scaling_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d17O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d17O_scaling_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d17O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d17O_of_wg_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d17O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d17O_of_wg_{sanitize(s)}"].stderr
			for s in sessions.index
		})


		samples = pandas.DataFrame()
		samples['N'] = self.data.value_counts(subset='Sample')

		samples[d18O_key_in] = (
			self.data
			.groupby('Sample')
			.agg({d18O_key_in: 'mean'})
		)

		samples[d17O_key_in] = (
			self.data
			.groupby('Sample')
			.agg({d17O_key_in: 'mean'})
		)

		samples[d18O_key_out] = samples.index.map({
			s: anchors[s][d18O_key_out] if s in anchors else fitresult.params[f"d18O_{sanitize(s)}"].value
			for s in samples.index
		})
		samples['SE_'+d18O_key_out] = samples.index.map({
			s: None if s in anchors else fitresult.params[f"d18O_{sanitize(s)}"].stderr
			for s in samples.index
		})
		samples['95CL_'+d18O_key_out] = samples['SE_'+d18O_key_out] * self.standardization[D17O_key_out]['t95']
		samples['SD_'+d18O_key_out] = (
			self.data
			.groupby('Sample')
			.agg({d18O_key_out+'_corrected': 'std'})
		)

		samples[d17O_key_out] = 0.

		samples[D17O_key_out] = samples.index.map({
			s: anchors[s][D17O_key_out] if s in anchors else fitresult.params[f"D17O_{sanitize(s)}"].value
			for s in samples.index
		})
		samples['SE_'+D17O_key_out] = samples.index.map({
			s: None if s in anchors else fitresult.params[f"D17O_{sanitize(s)}"].stderr
			for s in samples.index
		})
		samples['95CL_'+D17O_key_out] = samples['SE_'+D17O_key_out] * self.standardization[D17O_key_out]['t95']
		samples['SD_'+D17O_key_out] = (
			self.data
			.groupby('Sample')
			.agg({D17O_key_out+'_corrected': 'std'})
		)

		samples[d17O_key_out] = 1e3 * (
				numpy.exp((samples[D17O_key_out] + gamma_17) / 1e3)
				* (1 + samples[d18O_key_out] / 1e3)**lambda_17
				- 1
			)

		# TODO: propagate uncertainties to d17O_key_out using the code below
		#
		# try:
		# 	i = fitresult.var_names.index('d18O_'+_s)
		# 	j = fitresult.var_names.index('D17O_'+_s)
		# 	CM = fitresult.covar[ix_([i,j],[i,j])]
		#
		# 	J = array([
		# 		[
		# 			exp(out['samples'][s][D17O_key_out]/1e3)
		# 			* ISOPARAMETERS.lambda_17
		# 			* (1+out['samples'][s][d18O_key_out]/1e3) ** (ISOPARAMETERS.lambda_17 - 1)
		# 			],
		# 		[
		# 			exp(out['samples'][s][D17O_key_out]/1e3)
		# 			* (1+out['samples'][s][d18O_key_out]/1e3) ** ISOPARAMETERS.lambda_17
		# 			],
		# 		])
		#
		# 	out['samples'][s]['SE_'+d17O_key_out] = float(J.T @ CM @ J)**.5
		# 	out['samples'][s]['95CL_'+d17O_key_out] = out['samples'][s]['SE_'+d17O_key_out] * out['t95']
		# except ValueError:
		# 	pass

		if hasattr(self, 'samples'):
			columns = self.samples.columns.tolist()
			for c in samples.columns:
				if c not in columns:
					columns.append(c)
			self.samples = pandas.concat([self.samples, samples], axis = 'columns')
		else:
			self.samples = samples.copy()


		df_anchors = pandas.DataFrame([{
			d18O_key_out: anchors[s][d18O_key_out],
			D17O_key_out: anchors[s][D17O_key_out],
			} for s in anchors], index = anchors.keys()
		)

		if hasattr(self, 'anchors'):
# 			self.anchors = self.anchors.combine_first(df_anchors)
			self.anchors = pandas.concat([self.anchors, df_anchors], axis = 'columns')
		else:
			self.anchors = df_anchors.copy()
			self.anchors.index.name = 'Sample'

		# TODO: compute anchor mismatch ie mean and SE of residuals

	def plot_residuals(self,
		key,
		ylabel = None,
		subplots = (221, 325, (2,20,20)),
		subplots_adjust = [.1, .1, .98, .9, -0.75, -0.25],
		min_width_for_horizontal_session_labels = 0.1,
	):

		anchor_key = self.stdz_lookup[key] + '_anchor'

		ax1 = pyplot.subplot(*(subplots[0] if hasattr(subplots[0], '__iter__') else (subplots[0],)))
		ax2 = pyplot.subplot(*(subplots[1] if hasattr(subplots[1], '__iter__') else (subplots[1],)), sharex = ax1)
		ax3 = pyplot.subplot(*(subplots[2] if hasattr(subplots[2], '__iter__') else (subplots[2],)), sharey = ax1)

		pyplot.subplots_adjust(*subplots_adjust)

		X = numpy.arange(len(self.data)) + 0.5
		_s = ''
		sessions = []
		for x, s in zip(X, self.data['Session']):
			if s != _s:
				sessions.append(x - 0.5)
				sessions.append(s)
				_s = s
		sessions.append(x + 0.5)

		for x1, x2, s in zip(sessions[:-2:2], sessions[2::2], sessions[1::2]):
			if x1:
				ax1.axvline(x1, color = [.75]*3, lw = 0.8)
				ax2.axvline(x1, color = [.75]*3, lw = 0.8)
			ax1.text(
				(x1+x2)/2,
				1.02,
				s,
				transform = ax1.get_xaxis_transform(),
				**(
					dict(
						va = 'bottom',
						ha = 'center',
					)
					if ((x2-x1)/len(self.data)) > min_width_for_horizontal_session_labels
					else dict(
						va = 'bottom',
						ha = 'left',
						rotation = 30,
					)
				)
			)

		ax1.plot(X, self.data.where( self.data[anchor_key])[key + '_residual'], 'ro', mew = 0, ms = 5)
		ax1.plot(X, self.data.where(~self.data[anchor_key])[key + '_residual'], 'ko', mew = 0, ms = 5)

		samples = self.data.sort_values([anchor_key, 'Sample'], ascending=[False, True])['Sample'].unique()
		y_samples = {s:k+1 for k,s in enumerate(samples)} | {None: None}

		ax2.plot(X, self.data.where( self.data[anchor_key])['Sample'].map(y_samples), 'ro', mew = 0, ms = 5)
		ax2.plot(X, self.data.where(~self.data[anchor_key])['Sample'].map(y_samples), 'ko', mew = 0, ms = 5)


		ax2.axis([0, len(self.data), 0, len(samples)+1])

		ax2.set_yticks(numpy.arange(len(samples))+1)
		ax2.set_xticks(X)
		ax2.xaxis.set_ticklabels(
			self.data.index,
			rotation = -45,
			ha = 'left',
			va = 'top',
			size = 5,
			)
		ax2.yaxis.set_ticklabels(
			samples,
			size = 9,
			)
		ax1.tick_params(
			top = False,
			labeltop = False,
			bottom = False,
			labelbottom = False,
		)
		ax2.tick_params(
			left = False,
			labelleft = False,
			right = False,
			labelright = True,
		)

		ymin, ymax = ax1.axis()[-2:]
		_ = (ymax - ymin)/3
		ymin -= _
		ymax += _

		yi = numpy.linspace(ymin, ymax, 201)
		xi = gaussian_kde(self.data[key + '_residual']).evaluate(yi)
		ax3.fill_betweenx(yi, xi, xi*0, fc = (0,0,0,.15), lw = 1, ec = (.75,.75,.75))
		ax3.axis([0, None, ymin, ymax])

		SD = (
			(self.data[key + '_residual']**2).sum()
			/ (
				len(self.data)
				- 2 * len(self.data.Session.unique())
				- len(self.data[~self.data[anchor_key]]['Sample'].unique())
			)
		)**.5

		ax3.text(0, (ymin+ymax)/2, f'   SD = {SD*1000:.1f} ppm', va = 'center', size = 10)

		ax3.spines['right'].set_visible(False)
		ax3.spines['top'].set_visible(False)
		ax3.spines['bottom'].set_visible(False)

		ax3.tick_params(
			left = False,
			labelleft = False,
			bottom = False,
			labelbottom = False,
		)
		ax3.grid(False)

		if ylabel is None:
			 ylabel = key + ' residuals   [‰]'
		ax1.set_ylabel(ylabel)

		return ax1, ax2

	def plot_D17O_residuals(self,
		d18O_key_out = 'd18O_VSMOW',
		D17O_key_out = 'D17O_VSMOW',
		margins = 0.25,
		cbins = [-1, 0, 0.25, 0.5, 0.75, 0.95],
		robust_cov_estimator = True,
		kw_plot = dict(ls = 'None', marker = '+', mec = 'k', alpha = 0.5, mew = 0.5, ms = 4),
		kw_ellipse = dict(ec = 'k', lw = 1, ls = (0, (6, 2, 2, 2))),
	):
		rchisq_key_out = self.stdz_lookup[D17O_key_out]

		x = self.data[d18O_key_out+'_residual'].array
		y = self.data[D17O_key_out+'_residual'].array

		xmin, xmax = x.min(), x.max()
		xmin, xmax = (
			xmin - (xmax - xmin) * margins,
			xmax + (xmax - xmin) * margins,
		)

		ymin, ymax = y.min(), y.max()
		ymin, ymax = (
			ymin - (ymax - ymin) * margins,
			ymax + (ymax - ymin) * margins,
		)

		xx, yy = numpy.mgrid[xmin:xmax:1000j, ymin:ymax:1000j]
		positions = numpy.vstack([xx.ravel(), yy.ravel()])
		values = numpy.vstack([x, y])
		kernel = scipy.stats.gaussian_kde(values)

		f = numpy.reshape(kernel(positions).T, xx.shape)
		f /= f.max()
		f = 1 - f

		ax = pyplot.gca()

		ax.set_xlim(xmin, xmax)
		ax.set_ylim(ymin, ymax)

		cfset = ax.contourf(xx, yy, f, levels = cbins, cmap = 'Blues_r', alpha = 1)
		cset = ax.contour(xx, yy, f, levels = cbins, colors = 'k')
		ax.plot(x, y, **kw_plot)

		ax.clabel(cset, inline = 1, fontsize = 10)
		ax.set_xlabel('δ$^{18}$O   [‰]')
		ax.set_ylabel('Δ$^{17}$O   [‰]')

		## Or kernel density estimate plot instead of the contourf plot
		# ax.imshow(numpy.rot90(f), cmap = 'Blues', extent = [xmin, xmax, ymin, ymax])		# Contour plot

		CM = (
			estimate_covariance(numpy.array([x, y]).T, robust = robust_cov_estimator)
			* self.standardization[rchisq_key_out]['rchisq']
		)

		w, h, r = cov_ellipse(CM)
		ax.add_patch(
			Ellipse(
				xy = (0, 0),
				width = w.item(),
				height = h.item(),
				angle = r.item(),
				fc = 'None',
				**kw_ellipse,
			)
		)
		ax.plot(
			[], [],
			label = '95 % confidence from\ncovariance estimate',
			color = kw_ellipse['ec'],
			lw = kw_ellipse['lw'],
			ls = kw_ellipse['ls'],
		)

		ax.legend(fontsize = 7, handlelength = 2.5)

		return ax



	def standardize_water_eq_CO2_D17O(
		self,
		water_anchors,
		d18O_key_in = 'd628',
		d17O_key_in = 'd627',
		d18O_key_out = 'd18O_VSMOW',
		d17O_key_out = 'd17O_VSMOW',
		alpha18_CO2g_H2O = 1.041036, # Barkan & Luz (2012)
		theta17_CO2g_H2O = 0.5229,   # Barkan & Luz (2012)
		R13_VPDB = 0.01118,
		R18_VSMOW = 0.00200520,
		R17_VSMOW = 0.00038475,
		lambda_17 = 0.528,
		method = 'correct_observations', # 'correct_observations' | 'predict_observations'
		residuals_17 = 'D17', # 'D17' | 'd17'
		constraints = {},
		robust_cov_estimator = True,
		relative_sigma_precision = 1e-3,
		max_iter = 10,
	):

		D17O_key_out = d17O_key_out.replace('d17', 'D17')

		water_anchors = water_anchors.copy()
		for s in water_anchors:
			if d17O_key_out not in water_anchors[s]:
				water_anchors[s][d17O_key_out] = 1e3 * (
					numpy.exp(water_anchors[s][D17O_key_out] / 1e3)
					* (1 + water_anchors[s][d18O_key_out] / 1e3)**lambda_17
					- 1
				)
			if D17O_key_out not in water_anchors[s]:
				water_anchors[s][D17O_key_out] = 1e3 * (
					numpy.log(1+water_anchors[s][d17O_key_out] / 1e3)
					- lambda_17 * numpy.log(1+water_anchors[s][d18O_key_out] / 1e3)
				)

		# TODO: check consistency between d17O and D17O in anchors

		fitparams = lmfit.Parameters()

		for s in self.data.Session.unique():
			fitparams.add('d18O_scaling_' + sanitize(s), value = 1.0)
			fitparams.add(  'd18O_of_wg_' + sanitize(s), value = 0.0)
			fitparams.add('d17O_scaling_' + sanitize(s), value = 1.0)
			fitparams.add(  'd17O_of_wg_' + sanitize(s), value = 0.0)
		for s in self.data.Sample.unique():
			if s not in water_anchors:
				fitparams.add('d18O_' + sanitize(s), value = 0.0)
				fitparams.add('D17O_' + sanitize(s), value = 0.0)

		for p in fitparams:
			if p in constraints:
				fitparams[p].expr = constraints[p]

		def observations():
			return (
				self.data[d18O_key_in],
				self.data[d17O_key_in],
				1e3 * (
					numpy.log(1+self.data[d17O_key_in]/1e3)
					- lambda_17 * numpy.log(1+self.data[d18O_key_in]/1e3)
				),
			)

		def _CO2eq_fun(r):
			return D17_CO2_H2O_eq(
				h2o_to_co2_molecular_ratio = float(r['H2O/CO2']),
				d18i_co2 = float(r['eqCO2_d18O_VSMOW']),
				D17i_co2 = float(r['eqCO2_D17O_VSMOW']),
				d18i_h2o = float(r[f'water_{d18O_key_out}_true']),
				D17i_h2o = float(r[f'water_{D17O_key_out}_true']),
				d13_co2 = float(r[f'eqCO2_d13C_VPDB']),
				alpha18_CO2g_H2O = alpha18_CO2g_H2O,
				alpha17_CO2g_H2O = alpha18_CO2g_H2O ** theta17_CO2g_H2O,
				R13_VPDB = R13_VPDB,
				R18_VSMOW = R18_VSMOW,
				R17_VSMOW = R17_VSMOW,
				lambda_17 = lambda_17,
				D727i = 0.1043,
				D728i = 0.1951,
				D828i = 0.3743,
				D637i = 0.5101,
				D638i = 0.9714,
				D737i = 1.1188,
				D738i = 1.6662,
				D838i = 2.3046,
				D727=0.1043,
				D728=0.1951,
				D828=0.3743,
				D637=0.5101,
				D638=0.9714,
				D737=1.1188,
				D738=1.6662,
				D838=2.3046,
				stochastic = False,
				xtol = 1e-12,
			)[:2]


		def truevalues(p):

			self.data[f'water_{d18O_key_out}_true'] = self.data.Sample.map({
				s: water_anchors[s][d18O_key_out]
				if s in water_anchors
				else float(p[f"d18O_{sanitize(s)}"])
				for s in self.data.Sample.unique()
			})

			self.data[f'water_{D17O_key_out}_true'] = self.data.Sample.map({
				s: water_anchors[s][D17O_key_out]
				if s in water_anchors
				else float(p[f"D17O_{sanitize(s)}"])
				for s in self.data.Sample.unique()
			})

			_ = self.data.apply(_CO2eq_fun, axis = 1, result_type = 'expand')

			d18co2, D17co2 = _[0], _[1]

			d17co2 = 1e3 * (
				numpy.exp((D17co2)/1e3)
				* (1+d18co2/1e3)**lambda_17
				- 1
			)

			return (d18co2, d17co2, D17co2)

		def correct_observations(p):
			d18O_scalinp_map = {
				s: float(p[f"d18O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18O_of_wg_map = {
				s: float(p[f"d18O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18corrected = 1e3 * (
				(1 + self.data[d18O_key_in] / 1e3 / self.data.Session.map(d18O_scalinp_map))
				* (1 + self.data.Session.map(d18O_of_wg_map) / 1e3)
				- 1
			)

			d17O_scalinp_map = {
				s: float(p[f"d17O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17O_of_wg_map = {
				s: float(p[f"d17O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17corrected = 1e3 * (
				(1 + self.data[d17O_key_in] / 1e3 / self.data.Session.map(d17O_scalinp_map))
				* (1 + self.data.Session.map(d17O_of_wg_map) / 1e3)
				- 1
			)

			D17corrected = 1e3 * (
				numpy.log(1+d17corrected/1e3)
				- lambda_17 * numpy.log(1+d18corrected/1e3)
			)

			return (d18corrected, d17corrected, D17corrected)

		def predict_observations(p):
			d18true, d17true, D17true = truevalues(p)

			d18O_scalinp_map = {
				s: float(p[f"d18O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18O_of_wg_map = {
				s: float(p[f"d18O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d18predicted = (
				(
					(1 + d18true / 1e3)
					/ (1 + self.data.Session.map(d18O_of_wg_map) / 1e3)
					- 1
				)
				* 1e3
				* self.data.Session.map(d18O_scalinp_map)
			)

			d17O_scalinp_map = {
				s: float(p[f"d17O_scaling_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17O_of_wg_map = {
				s: float(p[f"d17O_of_wg_{sanitize(s)}"])
				for s in self.data.Session.unique()
				}
			d17predicted = (
				(
					(1 + d17true / 1e3)
					/ (1 + self.data.Session.map(d17O_of_wg_map) / 1e3)
					- 1
				)
				* 1e3
				* self.data.Session.map(d17O_scalinp_map)
			)

			D17predicted = 1e3 * (
				numpy.log(1+d17predicted/1e3)
				- lambda_17 * numpy.log(1+d18predicted/1e3)
			)

			return (d18predicted, d17predicted, D17predicted)

		def prediction_residuals(p, residuals_17 = residuals_17, ChM = None):

			d18predicted, d17predicted, D17predicted = predict_observations(p)
			d18obs, d17obs, D17obs = observations()

			if residuals_17 == 'D17':
				R = numpy.array([
					d18predicted - d18obs,
					D17predicted - D17obs,
				])

			elif residuals_17 == 'd17':
				R = numpy.array([
					d18predicted - d18obs,
					d17predicted - d17obs,
				])

			if ChM is None:
				return R.flatten()
			else:
				return (ChM @ R).flatten()

		def correction_residuals(p, residuals_17 = residuals_17, ChM = None):

			d18true, d17true, D17true = truevalues(p)
			d18corrected, d17corrected, D17corrected = correct_observations(p)

			if residuals_17 == 'D17':
				R = numpy.array([
					d18corrected - d18true,
					D17corrected - D17true,
				])

			elif residuals_17 == 'd17':
				R = numpy.array([
					d18corrected - d18true,
					d17corrected - d17true,
				])
#
			if ChM is None:
				return R.flatten()
			else:
				_ = numpy.linalg.norm(ChM @ R, axis = 0)
				return (ChM @ R).flatten()

		residuals = {
			'correct_observations': correction_residuals,
			'predict_observations': prediction_residuals,
		}[method]

		fitresult = lmfit.minimize(
			residuals, fitparams, method = 'least_squares', scale_covar = True
		)

		_s1, _s2, _r = 0, 0, 10

		for _ in range(max_iter):
			CM = estimate_covariance(
				residuals(fitresult.params).reshape((2, self.data.shape[0])).T,
				robust = robust_cov_estimator,
				)

# 			print(f'{CM[0,0]**.5:.4f}, {CM[1,1]**.5:.4f}, {CM[0,1]/CM[0,0]**.5/CM[1,1]**.5:.3f}')

			if (
				abs(_s1/CM[0,0]**.5 - 1) < relative_sigma_precision
				and abs(_s2/CM[1,1]**.5 - 1) < relative_sigma_precision
				and abs(_r - CM[0,1]/CM[0,0]**.5/CM[1,1]**.5) < relative_sigma_precision
			):
				break

			_s1, _s2, _r = CM[0,0]**.5, CM[1,1]**.5, CM[0,1]/CM[0,0]**.5/CM[1,1]**.5

			ChM = scipy.linalg.cholesky(scipy.linalg.inv(CM))
			fitresult = lmfit.minimize(
				residuals,
				fitparams,
				method = 'least_squares',
				scale_covar = True, # ensures that we are insensitive to CM scaling
				kws = {'ChM': ChM},
				)
		else:
			raise ValueError(f'D17O Standardization failed to converge after {_+1} iterations.')

		self.data['water_'+D17O_key_out+'_anchor'] = self.data['Sample'].isin(water_anchors)

		self.data['CO2_'+d18O_key_out+'_corrected'] = correct_observations(fitresult.params)[0]
		self.data['CO2_'+d18O_key_out+'_true'] = truevalues(fitresult.params)[0]
		self.data['CO2_'+d18O_key_out+'_residual'] = self.data['CO2_'+d18O_key_out+'_corrected'] - self.data['CO2_'+d18O_key_out+'_true']

		self.data['CO2_'+D17O_key_out+'_corrected'] = correct_observations(fitresult.params)[2]
		self.data['CO2_'+D17O_key_out+'_true'] = truevalues(fitresult.params)[2]
		self.data['CO2_'+D17O_key_out+'_residual'] = self.data['CO2_'+D17O_key_out+'_corrected'] - self.data['CO2_'+D17O_key_out+'_true']


		if not hasattr(self, 'standardization'):
			self.standardization = dict()

		self.standardization['water_'+D17O_key_out] = dict(
			fitresult = fitresult,
			report = lmfit.fit_report(fitresult),
			rchisq = fitresult.redchi,
			Nf = fitresult.nfree // 2,
			t95 = scipy.stats.t.ppf(1 - 0.05 / 2, fitresult.nfree // 2)
		)

		if not hasattr(self, 'stdz_lookup'):
			self.stdz_lookup = dict()

		self.stdz_lookup['CO2_'+d18O_key_out]   = 'water_'+D17O_key_out
		self.stdz_lookup['CO2_'+D17O_key_out]   = 'water_'+D17O_key_out
		self.stdz_lookup['water_'+d18O_key_out] = 'water_'+D17O_key_out
		self.stdz_lookup['water_'+D17O_key_out] = 'water_'+D17O_key_out

		if not hasattr(self, 'sessions'):
			self.sessions = dict()

		self.sessions['water_'+D17O_key_out] = pandas.DataFrame()
		sessions = self.sessions['water_'+D17O_key_out]

		sessions['N'] = self.data.value_counts(subset='Session')

		NaNu = defaultdict(
			lambda : 0,
			self.data.groupby(['Session', 'water_'+D17O_key_out+'_anchor']).size().to_dict()
			)
		sessions['Na'] = sessions.index.map({
			s: NaNu[(s, True)] for s in sessions.index
		})
		sessions['Nu'] = sessions.index.map({
			s: NaNu[(s,  False)] for s in sessions.index
		})

		sessions[f'{d18O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d18O_scaling_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d18O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d18O_scaling_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d18O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d18O_of_wg_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d18O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d18O_of_wg_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d17O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d17O_scaling_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d17O_key_in}_scaling'] = sessions.index.map({
			s: fitresult.params[f"d17O_scaling_{sanitize(s)}"].stderr
			for s in sessions.index
		})
		sessions[f'{d17O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d17O_of_wg_{sanitize(s)}"].value
			for s in sessions.index
		})
		sessions[f'SE_{d17O_key_out}_of_wg'] = sessions.index.map({
			s: fitresult.params[f"d17O_of_wg_{sanitize(s)}"].stderr
			for s in sessions.index
		})


		samples = pandas.DataFrame()
		samples['N'] = self.data.value_counts(subset='Sample')

		samples[d18O_key_in] = (
			self.data
			.groupby('Sample')
			.agg({d18O_key_in: 'mean'})
		)

		samples[d17O_key_in] = (
			self.data
			.groupby('Sample')
			.agg({d17O_key_in: 'mean'})
		)

		samples['water_'+d18O_key_out] = samples.index.map({
			s: water_anchors[s][d18O_key_out] if s in water_anchors else fitresult.params[f"d18O_{sanitize(s)}"].value
			for s in samples.index
		})
		samples['SE_'+'water_'+d18O_key_out] = samples.index.map({
			s: None if s in water_anchors else fitresult.params[f"d18O_{sanitize(s)}"].stderr
			for s in samples.index
		})
		samples['95CL_'+'water_'+d18O_key_out] = samples['SE_'+'water_'+d18O_key_out] * self.standardization['water_'+D17O_key_out]['t95']
		samples['SD_'+'water_'+d18O_key_out] = (
			self.data
			.groupby('Sample')
			.agg({'CO2_'+d18O_key_out+'_corrected': 'std'})
		)

		samples.sort_index(inplace = True)

		samples['water_'+D17O_key_out] = samples.index.map({
			s: water_anchors[s][D17O_key_out] if s in water_anchors else fitresult.params[f"D17O_{sanitize(s)}"].value
			for s in samples.index
		})
		samples['SE_'+'water_'+D17O_key_out] = samples.index.map({
			s: None if s in water_anchors else fitresult.params[f"D17O_{sanitize(s)}"].stderr
			for s in samples.index
		})
		samples['95CL_'+'water_'+D17O_key_out] = samples['SE_'+'water_'+D17O_key_out] * self.standardization['water_'+D17O_key_out]['t95']
		samples['SD_'+'water_'+D17O_key_out] = (
			self.data
			.groupby('Sample')
			.agg({'CO2_'+D17O_key_out+'_corrected': 'std'})
		)

		samples['water_'+d17O_key_out] = 1e3 * (
				numpy.exp((samples['water_'+D17O_key_out]) / 1e3)
				* (1 + samples['water_'+d18O_key_out] / 1e3)**lambda_17
				- 1
			)

		# TODO: propagate uncertainties to d17O_key_out using the code below
		#
		# try:
		# 	i = fitresult.var_names.index('d18O_'+_s)
		# 	j = fitresult.var_names.index('D17O_'+_s)
		# 	CM = fitresult.covar[ix_([i,j],[i,j])]
		#
		# 	J = array([
		# 		[
		# 			exp(out['samples'][s][D17O_key_out]/1e3)
		# 			* ISOPARAMETERS.lambda_17
		# 			* (1+out['samples'][s][d18O_key_out]/1e3) ** (ISOPARAMETERS.lambda_17 - 1)
		# 			],
		# 		[
		# 			exp(out['samples'][s][D17O_key_out]/1e3)
		# 			* (1+out['samples'][s][d18O_key_out]/1e3) ** ISOPARAMETERS.lambda_17
		# 			],
		# 		])
		#
		# 	out['samples'][s]['SE_'+d17O_key_out] = float(J.T @ CM @ J)**.5
		# 	out['samples'][s]['95CL_'+d17O_key_out] = out['samples'][s]['SE_'+d17O_key_out] * out['t95']
		# except ValueError:
		# 	pass

		if hasattr(self, 'samples'):
			columns = self.samples.columns.tolist()
			for c in samples.columns:
				if c not in columns:
					columns.append(c)
			self.samples = pandas.concat([self.samples, samples], axis = 'columns')
		else:
			self.samples = samples.copy()


		df_anchors = pandas.DataFrame([{
			'water_'+d18O_key_out: water_anchors[s][d18O_key_out],
			'water_'+D17O_key_out: water_anchors[s][D17O_key_out],
			} for s in water_anchors], index = water_anchors.keys()
		)

		if hasattr(self, 'anchors'):
# 			self.anchors = self.anchors.combine_first(df_anchors)
			self.anchors = pandas.concat([self.anchors, df_anchors], axis = 'columns')
		else:
			self.anchors = df_anchors.copy()
			self.anchors.index.name = 'Sample'

		# TODO: compute anchor mismatch ie mean and SE of residuals


def D17_CO2_H2O_eq(
	h2o_to_co2_molecular_ratio,
	d18i_co2,
	D17i_co2,
	d18i_h2o,
	D17i_h2o,
	d13_co2,
	alpha18_CO2g_H2O,
	alpha17_CO2g_H2O,
	R13_VPDB,
	R18_VSMOW,
	R17_VSMOW,
	lambda_17,
	D727i = 0.1043,
	D728i = 0.1951,
	D828i = 0.3743,
	D637i = 0.5101,
	D638i = 0.9714,
	D737i = 1.1188,
	D738i = 1.6662,
	D838i = 2.3046,
	D727=0.1043,
	D728=0.1951,
	D828=0.3743,
	D637=0.5101,
	D638=0.9714,
	D737=1.1188,
	D738=1.6662,
	D838=2.3046,
	stochastic = False,
	xtol = 1e-12,
):
	if stochastic:
		rK637, rK637i = 1, 1
		rK638, rK638i = 1, 1
		rK727, rK727i = 1, 1
		rK728, rK728i = 1, 1
		rK737, rK737i = 1, 1
		rK738, rK738i = 1, 1
		rK828, rK828i = 1, 1
		rK838, rK838i = 1, 1
	else:
		rK637, rK637i = 1 + D637 / 1e3, 1 + D637i / 1e3
		rK638, rK638i = 1 + D638 / 1e3, 1 + D638i / 1e3
		rK727, rK727i = 1 + D727 / 1e3, 1 + D727i / 1e3
		rK728, rK728i = 1 + D728 / 1e3, 1 + D728i / 1e3
		rK737, rK737i = 1 + D737 / 1e3, 1 + D737i / 1e3
		rK738, rK738i = 1 + D738 / 1e3, 1 + D738i / 1e3
		rK828, rK828i = 1 + D828 / 1e3, 1 + D828i / 1e3
		rK838, rK838i = 1 + D838 / 1e3, 1 + D838i / 1e3

	x_co2 = 1 / (1 + h2o_to_co2_molecular_ratio)
	x_h2o = 1 - x_co2

	r13i_co2 = R13_VPDB * (1 + d13_co2 / 1e3)  # initial [13]/[12]
	r18i_co2 = R18_VSMOW * (1 + d18i_co2 / 1e3)  # initial [18]/[16]
	r17i_co2 = (
		numpy.exp(D17i_co2 / 1e3) * R17_VSMOW * (1 + d18i_co2 / 1e3) ** lambda_17
	)  # initial [17]/[16]

	r628i_co2 = 2 * r18i_co2  # initial [628]/[626]
	r627i_co2 = 2 * r17i_co2  # initial [627]/[626]
	r828i_co2 = rK828i * r18i_co2**2  # initial [828]/[626]
	r727i_co2 = rK727i * r17i_co2**2  # initial [727]/[626]
	r728i_co2 = rK728i * 2 * r17i_co2 * r18i_co2  # initial [728]/[626]
	r636i_co2 = r13i_co2  # initial [628]/[626]
	r638i_co2 = rK638i * 2 * r18i_co2 * r13i_co2  # initial [628]/[626]
	r637i_co2 = rK637i * 2 * r17i_co2 * r13i_co2  # initial [627]/[626]
	r838i_co2 = rK838i * r18i_co2**2 * r13i_co2  # initial [828]/[626]
	r737i_co2 = rK737i * r17i_co2**2 * r13i_co2  # initial [727]/[626]
	r738i_co2 = (
		rK738i * 2 * r17i_co2 * r18i_co2 * r13i_co2
	)  # initial [727]/[626]                                  # initial [728]/[626]

	x626i_co2 = x_co2 / (
		1
		+ r628i_co2
		+ r627i_co2
		+ r828i_co2
		+ r727i_co2
		+ r728i_co2
		+ r636i_co2
		+ r638i_co2
		+ r637i_co2
		+ r838i_co2
		+ r737i_co2
		+ r738i_co2
	)  # initial fraction of [626] over all molecules
	x628i_co2 = x626i_co2 * r628i_co2  # initial fraction of [628] over all molecules
	x627i_co2 = x626i_co2 * r627i_co2  # initial fraction of [627] over all molecules
	x828i_co2 = x626i_co2 * r828i_co2  # initial fraction of [828] over all molecules
	x727i_co2 = x626i_co2 * r727i_co2  # initial fraction of [727] over all molecules
	x728i_co2 = x626i_co2 * r728i_co2  # initial fraction of [728] over all molecules
	x636i_co2 = x626i_co2 * r636i_co2  # initial fraction of [628] over all molecules
	x638i_co2 = x626i_co2 * r638i_co2  # initial fraction of [628] over all molecules
	x637i_co2 = x626i_co2 * r637i_co2  # initial fraction of [627] over all molecules
	x838i_co2 = x626i_co2 * r838i_co2  # initial fraction of [828] over all molecules
	x737i_co2 = x626i_co2 * r737i_co2  # initial fraction of [727] over all molecules
	x738i_co2 = x626i_co2 * r738i_co2  # initial fraction of [728] over all molecules

	r18i_h2o = R18_VSMOW * (1 + d18i_h2o / 1e3)  # initial [181]/[161]
	r17i_h2o = (
		numpy.exp(D17i_h2o / 1e3) * R17_VSMOW * (1 + d18i_h2o / 1e3) ** lambda_17
	)  # initial [171]/[161]
	x16i_h2o = x_h2o / (
		1 + r18i_h2o + r17i_h2o
	)  # initial fraction of [161] over all molecules
	x17i_h2o = x16i_h2o * r17i_h2o
	x18i_h2o = x16i_h2o * r18i_h2o

	def equations(p):
		(
			x626_co2,
			x627_co2,
			x628_co2,
			x727_co2,
			x828_co2,
			x728_co2,
			x636_co2,
			x637_co2,
			x638_co2,
			x737_co2,
			x838_co2,
			x738_co2,
			x16_h2o,
			x17_h2o,
			x18_h2o,
		) = numpy.exp(p)
		return (
			# conservation of 13C
			(
				x636_co2
				+ x637_co2
				+ x638_co2
				+ x737_co2
				+ x738_co2
				+ x838_co2
				- x636i_co2
				- x637i_co2
				- x638i_co2
				- x737i_co2
				- x738i_co2
				- x838i_co2
			),
			# conservation of 16O
			(
				x16_h2o
				+ 2 * x626_co2
				+ x627_co2
				+ x628_co2
				+ 2 * x636_co2
				+ x637_co2
				+ x638_co2
				- x16i_h2o
				- 2 * x626i_co2
				- x627i_co2
				- x628i_co2
				- 2 * x636i_co2
				- x637i_co2
				- x638i_co2
			),
			# conservation of 17O
			(
				x17_h2o
				+ 2 * x727_co2
				+ x627_co2
				+ x728_co2
				+ 2 * x737_co2
				+ x637_co2
				+ x738_co2
				- x17i_h2o
				- 2 * x727i_co2
				- x627i_co2
				- x728i_co2
				- 2 * x737i_co2
				- x637i_co2
				- x738i_co2
			),
			# conservation of 18O
			(
				x18_h2o
				+ 2 * x828_co2
				+ x628_co2
				+ x728_co2
				+ 2 * x838_co2
				+ x638_co2
				+ x738_co2
				- x18i_h2o
				- 2 * x828i_co2
				- x628i_co2
				- x728i_co2
				- 2 * x838i_co2
				- x638i_co2
				- x738i_co2
			),
			# alpha18 = ([628]/[626]/2) / ([181]/[616])
			alpha18_CO2g_H2O * x18_h2o / x16_h2o - x628_co2 / 2 / x626_co2,
			# alpha17 = ([627]/[626]/2) / ([171]/[616])
			alpha17_CO2g_H2O * x17_h2o / x16_h2o - x627_co2 / 2 / x626_co2,
			# conservation of H2O
			x16_h2o + x17_h2o + x18_h2o - x_h2o,
			# clumped isotopes
			x828_co2 / x626_co2 - rK828 * (x628_co2 / x626_co2 / 2) ** 2,
			x727_co2 / x626_co2 - rK727 * (x627_co2 / x626_co2 / 2) ** 2,
			x728_co2 / x626_co2 - rK728 * x627_co2 / x626_co2 * x628_co2 / x626_co2 / 2,
			x638_co2 / x626_co2 - rK638 * (x628_co2 / x626_co2) * x636_co2 / x626_co2,
			x637_co2 / x626_co2 - rK637 * (x627_co2 / x626_co2) * x636_co2 / x626_co2,
			x838_co2 / x626_co2
			- rK838 * (x628_co2 / x626_co2 / 2) ** 2 * x636_co2 / x626_co2,
			x737_co2 / x626_co2
			- rK737 * (x627_co2 / x626_co2 / 2) ** 2 * x636_co2 / x626_co2,
			x738_co2 / x626_co2
			- rK738
			* x636_co2
			/ x626_co2
			* x627_co2
			/ x626_co2
			* x628_co2
			/ x626_co2
			/ 2,
		)

	(
		x626_co2,
		x627_co2,
		x628_co2,
		x727_co2,
		x828_co2,
		x728_co2,
		x636_co2,
		x637_co2,
		x638_co2,
		x737_co2,
		x838_co2,
		x738_co2,
		x16_h2o,
		x17_h2o,
		x18_h2o,
	) = numpy.exp(
		fsolve(
			equations,
			numpy.log(
				numpy.array(
					[
						x626i_co2,
						x627i_co2,
						x628i_co2,
						x727i_co2,
						x828i_co2,
						x728i_co2,
						x636i_co2,
						x637i_co2,
						x638i_co2,
						x737i_co2,
						x838i_co2,
						x738i_co2,
						x16i_h2o,
						x17i_h2o,
						x18i_h2o,
					]
				)
			),
			xtol = xtol,
		)
	)

	d18e_h2o = 1e3 * ((x18_h2o / x16_h2o / R18_VSMOW) - 1)
	D17e_h2o = 1e3 * (
		numpy.log(x17_h2o / x16_h2o / R17_VSMOW)
		- lambda_17 * numpy.log(x18_h2o / x16_h2o / R18_VSMOW)
	)
	d18e_co2 = 1e3 * ((x628_co2 / x626_co2 / 2 / R18_VSMOW) - 1)
	D17e_co2 = 1e3 * (
		numpy.log(x627_co2 / x626_co2 / 2 / R17_VSMOW)
		- lambda_17 * numpy.log(x628_co2 / x626_co2 / 2 / R18_VSMOW)
	)

	return d18e_co2, D17e_co2, d18e_h2o, D17e_h2o

if __name__ == '__main__':

	from random import gauss

	N = 100
	_df = pandas.DataFrame([
		{
			'UID': f'X{_+1:03.0f}',
			'Session': f'2024-{_//20+1:02.0f}',
			'Sample': ['IAEA603', 'IAEA612', 'NBS18', 'FOO', 'BAR'][_ % 5],
			'd636': [2.00, -35.00, -6., 35., -20.][_ % 5] + 0.1 * gauss(),
			'd628': [-2.00, -11.00, -23.0, 0., 20.][_ % 5] + 0.1 * gauss(),
			'D627': 0.01 * gauss(),
		}
		for _ in range(N)
	])

	_df['d627']	= (numpy.exp(_df['D627']/1e3) * (1+_df['d628']/1e3)**0.528 - 1) * 1e3

	df = Dataset(_df)

	df.standardize('d636', 'd13C_VPDB', dict(IAEA603 = 2.46, IAEA612 = -36.722))
	df.standardize('d628', 'd18O_VPDB', dict(IAEA603 = -2.37, NBS18 = -23.01))

	df.standardize_D17O(
		'd628', 'd627', 'd18O_VSMOW', 'd17O_VSMOW',
		dict(
			IAEA603 = dict(d18O_VSMOW = -2, D17O_VSMOW = 0.),
			NBS18 = dict(d18O_VSMOW = -23., D17O_VSMOW = 0.),
		),
		gamma_17 = -0.1273,
	)

	df.plot_D17O_residuals()

# 	df.standardize_D17O(
# 		'd628', 'd627', 'd18O_IAEA603', 'd17O_IAEA603',
# 		dict(
# 			IAEA603 = dict(d18O_IAEA603 = 0., D17O_IAEA603 = 0.000),
# 			NBS18 = dict(d18O_IAEA603 = -21., D17O_IAEA603 = 0.027),
# 		),
# 	)
#
# 	df.plot_D17O_residuals(d18O_key_out = 'd18O_IAEA603', D17O_key_out = 'D17O_IAEA603')

	pyplot.show()

	df.to_csv('foo.csv', float_format = '%.6f')

	print(df.data)
	for k in df.sessions:
		print()
		print(f'[{k}]')
		print(df.sessions[k])
	for k in df.standardization:
		print()
		print(f'[{k}]')
		print({_: df.standardization[k][_] for _ in df.standardization[k] if _ != 'report'})
	print()
	print(df.samples)
	print()
	print(df.anchors)
