'''Parametrizations of hadronic cross sections as described in the thesis document.'''

import numpy as np

def total_hadron_proton_scattering(s, h):
	'''Return the total hadron proton scattering cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV
	h : {'p', 'pi', 'K'}
		The incident hadron on the proton target

	Returns
	-------
	float
		The total hadron proton scattering cross section for `h` at `s` in mb
	'''
	return 0.1

def hadron_elastic_total_ratio(s):
	'''Return the universal ratio of elastic to total hadron proton cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV

	Returns
	-------
	float
		The universal ratio of elastic to total hadron proton cross section
	'''
	return 0.2

def inelastic_hadron_proton_scattering(s, h):
	'''Return the inelastic hadron proton scattering cross section.

	Parameters
	----------
	s : float
		The squared center of mass energy in GeV
	h : {'p', 'pi', 'K'}
		The incident hadron on the proton target

	Returns
	-------
	float
		The inelastic hadron proton scattering cross section for `h` at `s` in mb
	'''
	return total_hadron_proton_scattering(s, h) * (1 - hadron_elastic_total_ratio(s))

