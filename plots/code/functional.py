'''
Vectorization of given parametrizations.

	Functions
	---------
	hadron_proton_cooling_factor
		Returns the cooling factor for hadrons scattered by protons

	proton_proton_optical_depth
		Returns the effective optical depth for protons hitting protons

	total_hadron_proton_scattering
		Returns the total hadron-proton scattering cross section

	hadron_elastic_total_ratio
		Returns the universal ratio of elastic to total hadron-proton cross section

	inelastic_hadron_proton_scattering
		Returns the inelastic hadron-proton scattering cross section

	charm_quark_differential_production
		Returns the charm quark differential cross section for production in proton-proton collisions

	charmed_hadron_differential_production
		Returns the charmed hadron differential cross section for production in proton-proton collisions

	meson_production
		Returns the proton-proton to pion or kaon singular production spectrum

	meson_decay_neutrinos
		Returns the pion or kaon to neutrino singular decay spectrum

	charmed_hadron_production
		Returns the proton-proton to charmed hadron singular production spectrum

	charmed_hadron_decay_neutrinos
		Returns the charmed hadron to neutrino singular decay spectrum

	charmed_hadron_fragmentation_function
		Returns the charmed hadrons from charm quarks fragmentation function

'''

import numpy as np

import parametrizations.collisions_decay as cd
import parametrizations.cross_sections as cr
import parametrizations.distribution_spectra as ds
import parametrizations.fragmentation_function as ff


hadron_proton_cooling_factor = np.vectorize(cd.hadron_proton_cooling_factor)
proton_proton_optical_depth = np.vectorize(cd.proton_proton_optical_depth)

total_hadron_proton_scattering = np.vectorize(cr.total_hadron_proton_scattering)
hadron_elastic_total_ratio = np.vectorize(cr.hadron_elastic_total_ratio)
inelastic_hadron_proton_scattering = np.vectorize(cr.inelastic_hadron_proton_scattering)
charm_quark_differential_production = np.vectorize(cr.charm_quark_differential_production)
charmed_hadron_differential_production = np.vectorize(cr.charmed_hadron_differential_production)

meson_production = np.vectorize(ds.meson_production)
meson_decay_neutrinos = np.vectorize(ds.meson_decay_neutrinos)
charmed_hadron_production = np.vectorize(ds.charmed_hadron_production)
charmed_hadron_decay_neutrinos = np.vectorize(ds.charmed_hadron_decay_neutrinos)

charmed_hadron_fragmentation_function = np.vectorize(ff.charmed_hadron_fragmentation_function)
