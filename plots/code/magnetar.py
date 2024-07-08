'''
Object oriented implementation of magnetar class as described in the thesis document.

	Classes
	-------
	magnetar

'''
import numpy as np
import scipy.constants as con

import calculations as cal


class magnetar(R, B, O, chi, I):
	'''
	Collection of parameters and associated methods with the magnetar model.

	Attributes
	----------
	R : float
		The stellar radius in cm
	B : float
		The polar magnetic field strength in G
	O : float
		The angular frequency in rad / s
	chi : float
		The relative dipole tilt to rotational axis angle in rad
	I : float
		The moment of inertia in g * cm**2
	'''
