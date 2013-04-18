
import numpy as np
# Converts polar to cartesian coordinates using numpy types rather than cmath
def rect(r, phi):

	x = r * np.cos(phi)
	y = r * np.sin(phi)

	return x + (1j * y)


def cphase(a):

	return np.arctan2(a.imag, a.real )