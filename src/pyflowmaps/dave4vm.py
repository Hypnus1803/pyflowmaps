"""
DAVE4VM core solver integrated into pyflowmaps.

This module provides the core pixel-wise DAVE4VM computation. It expects precomputed
magnetic field components and their spatial derivatives, plus a window size and pixel scales.

Units
-----
- Magnetic components (`bx`, `by`, `bz`) in Gauss.
- Spatial derivatives (`bxx`, `bxy`, `byx`, `byy`, `bzx`, `bzy`) in 1/km.
- Time derivative `bzt` in Gauss/s; include `dt` in seconds if needed upstream.
- `dx`, `dy` in km/pixel.
- Output velocities (`U0`, `V0`, `W0`) in km/s; local coefficients (`UX`, `UY`, `VX`, `VY`, `WX`, `WY`).
"""

import warnings
import numpy as np
from .math_tools import the_matrix
from scipy.linalg import lstsq


_REQUIRED_KEYS = {'bx', 'by', 'bz', 'bxx', 'bxy', 'byx', 'byy', 'bzx', 'bzy', 'bzt'}


def pydave4vm(magnetic_vector,window_size):
	"""Core DAVE4VM solver producing local coefficients over pixels.

	Parameters
	----------
	magnetic_vector : dict
		Keys: 'bx','by','bz','bxx','bxy','byx','byy','bzx','bzy','bzt'. Arrays shape (ny, nx).
	window_size : int
		Odd window size in pixels for local averaging/inversion.

	Returns
	-------
	dict
		Velocity components and local coefficients over the input grid.
	"""
	# --- Input validation ---
	missing = _REQUIRED_KEYS - set(magnetic_vector.keys())
	if missing:
		raise ValueError('magnetic_vector is missing required keys: {}'.format(sorted(missing)))

	ref_shape = magnetic_vector['bz'].shape
	for key in _REQUIRED_KEYS:
		if magnetic_vector[key].shape != ref_shape:
			raise ValueError(
				'Shape mismatch: magnetic_vector[\'{}\'] has shape {}, expected {}'.format(
					key, magnetic_vector[key].shape, ref_shape))

	dx = magnetic_vector['dx']
	dy = magnetic_vector['dy']

	if dx <= 0 or dy <= 0:
		raise ValueError('dx and dy must be positive, got dx={}, dy={}'.format(dx, dy))

	if not isinstance(window_size, (int, np.integer)) or window_size < 1:
		raise ValueError('window_size must be a positive integer, got {}'.format(window_size))
	if window_size % 2 == 0:
		warnings.warn('window_size should be odd for symmetric averaging; got {}. Adding 1.'.format(window_size))
		window_size += 1

	U0 = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	V0 = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	UX = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	VY = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	UY = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	VX = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	W0 = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	WX = np.zeros_like(magnetic_vector['bz'], dtype='float32')
	WY = np.zeros_like(magnetic_vector['bz'], dtype='float32')

	x, y = np.meshgrid(np.arange(-(window_size//2),window_size//2+1)*dx,
					   np.arange(-(window_size//2),window_size//2+1)*dy)
	psf = np.ones_like(x)/x.size
	psfx = psf * x
	psfy = psf * y
	psfxx = psf * x * x
	psfyy = psf * y * y
	psfxy = psf * x * y

	AM = the_matrix(magnetic_vector['bx'],  magnetic_vector['bxx'],
					magnetic_vector['bxy'], magnetic_vector['by'],
					magnetic_vector['byx'], magnetic_vector['byy'],
					magnetic_vector['bz'],  magnetic_vector['bzx'],
					magnetic_vector['bzy'], magnetic_vector['bzt'],
					psf, psfx, psfy, psfxx, psfyy, psfxy)

	AM = np.reshape(AM, (10, 10, *magnetic_vector['bz'].shape))

	trc = np.trace(AM)

	yi, xi = np.where(trc > 1.0)

	for j, i in zip(yi, xi):
		#AA = AM[:, :, j, i]
		GA = AM[0:9, 0:9,j,i]
		FA = -1 * AM[9, 0:9,j,i]
		try:
			vector = np.linalg.solve(GA, FA)
		except np.linalg.LinAlgError:
			vector = lstsq(GA, FA, lapack_driver='gelsy')[0]

			
		U0[j, i] = vector[0]
		V0[j, i] = vector[1]
		UX[j, i] = vector[2]
		VY[j, i] = vector[3]
		UY[j, i] = vector[4]
		VX[j, i] = vector[5]
		W0[j, i] = vector[6]
		WX[j, i] = vector[7]
		WY[j, i] = vector[8]

	vel4vm = {'U0': U0, 'UX': UX, 'UY': UY,
			  'V0': V0, 'VX': VX, 'VY': VY,
			  'W0': W0, 'WX': WX, 'WY': WY}

	return vel4vm



