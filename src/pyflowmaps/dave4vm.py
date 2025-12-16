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

import numpy as np
from .math_tools import the_matrix


def pydave4vm(magnetic_vector,window_size,dx,dy):
	"""Core DAVE4VM solver producing local coefficients over pixels.

	Parameters
	----------
	magnetic_vector : dict
		Keys: 'bx','by','bz','bxx','bxy','byx','byy','bzx','bzy','bzt'. Arrays shape (ny, nx).
	window_size : int
		Odd window size in pixels for local averaging/inversion.
	dx, dy : float
		Pixel scales in km/pixel.

	Returns
	-------
	dict
		Velocity components and local coefficients over the input grid.
	"""
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
		AA = AM[:, :, j, i]
		GA = AA[0:9, 0:9]
		FA = -1 * np.reshape(AA[9, 0:9], 9)
		vector = np.linalg.solve(GA, FA)
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



