#!/usr/bin/python
# -*- coding: utf8 -*-


#Class flowmaker
import numpy as np
from .math_tools import fivepoint, qfit2, crossD, fft_convolve2d_precomputed
from scipy.ndimage import correlate, correlate1d
from astropy.convolution import Box2DKernel, Gaussian1DKernel, Gaussian2DKernel
from collections import namedtuple
try:
	from sunpy.image.resample import resample as _sunpy_resample
	_HAS_SUNPY = True
except Exception:
	_HAS_SUNPY = False

try:
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	from sunpy.coordinates import Helioprojective
	_HAS_SUNPY_COORDS = True
except Exception:
	_HAS_SUNPY_COORDS = False

def _resample(image, new_shape, method='nearest', center=None, minusone=None):
	"""Resize 2D array to new_shape using SunPy when available, else SciPy.

	Parameters
	----------
	image : 2D array
		Input image.
	new_shape : tuple(int, int)
		Target shape (ny, nx).
	method : {'nearest', 'linear'}
		Interpolation method. 'nearest' maps to order=0 in SciPy fallback.
	center, minusone : any
		Accepted for compatibility; ignored in SciPy fallback.
	"""
	if _HAS_SUNPY:
		return _sunpy_resample(image, new_shape, method=method, minusone=(False if minusone is None else minusone))
	else:
		from scipy.ndimage import zoom
		zoom_factors = (new_shape[0] / image.shape[0], new_shape[1] / image.shape[1])
		order = 0 if method == 'nearest' else 1
		return zoom(image, zoom_factors, order=order)


def _spatially_varying_convolve(data, std_base, mu_map, theta_map):
	"""Gaussian convolution with per-pixel varying elliptical kernel.

	The image is convolved multiple times, once for each bin in a
	discretised (mu, theta) grid.  Only pixels whose local geometry
	falls in a given bin take their values from the corresponding
	convolution, producing a smooth per-pixel correction.

	Performance: the image FFT is computed **once** and reused for every
	(mu, theta) bin.  Each bin only requires a small kernel FFT + one
	element-wise multiply + one inverse FFT, making this ~10-50x faster
	than spatial-domain convolution for 700x700 images.

	Parameters
	----------
	data : 2D array
		Input image to convolve.
	std_base : float
		Base Gaussian standard deviation (pixels) at mu = 1.
	mu_map : 2D array
		Per-pixel cos(heliocentric angle).
	theta_map : 2D array
		Per-pixel position angle (radians) of radial direction.

	Returns
	-------
	result : 2D array
		Convolved image.
	"""
	from scipy.fft import fft2, next_fast_len

	ny, nx = data.shape

	# Fast path: foreshortening negligible everywhere
	if mu_map.min() > 0.95:
		kernel = Gaussian1DKernel(stddev=std_base).array
		return np.rot90(correlate1d(np.rot90(
			correlate1d(data, kernel), 3), kernel), 1)

	result = np.zeros_like(data)

	mu_lo, mu_hi = float(mu_map.min()), float(mu_map.max())
	theta_lo, theta_hi = float(theta_map.min()), float(theta_map.max())

	# Adaptive bin counts: ~0.05 in mu, ~5 deg in theta (coarser for speed)
	n_mu = min(max(1, int(np.ceil((mu_hi - mu_lo) / 0.05))), 20)
	n_th = min(max(1, int(np.ceil((theta_hi - theta_lo) / np.radians(5)))), 20)

	mu_edges = np.linspace(mu_lo, mu_hi + 1e-10, n_mu + 1)
	th_edges = np.linspace(theta_lo, theta_hi + 1e-10, n_th + 1)

	mu_idx = np.clip(np.digitize(mu_map, mu_edges) - 1, 0, n_mu - 1)
	th_idx = np.clip(np.digitize(theta_map, th_edges) - 1, 0, n_th - 1)

	# Build a reference kernel to know the max kernel size
	ref_kernel = Gaussian2DKernel(
		x_stddev=std_base, y_stddev=max(std_base * mu_lo, 0.5)
	).array
	ky, kx = ref_kernel.shape

	# Pre-compute padded FFT of the image (done ONCE)
	py = next_fast_len(ny + ky - 1)
	px = next_fast_len(nx + kx - 1)
	pad_shape = (py, px)
	data_fft = fft2(data, s=pad_shape)

	for mi in range(n_mu):
		for ti in range(n_th):
			mask = (mu_idx == mi) & (th_idx == ti)
			if not np.any(mask):
				continue
			mu_avg = 0.5 * (mu_edges[mi] + mu_edges[mi + 1])
			th_avg = 0.5 * (th_edges[ti] + th_edges[ti + 1])
			kernel_2d = Gaussian2DKernel(
				x_stddev=std_base, y_stddev=std_base * mu_avg, theta=th_avg
			).array
			convolved = fft_convolve2d_precomputed(
				data_fft, kernel_2d, pad_shape, (ny, nx))
			result[mask] = convolved[mask]

	return result


def _deproject_velocities(vx, vy, mu_map, theta_map):
	"""Deproject apparent image-plane velocities to solar-surface velocities.

	The LCT measures displacements in the image plane.  The radial
	component (towards/away from disk centre) is foreshortened by mu.
	This routine corrects that per pixel.

	Parameters
	----------
	vx, vy : 2D arrays
		Apparent velocity components in image-plane coordinates.
	mu_map : 2D array
		Per-pixel cos(heliocentric angle).
	theta_map : 2D array
		Per-pixel position angle (rad, CCW from +x) of radial direction.

	Returns
	-------
	vx_corr, vy_corr : 2D arrays
		Deprojected velocity components.
	"""
	cos_t = np.cos(theta_map)
	sin_t = np.sin(theta_map)

	# Decompose into radial / tangential
	v_rad = vx * cos_t + vy * sin_t
	v_tan = -vx * sin_t + vy * cos_t

	# Correct radial component for foreshortening
	v_rad = v_rad / mu_map

	# Reconstruct Cartesian components
	vx_corr = v_rad * cos_t - v_tan * sin_t
	vy_corr = v_rad * sin_t + v_tan * cos_t

	return vx_corr, vy_corr


def compute_mu_theta_maps(shape, scale, crpix, obstime=None, observer='earth',
						  rsun=None):
	"""Compute per-pixel mu and theta maps from image geometry.

	Works in two modes:

	* **Pure numpy** (no SunPy needed): provide ``rsun`` as a float
	  (apparent solar radius in arcsec).  The maps are computed with
	  simple Euclidean geometry on the projected disk.
	* **SunPy mode**: provide ``obstime`` (and optionally ``observer``).
	  SunPy's `Helioprojective` frame computes the accurate apparent
	  solar radius from the observer ephemeris.

	Parameters
	----------
	shape : tuple (ny, nx)
		Image shape in pixels.
	scale : float
		Plate scale in arcsec/pixel.
	crpix : tuple (crpix_x, crpix_y)
		Pixel location of Sun centre.  FITS convention (1-based).
	obstime : str, `~astropy.time.Time`, or None
		Observation time.  When given **and** SunPy is installed, the
		apparent solar radius and observer geometry are computed from
		ephemeris.  When ``None``, the pure-numpy path is used and
		``rsun`` must be provided.
	observer : str or coordinate, optional
		Observer location (only used in SunPy mode).  Default ``'earth'``.
	rsun : float, `~astropy.units.Quantity`, or None
		Apparent solar angular radius in arcsec.  Required when
		``obstime`` is ``None``.  In SunPy mode, if ``None`` the value
		is obtained from the ephemeris; if given it overrides it.

	Returns
	-------
	mu_map : 2D ndarray of shape ``(ny, nx)``
		Per-pixel cos(heliocentric angle).  Off-disk pixels are NaN.
	theta_map : 2D ndarray of shape ``(ny, nx)``
		Per-pixel position angle (radians, CCW from +x) toward disk
		centre.

	Examples
	--------
	Pure numpy (no SunPy):

	>>> mu, theta = compute_mu_theta_maps((512, 512), 0.6, (256.5, 256.5),
	...                                   rsun=960.0)

	With SunPy ephemeris:

	>>> mu, theta = compute_mu_theta_maps((512, 512), 0.6, (256.5, 256.5),
	...                                   obstime='2023-06-15T12:00:00')
	"""
	ny, nx = shape
	scale_val = float(scale.value if hasattr(scale, 'value') else scale)

	# Build helioprojective Tx, Ty in arcsec
	ix = np.arange(nx, dtype=np.float64)
	iy = np.arange(ny, dtype=np.float64)
	ix_grid, iy_grid = np.meshgrid(ix, iy)

	tx = (ix_grid - (crpix[0] - 1)) * scale_val   # arcsec
	ty = (iy_grid - (crpix[1] - 1)) * scale_val   # arcsec

	# Determine apparent solar radius
	if obstime is not None and _HAS_SUNPY_COORDS:
		# --- SunPy path: accurate ephemeris ---
		frame_kw = dict(obstime=obstime, observer=observer)
		if rsun is not None:
			frame_kw['rsun'] = rsun if hasattr(rsun, 'unit') else rsun * u.arcsec
		hp = Helioprojective(**frame_kw)
		# angular_radius needs a SkyCoord with data
		sc = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=hp)
		rsun_arcsec = float(sc.angular_radius.to(u.arcsec).value)
	else:
		# --- Pure numpy path ---
		if rsun is None:
			raise ValueError(
				'rsun (apparent solar radius in arcsec) is required when '
				'obstime is not provided.  Typical value ~960 arcsec.')
		rsun_arcsec = float(rsun.value if hasattr(rsun, 'value') else rsun)

	rho = np.sqrt(tx**2 + ty**2)
	rho_frac = rho / rsun_arcsec

	on_disk = rho_frac < 1.0
	mu_map = np.where(on_disk, np.sqrt(1.0 - rho_frac**2), np.nan)
	theta_map = np.arctan2(-ty, -tx)

	return mu_map, theta_map


def compute_mu_theta_from_map(smap):
	"""Compute per-pixel mu and theta maps from a SunPy Map.

	This is a convenience wrapper around `compute_mu_theta_maps` that
	extracts all necessary parameters (plate scale, reference pixel,
	observation time, observer) directly from a ``sunpy.map.GenericMap``
	object, ensuring consistency with the FITS header metadata.

	Parameters
	----------
	smap : `~sunpy.map.GenericMap`
		A SunPy Map (e.g. loaded from a FITS file via ``sunpy.map.Map()``).

	Returns
	-------
	mu_map : 2D ndarray
		Per-pixel cos(heliocentric angle).  Off-disk pixels are NaN.
	theta_map : 2D ndarray
		Per-pixel position angle (radians, CCW from +x).

	Raises
	------
	ImportError
		If SunPy is not available.

	Example
	-------
	>>> import sunpy.map
	>>> smap = sunpy.map.Map('hmi_continuum.fits')
	>>> mu, theta = compute_mu_theta_from_map(smap)
	>>> vx, vy = pyflowmaker(cube, fwhm, mu=mu, theta=theta)
	"""
	if not _HAS_SUNPY_COORDS:
		raise ImportError(
			'compute_mu_theta_from_map requires sunpy.  '
			'Install with: pip install sunpy')

	ny, nx = smap.data.shape

	# Build pixel grids (0-based)
	ix = np.arange(nx, dtype=np.float64)
	iy = np.arange(ny, dtype=np.float64)
	ix_grid, iy_grid = np.meshgrid(ix, iy)

	# Use the Map's WCS for pixel → world conversion
	coords = smap.pixel_to_world(ix_grid.ravel() * u.pix,
								 iy_grid.ravel() * u.pix)

	# Helioprojective Tx, Ty in arcsec
	tx = coords.Tx.to(u.arcsec).value.reshape(ny, nx)
	ty = coords.Ty.to(u.arcsec).value.reshape(ny, nx)

	# Angular distance from disk centre and solar angular radius
	rho = np.sqrt(tx**2 + ty**2)
	rsun_ang = coords[0].angular_radius.to(u.arcsec).value

	rho_frac = rho / rsun_ang
	on_disk = rho_frac < 1.0

	mu_map = np.where(on_disk, np.sqrt(1.0 - rho_frac**2), np.nan)
	theta_map = np.arctan2(-ty, -tx)

	return mu_map, theta_map


__all__ = ['pyflowmaker', 'compute_mu_theta_maps', 'compute_mu_theta_from_map']
__authors__ = ["Jose Ivan Campos Rozo, Santiago Vargas Dominguez"]
__email__ = "hypnus1803@gmail.com"

method = ['square','absolute','corr']
interpolation = ['fivepoint','qfit2','crossD']
window = ['gaussian','boxcar']

VelocityPair = namedtuple('VelocityPair', 'vx vy')


def pyflowmaker(mc,fwhm,reb=1, lag=1, method='square', interpolation = 'fivepoint', window = 'gaussian',
				mu=1.0, theta=0.0):
	"""
	Compute flow maps and returns X and Y components for the proper
	motion map.

	Parameters
	----------

	mc :
        A `~numpy array` of shape ``(nt, ny,nx)`` where ``nt`` is the number of
        layers or images to calculate the flow maps.
	fwhm :
	    Apodization window size in pixels.
	lag :
	    time-lag between 'references' and 'life' subseries. Default 1.
	reb :
	    rebinning factor to change scale/range of November's method. Default 1
	method : {'square' | 'absolute' | 'corr'}
	      Methods to calculate local diferences.
		   * square - Sum of squares of the local differences.
		   * absolute - Apply an absolute differences algorithm.
		   * cross -  uses a multiplicative algorithm over the differences.
		  Default is ``square`` method.
	interpolation : {'fivepoint' | 'qfit2' | 'crossD'}
	      Mode to calculate the extrema values in 3x3 matrix.
		   * fivepoint - Measure the position of minimum or maximum in a 3x3 matrix.
		   * qfit2 - Uses 9 points fitting to measure the extrema value in a 3x3 matrix.
		   * crossD - Measure the extrema using cross derivative interpolation in 3x3 matrix.
		  Default is ``fivepoint`` mode.
	window : {'gaussian' | 'boxcar'}
	      Kind of apodization window during the calculation fo the flow maps.
		   * gaussian - 2-dimensional correlation with a gaussian filter kernel.
		     When ``mu < 1.0``, an elliptical Gaussian kernel is used to account
		     for foreshortening near the solar limb (see ``mu`` and ``theta``).
		   * boxcar - 2-dimensional correlation with a box filter kernel.
	mu : float or 2D array, optional
	    Cosine of the heliocentric angle (0 < mu <= 1).  Can be a scalar
	    (single value for the entire FOV) or a 2D array of shape ``(ny, nx)``
	    giving per-pixel values.  When a 2D array is provided, a spatially
	    varying elliptical Gaussian kernel is used and the resulting
	    velocities are deprojected per pixel.  Default 1.0.
	theta : float or 2D array, optional
	    Position angle of the radial direction in radians, measured
	    counterclockwise from the positive x-axis.  Must match the type of
	    ``mu`` (scalar or 2D array of the same shape).  Default 0.0.

	Returns
	-------
	vx, vy : 2D velocity components (vx, vy) for the proper displacements map.



	Example:
	--------
			>>> vx,vy=flowmaker(cube,1,8*u.pix,1*u.pix)
			>>> # Near-limb observation with scalar mu=0.7, radial direction at 45 degrees
			>>> vx,vy=flowmaker(cube,1,8*u.pix,1*u.pix, mu=0.7, theta=np.pi/4)
			>>> # Per-pixel mu/theta maps (2D arrays of shape (ny, nx))
			>>> vx,vy=flowmaker(cube,1,8*u.pix,1*u.pix, mu=mu_map, theta=theta_map)

	"""

	# Normalize and validate parameters

	if method not in ('square', 'absolute', 'corr'):
		raise ValueError('Acceptable method keywords are "absolute" | "corr" | "square"; method "square" is the default.')
	if interpolation not in ('fivepoint', 'qfit2', 'crossD'):
		raise ValueError('Acceptable mode keywords are "fivepoint" | "qfit2" | "crossD"; mode "fivepoint" is the default.')
	if window not in ('gaussian', 'boxcar'):
		raise ValueError('Acceptable window keywords are "boxcar" | "gaussian"; window "gaussian" is the default.')

	shf=1
	std1=fwhm/(2*np.sqrt(2*np.log(2)))
	std2=std1/reb

	dims = np.ndim(mc)
	n_im = len(mc)
	yy = mc.shape[1]
	yy_r = int(yy/reb)
	xx =  mc.shape[2]
	xx_r = int(xx/reb)

	# --- Per-pixel or scalar geometry ---
	mu_is_map = isinstance(mu, np.ndarray) and mu.ndim == 2
	theta_is_map = isinstance(theta, np.ndarray) and theta.ndim == 2

	if mu_is_map != theta_is_map:
		raise ValueError('mu and theta must both be scalars or both be 2D arrays.')

	if mu_is_map:
		if mu.shape != (yy, xx):
			raise ValueError(
				'mu map shape {} does not match image shape ({}, {}).'.format(mu.shape, yy, xx))
		if theta.shape != (yy, xx):
			raise ValueError(
				'theta map shape {} does not match image shape ({}, {}).'.format(theta.shape, yy, xx))
		if np.any(mu <= 0) or np.any(mu > 1.0):
			raise ValueError('All mu values must be in the range (0, 1].')
		if reb != 1:
			from scipy.ndimage import zoom as _zoom
			_zf = (yy_r / yy, xx_r / xx)
			mu_r = _zoom(mu, _zf, order=1)
			theta_r = _zoom(theta, _zf, order=1)
		else:
			mu_r = mu
			theta_r = theta
	else:
		if not (0 < mu <= 1.0):
			raise ValueError('mu must be in the range (0, 1], got {}'.format(mu))


	n = int(n_im-lag)

	n_p = xx_r*yy_r

	cc=np.zeros((3,3,yy_r,xx_r))


	for k in range(n):
		if reb != 1:
			map_a = _resample(mc[k,:,:], (yy_r, xx_r), method='nearest', minusone=False)
			map_b = _resample(mc[k+lag], (yy_r, xx_r), method='nearest', minusone=False)
		else:
			map_a = mc[k,:,:]
			map_b = mc[k+lag,:,:]
		map_a = map_a - np.nansum(map_a)/n_p
		map_b = map_b - np.nansum(map_b )/n_p

		for i in range(-1,2):
			for j in range(-1,2):

				if method == 'absolute':
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+abs(np.roll(map_a,(i*shf,j*shf),axis=(1,0))-np.roll(map_b,(-i*shf,-j*shf),axis=(1,0)))

				elif method == 'corr':
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:] + np.roll(map_a,(i*shf,j*shf),axis=(1,0)) - np.roll(map_b,(-i*shf,-j*shf),axis=(1,0))

				elif method == 'square':
					dumb = np.roll(map_a,(i*shf,j*shf),axis=(1,0)) - np.roll(map_b,(-i*shf,-j*shf),axis=(1,0))
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+dumb*dumb
					del dumb
				else:
					raise ValueError('Aceptable method keywords are "absolute" | "corr" | "square"; method "square" is the default.')
		del map_a, map_b

	#X boundaries
	cc[:,:,:,0]=cc[:,:,:,1]
	cc[:,:,:,xx_r-1]=cc[:,:,:,xx_r-2]
	
	#Y boundaries
	cc[:,:,0,:]=cc[:,:,1,:]
	cc[:,:,yy_r-1,:]=cc[:,:,yy_r-2,:]

	for i in range(3):
		for j in range(3):
			if window == 'boxcar':
				boxcar = Box2DKernel(fwhm/reb).array
				cc[j,i,:,:] = correlate(cc[j,i,:,:],boxcar)
			elif window == 'gaussian':
				if mu_is_map:
					cc[j,i,:,:] = _spatially_varying_convolve(cc[j,i,:,:], std2, mu_r, theta_r)
				elif mu < 1.0:
					kernel_2d = Gaussian2DKernel(x_stddev=std2, y_stddev=std2*mu, theta=theta).array
					cc[j,i,:,:] = correlate(cc[j,i,:,:], kernel_2d)
				else:
					kernel = Gaussian1DKernel(stddev=std2).array
					cc[j,i,:,:] = np.rot90(correlate1d(np.rot90(correlate1d(cc[j,i,:,:],kernel),3),kernel),1)
			else:
				raise ValueError('Acceptable window keywords are "boxcar" | "gaussian"; window "gaussian" is the default.')

	if interpolation =='qfit2':
		vx,vy=qfit2(cc)
	elif interpolation =='crossD':
		vx,vy=crossD(cc)
	elif interpolation == 'fivepoint':
		vx,vy=fivepoint(cc)
	else:
		raise ValueError('Acceptable mode keywords are "fivepoint" | "qfit2" | "crossD"; mode "fivepoint" is the default.')

	vx = 2.*shf*vx
	vy = 2.*shf*vy

	if reb != 1:
		vx = _resample(vx, (yy, xx), method='nearest')*reb
		vy = _resample(vy, (yy, xx), method='nearest')*reb

	vx = vx*reb
	vy = vy*reb

	# Deproject velocities when per-pixel geometry maps are provided
	if mu_is_map:
		vx, vy = _deproject_velocities(vx, vy, mu, theta)

	return VelocityPair(vx,vy)
from sunpy.coordinates import frames