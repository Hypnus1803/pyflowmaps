#!/usr/bin/python
# -*- coding: utf8 -*-


#Class flowmaker
import numpy as np
from .math_tools import fivepoint, qfit2, crossD
from scipy.ndimage import correlate, correlate1d
from astropy.convolution import Box2DKernel, Gaussian1DKernel, Gaussian2DKernel
from collections import namedtuple
try:
	from sunpy.image.resample import resample as _sunpy_resample
	_HAS_SUNPY = True
except Exception:
	_HAS_SUNPY = False

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

__all__ = ['pyflowmaker']
__authors__ = ["Jose Ivan Campos Rozo, Santiago Vargas Dominguez"]
__email__ = "hypnus1803@gmail.com"

method = ['square','absolute','corr']
interpolation = ['fivepoint','qfit2','crossD']
window = ['gaussian','boxcar']

VelocityPair = namedtuple('VelocityPair', 'vx vy')


def pyflowmaker(mc,fwhm,reb=1, lag=1, method='square', interpolation = 'fivepoint', window = 'gaussian'):
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
		   * boxcar - 2-dimensional correlation with a box filter kernel.

	Returns
	-------
	vx, vy : 2D velocity components (vx, vy) for the proper displacements map.



	Example:
	--------
			>>> vx,vy=flowmaker(cube,1,8*u.pix,1*u.pix)

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

	return VelocityPair(vx,vy)
