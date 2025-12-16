#!/usr/bin/python
# -*- coding: utf8 -*-
import numpy as np
from scipy.fft import fft2, ifft2,ifftshift
from astropy.convolution import convolve, Box2DKernel
from scipy.ndimage import gaussian_filter

__all__ = ['fivepoint','qfit2','crossD','divergence','vorticity','fft_differentiation', 'fft_poisson', 'smooth',
		   'odiff','the_matrix','pointingflux','neutralline','v_perp']
__authors__ = ["Jose Ivan Campos Rozo"]
__email__ = ["hypnus1803@gmail.com"]



def fivepoint(cc):
	"""
		Measure the position of minimum or maximum in a 3x3 matrix.
		Written by Roberto Luis Molowny Horas, Institute of Theoretical
		Astrophysics, University of Oslo. August 1991.
		This version is based on the IDL version.
		Adapted in Python by J. I. Campos-Rozo

		Parameters
		----------
				 CC : array-like
					  Cross correlation function. It must have dimensions like
					  CC[3,3,:,:], or CC[3,3,:], or CC[3,3]
		Returns
		-------
				 x,y : Position of the minimum using a 2nd order polynomial fit.
					   The center is taken on CC[1,1,:,:].
	"""

	dim = np.ndim(cc)

	if dim < 2 or dim>4:
		raise ValueError('Wrong input array dimensions')
	if cc.shape[0] != 3 or cc.shape[1] != 3:
		raise ValueError('Array-shape must be CC[3,3,:,:], or CC[3,3,:], or CC[3,3]')

	if dim == 4:
		y = 2*cc[1,1,:,:]
		x = (cc[1,0,:,:]-cc[1,2,:,:])/(cc[1,2,:,:]+cc[1,0,:,:]-y)*0.5
		y=(cc[0,1,:,:]-cc[2,1,:,:])/(cc[2,1,:,:]+cc[0,1,:,:]-y)*0.5


	elif dim == 3:
		y=2.*cc[1,1,:]
		x=(cc[1,0,:]-cc[1,2,:])/(cc[1,2,:]+cc[1,0,:]-y)*0.5
		y=(cc[0,1,:]-cc[2,1,:])/(cc[2,1,:]+cc[0,1,:]-y)*0.5

	elif dim == 2:
		y=2.*cc[1,1]
		x=(cc[1,0]-cc[1,2])/(cc[1,2]+cc[1,0]-y)*0.5
		y=(cc[0,1]-cc[2,1])/(cc[2,1]+cc[0,1]-y)*0.5

	return x,y



def qfit2(cc):
	"""
		Measure the position of extrem value in a 3x3 matrix
		Fits 9 points in a 2D array following the interpolation
		polynomial function:
				f(x,y) = c1*f +c2*x*f + c3*y*f + c4*x**2*f + c5*y**2*f + c6
		Following FORTRAN program qfit2.fts on PE by L.J. November,
		NSO/SP, 1986.

		References:
		T. A. Darvann, 1991, Master Thesis, University of Oslo

		Taken from the IDL version written by R. Molowny Horas.
		Adapated for Python by J. I. Campos-Rozo
	"""


	dim = np.ndim(cc)


	if (dim < 2) or (dim > 4):
		raise ValueError('Wrong input array dimasions')
	if cc.shape[0] != 3 or cc.shape[1] != 3:
		raise ValueError('Array-shape must be CC[3,3,:,:], or CC[3,3,:], or CC[3,3]')

	if dim == 4 :
		a1=cc[0,0,:,:]+cc[0,2,:,:]+cc[2,0,:,:]+cc[2,2,:,:]
		a2 = a1+cc[0,1,:,:]+cc[2,1,:,:]
		a1 = a1+cc[1,0,:,:]+cc[1,2,:,:]
		a3 = cc[0,0,:,:]-cc[0,2,:,:]-cc[2,0,:,:]+cc[2,2,:,:]
		a4 = -cc[0,0,:,:]+cc[2,2,:,:]
		a5 = a4-cc[0,1,:,:]-cc[0,2,:,:]+cc[2,0,:,:]+cc[2,1,:,:]
		a4 = a4+cc[0,2,:,:]-cc[1,0,:,:]+cc[1,2,:,:]-cc[2,0,:,:]
		a1 = .5*a1-cc[0,1,:,:]-cc[1,1,:,:]-cc[2,1,:,:]
		a2 = .5*a2-cc[1,0,:,:]-cc[1,1,:,:]-cc[1,2,:,:]
	elif dim == 3:
		a1 = cc[0,0,:]+cc[0,2,:]+cc[2,0,:]+cc[2,2,:]
		a2 = a1+cc[0,1,:]+cc[2,1,:]
		a1 = a1+cc[1,0,:]+cc[1,2,:]
		a3 = cc[0,0,:]-cc[0,2,:]-cc[2,0,:]+cc[2,2,:]
		a4 = -cc[0,0,:]+cc[2,2,:]
		a5 = a4-cc[0,1,:]-cc[0,2,:]+cc[2,0,:]+cc[2,1,:]
		a4 = a4+cc[0,2,:]-cc[1,0,:]+cc[1,2,:]-cc[2,0,:]
		a1 = .5*a1-cc[0,1,:]-cc[1,1,:]-cc[2,1,:]
		a2 = .5*a2-cc[1,0,:]-cc[1,1,:]-cc[1,2,:]
	elif dim == 2:
		a1=cc[0,0]+cc[0,2]+cc[2,0]+cc[2,2]
		a2=a1+cc[0,1]+cc[2,1]
		a1=a1+cc[1,0]+cc[1,2]
		a3=cc[0,0]-cc[0,2]-cc[2,0]+cc[2,2]
		a4=-cc[0,0]+cc[2,2]
		a5=a4-cc[0,1]-cc[0,2]+cc[2,0]+cc[2,1]
		a4=a4+cc[0,2]-cc[1,0]+cc[1,2]-cc[2,0]
		a1=.5*a1-cc[0,1]-cc[1,1]-cc[2,1]
		a2=.5*a2-cc[1,0]-cc[1,1]-cc[1,2]

	dim_ = ((64./9)*a1*a2-a3**2)*1.5

	cx = (a3*a5-((8./3)*a2*a4))/dim_
	cy = (a3*a4-8./3*a1*a5)/dim_
	return cx, cy

def crossD(cc):
	"""
		Measure the position of a minimum or maximum in a 3x3 array. This function
		is an extension of the fivepoint function, but, taking into account the
		four points in the corners. The interpolation is made by a polynomial function like:
				f(x,y) = c1 +c2*x + c3*y + c4*x**2 + c5*y**2 + c6*x*y
		Written by Roberto Luis Molowny Horas, August 1991.
		Adapted for Python by J. I. Campos-Rozo
		References:
				Yi, Z. and Molowny H., R. 1992

		Parameters
		----------
				 CC : array-like
					  Cross correlation function. It must have dimensions like
					  CC[3,3,:,:], or CC[3,3,:], or CC[3,3]
		Returns
		-------
				 x,y : Position of the extrema using a 2nd order polynomial fit.
					   The center is taken on CC[1,1,:,:].

	"""
	dim = np.ndim(cc)


	if (dim < 2) or (dim > 4):
		raise ValueError('Wrong input array dimasions')
	if cc.shape[0] != 3 or cc.shape[1] != 3:
		raise ValueError('Array-shape must be CC[3,3,:,:], or CC[3,3,:], or CC[3,3]')

	if dim == 2:
		c4=cc[1,2]+cc[1,0]-cc[1,1]*2.
		c2=cc[1,2]-cc[1,0]
		c5=cc[2,1]+cc[0,1]-cc[1,1]*2.
		c3=cc[2,1]-cc[0,1]
		c6=(cc[2,2]-cc[2,0]-cc[0,2]+cc[0,0])/4.
	elif dim==3:
		c4=cc[1,2,:]+cc[1,0,:]-cc[1,1,:]*2.
		c2=cc[1,2,:]-cc[1,0,:]
		c5=cc[2,1,:]+cc[0,1,:]-cc[1,1,:]*2
		c3=cc[2,1,:]-cc[0,1,:]
		c6=(cc[2,2,:]-cc[2,0,:]-cc[0,2,:]+cc[0,0,:])/4.
	elif dim ==4:
		c4=cc[1,2,:,:]+cc[1,0,:,:]-cc[1,1,:,:]*2.
		c2=cc[1,2,:,:]-cc[1,0,:,:]
		c5=cc[2,1,:,:]+cc[0,1,:,:]-cc[1,1,:,:]*2
		c3=cc[2,1,:,:]-cc[0,1,:,:]
		c6=(cc[2,2,:,:]-cc[2,0,:,:]-cc[0,2,:,:]+cc[0,0,:,:])/4.
	determ=0.5/(c4*c5 - c6*c6)
	x=determ*(c6*c3 - c5*c2)
	y=determ*(c6*c2 - c4*c3)
	return x, y



def divergence(vx,vy):

	"""
		Make divergence between two 2-D velocity maps v_z = h_m*div(vx,vy),
		where $h_m$ is defined as mass-flux scale-heigth (November 1989, ApJ,344,494),
		and div(vx,vy) = $\frac{\partial v_x}{\partial x} + \frac{\partial v_y}{\partial y}$

		Parameters
		----------
				vx : map of velocities in x direction.
				vy : map of velocities in y direction.
		Results
		-------
				vz : map of velocities in z direction.

	"""
	du_x, du_y, dv_x, dv_y = _velocity_gradients(vx, vy)
	div = du_x + dv_y

	return div
def vorticity(vx,vy):

	"""
		Calculate the vorticity between two 2-D velocity maps v_z = h_m*vort(vx,vy),
		where $h_m$ is defined as mass-flux scale-heigth (November 1989, ApJ,344,494),
		and vort(vx,vy) = $\frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}$

		Parameters
		----------
				vx : map of velocities in x direction.
				vy : map of velocities in y direction.
		Results
		-------
				vz : map of velocities in z direction.

	"""
	du_x, du_y, dv_x, dv_y = _velocity_gradients(vx, vy)
	vort = dv_x - du_y

	return vort

def fft_differentiation(image, dx=1.0, dy = None):

	"""
		Return the differentiation in x, and y axis of an 2-dimensional array.
		The differentiation is computed using forward differentiation, that in
		the case of the frequency domain, should apply the shift theorem from
		the properties of the Discrete Fourier Transform.

		Parameters
		----------
				image : array_like
						A 2-dimensional array containing an image or samples of
						a scalar function.
				dx : constant sample distance for x-dimension. Default is 1.0
				dy : constant sample distance for y-dimension. If dy is not given,
					 dy = dx. Default is None.
		Results
		-------
				fft_differentiation : list with the differentiations of image with
				respect to each axis.

		Note
		----
				A full explanation about how the fft_differentiation function works
				is given in Stack Overflow.
				https://stackoverflow.com/questions/29189885/finding-the-derivative-of-a-2d-function-using-fft-properties
	"""


	if dy is None:
		dy = dx

	ny, nx = image.shape

	kx = np.linspace(-nx/2,nx/2,nx+1)/nx
	ky = np.linspace(-ny/2,ny/2,ny+1)/ny

	kx = kx[:-1]
	ky = ky[:-1]

	if np.any(np.isnan(image)):
		image = np.nan_to_num(image)

	data_wavenumberdomain = fft2(image)


	# Compute grid of wavenumbers
	KX, KY = np.meshgrid(ifftshift(kx),ifftshift(ky))

	# Compute 2D derivative in x-direction
	dx_wavenumberdomain =  1j*2*np.pi*KX*data_wavenumberdomain
	# Convert back to space domain
	dx_spacedomain = ifft2(dx_wavenumberdomain)

	dx_spacedomain = np.real(dx_spacedomain)/dx

	# Compute 2D derivative in y-direction
	dy_wavenumberdomain =  1j*2*np.pi*KY*data_wavenumberdomain
	# Convert back to space domain
	dy_spacedomain = ifft2(dy_wavenumberdomain)

	dy_spacedomain = np.real(dy_spacedomain)/dy

	return [dx_spacedomain, dy_spacedomain]




def fft_poisson(data,dx=1.0,dy=None):

	"""
		A 2-dimensional poisson equation solver using
		Discrete Fourier Transformation (DFT).

		Parameters
		----------
				data : A 2-dimensional array containing the data we want to solve.
				dx : constant sample distance for x-dimension. Default is 1.0
				dy : constant sample distance for y-dimension. If dy is not given,
					 dy = dx. Default is None.
		Reuslts
		-------
				fft_poisson : Solution of the poisson equation using forward differentiation,
				in the frequency domain (DFT).
	"""

	if dy is None:
		dy = dx

	ny, nx = data.shape

	kx = np.linspace(-nx/2,nx/2,nx+1)/nx
	ky = np.linspace(-ny/2,ny/2,ny+1)/ny

	kx = kx[:-1]
	ky = ky[:-1]

	if np.any(np.isnan(data)):
		data = np.nan_to_num(data)

	data_wavenumberdomain = np.ma.array(fft2(data))
	data_wavenumberdomain[0,0] = np.ma.masked

	# Compute grid of wavenumbers except in the point (0,0)
	KX, KY = np.ma.array(np.meshgrid(ifftshift(kx),ifftshift(ky)))
	KX[0,0] = np.ma.masked
	KY[0,0] = np.ma.masked


	# Compute 2D derivative
	data_wavenumberdomain_differentiated =  -data_wavenumberdomain/((4*np.pi**2)*(KX**2 + KY**2))

	# Convert back to space domain
	data_spacedomain_differentiated = ifft2(data_wavenumberdomain_differentiated )

	return (np.real(data_spacedomain_differentiated))*(dx*dy)

def smooth(im2d, w):

	ny, nx = im2d.shape

	if w%2 == 0:
		w = w+1

	mini = int((w-1)/2)
	max_x = int(nx - (w+1)/2)
	max_y = int(ny - (w+1)/2)

	kernel = Box2DKernel(w)
	smooth_0 = convolve(im2d,kernel,normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	smooth_f = im2d.copy()
	smooth_f[mini:max_y+1,mini:max_x+1] = smooth_0[mini:max_y+1,mini:max_x+1]

	return smooth_f

# -------------------------------------------------------------
# DAVE4VM helper routines (ported)
# -------------------------------------------------------------

def odiff(arr):
	"""
	Optimized finite-difference in x and y using weighted rolls.

	Returns
	-------
	dx, dy : 2D arrays
	"""
	c1, c2 = 0.12019, 0.74038

	# Rolling the array in all directions
	a1, a2 = np.roll(arr, -2, axis=1), np.roll(arr, -1, axis=1)
	a3, a4 = np.roll(arr, 1, axis=1), np.roll(arr, 2, axis=1)
	b1, b2 = np.roll(arr, -2, axis=0), np.roll(arr, -1, axis=0)
	b3, b4 = np.roll(arr, 1, axis=0), np.roll(arr, 2, axis=0)

	# Calculating the differentials dx and dy
	dx = -c1*a1 + c2*a2 - c2*a3 + c1*a4
	dy = -c1*b1 + c2*b2 - c2*b3 + c1*b4

	return dx, dy


def the_matrix(bx, bxx, bxy, by, byx, byy, bz, bzx, bzy,
			   bzt, psf, psfx, psfy, psfxx, psfyy, psfxy):
	"""
	Perform convolutions and construct the LKA matrix from which solutions are calculated.
	"""
	# Precompute products
	bz2 = bz * bz
	GGx = bz * bzx
	GGy = bz * bzy
	GGxx = bzx * bzx
	GGt = bzt * bz

	GGyy = bzy * bzy
	GGxy = bzx * bzy
	GGtx = bzt * bzx
	GGty = bzt * bzy
	GGtt = bzt * bzt

	# Convolutions
	G   = convolve(bz2, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	Gx  = convolve(GGx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGx = convolve(GGx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGx = convolve(GGx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	Gy  = convolve(GGy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGy = convolve(GGy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGy = convolve(GGy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	Ht  = convolve(GGt, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	Gxx = convolve(GGxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	Gyy = convolve(GGyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	Gxy = convolve(GGxy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	Gtx = convolve(GGtx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	Gty = convolve(GGty, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGxx = convolve(GGxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGyy = convolve(GGyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGxy = convolve(GGxy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	xGtx = convolve(GGtx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xGty = convolve(GGty, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGxx = convolve(GGxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGyy = convolve(GGyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	yGxy = convolve(GGxy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGtx = convolve(GGtx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yGty = convolve(GGty, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xxGxx = convolve(GGxx, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	xxGxy = convolve(GGxy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xxGyy = convolve(GGyy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xyGxx = convolve(GGxx, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	xyGxy = convolve(GGxy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	xyGyy = convolve(GGyy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yyGxx = convolve(GGxx, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yyGxy = convolve(GGxy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	yyGyy = convolve(GGyy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	Gtt = convolve(GGtt, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	Bxx = bx*bx
	Byy = by*by
	Bxy = bx*by
	Bzx = bz*bx
	Bzy = bz*by
	mbxbxx = bx*bxx
	mbxbyy = bx*byy
	mbxxbxx = bxx*bxx
	mbyybyy = byy*byy
	mbxxbyy = bxx*byy
	mbybxx = by*bxx
	mbybyy = by*byy
	mbzbxx = bz*bxx
	mbzbyy = bz*byy
	mbztbxx = bzt*bxx
	mbztbyy = bzt*byy
	mbzxbx = bzx*bx
	mbzxby = bzx*by
	mbzxbxx = bzx*bxx
	mbzxbyy = bzx*byy
	mbzybx = bzy*bx
	mbzyby = bzy*by
	mbzybxx = bzy*bxx
	mbzybyy = bzy*byy
	mbztbx = bzt*bx
	mbztby = bzt*by

	BxBx = convolve(Bxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	BxBy = convolve(Bxy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	ByBy = convolve(Byy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	BzBx = convolve(Bzx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	
	BzBy = convolve(Bzy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	BxBxx = convolve(mbxbxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	BxByy = convolve(mbxbyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')
	BxxBxx = convolve(mbxxbxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')

	ByyByy = convolve(mbyybyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.9
	BxxByy = convolve(mbxxbyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.10
	ByBxx = convolve(mbybxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.11
	ByByy = convolve(mbybyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.12

	BzBxx = convolve(mbzbxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.13
	BzByy = convolve(mbzbyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.14
	BztBxx = convolve(mbztbxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.15
	BztByy = convolve(mbztbyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.16

	BzxBx = convolve(mbzxbx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.17
	BzxBy = convolve(mbzxby, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.18
	BzxBxx = convolve(mbzxbxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.19
	BzxByy = convolve(mbzxbyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.20

	BzyBx = convolve(mbzybx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.21
	BzyBy = convolve(mbzyby, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.22
	BzyBxx = convolve(mbzybxx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.23
	BzyByy = convolve(mbzybyy, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.24

	BztBx = convolve(mbztbx, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.25
	BztBy = convolve(mbztby, psf, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.26
	xBzxBx = convolve(mbzxbx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.27
	xBzxBy = convolve(mbzxby, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.28

	xBzyBx = convolve(mbzybx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.29
	xBzyBy = convolve(mbzyby, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.30
	yBzyBx = convolve(mbzybx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.31
	yBzyBy = convolve(mbzyby, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.32

	yBzxBx = convolve(mbzxbx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.33
	yBzxBy = convolve(mbzxby, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.34
	yBxBxx = convolve(mbxbxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.35
	yBxByy = convolve(mbxbyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.36

	yByBxx = convolve(mbybxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.37
	yByByy = convolve(mbybyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.38
	xByBxx = convolve(mbybxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.39
	xByByy = convolve(mbybyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.40

	xBzxBxx = convolve(mbzxbxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.41
	xBzxByy = convolve(mbzxbyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.42
	yBzxBxx = convolve(mbzxbxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.43
	yBzxByy = convolve(mbzxbyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.44

	xBxxBxx = convolve(mbxxbxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.45
	xBxxByy = convolve(mbxxbyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.46
	xByyByy = convolve(mbyybyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.47
	yBxxBxx = convolve(mbxxbxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.48

	yBxxByy = convolve(mbxxbyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.49
	yByyByy = convolve(mbyybyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.50
	xBxBxx = convolve(mbxbxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.51
	xBxByy = convolve(mbxbyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.52

	xBzBxx = convolve(mbzbxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.53
	xBzByy = convolve(mbzbyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.54
	xBztBxx = convolve(mbztbxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.55
	xBztByy = convolve(mbztbyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.56

	yBztBxx = convolve(mbztbxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.57
	yBztByy = convolve(mbztbyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.58
	xyBxxBxx = convolve(mbxxbxx, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.59
	xyBxxByy = convolve(mbxxbyy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.60

	xyByyByy = convolve(mbyybyy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.61
	xyBzxBxx = convolve(mbzxbxx, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.62
	xyBzxByy = convolve(mbzxbyy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.63
	xyBzyBxx = convolve(mbzybxx, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.64

	xyBzyByy = convolve(mbzybyy, psfxy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.65
	yBzBxx = convolve(mbzbxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.66
	yBzByy = convolve(mbzbyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.67
	xBzyBxx = convolve(mbzybxx, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.68

	xBzyByy = convolve(mbzybyy, psfx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.69
	yBzyBxx = convolve(mbzybxx, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.70
	yBzyByy = convolve(mbzybyy, psfy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.71
	xxBxxBxx = convolve(mbxxbxx, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.72

	xxBxxByy = convolve(mbxxbyy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.73
	xxByyByy = convolve(mbyybyy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.74
	xxBzxBxx = convolve(mbzxbxx, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.75
	xxBzyBxx = convolve(mbzybxx, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.76

	xxBzxByy = convolve(mbzxbyy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.77
	xxBzyByy = convolve(mbzybyy, psfxx, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.78
	yyBxxBxx = convolve(mbxxbxx, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.79
	yyBxxByy = convolve(mbxxbyy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.80

	yyByyByy = convolve(mbyybyy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.81
	yyBzyBxx = convolve(mbzybxx, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.82
	yyBzyByy = convolve(mbzybyy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.83
	yyBzxBxx = convolve(mbzxbxx, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.84

	yyBzxByy = convolve(mbzxbyy, psfyy, normalize_kernel=False, preserve_nan=True, nan_treatment='fill')  # b.85


	A = np.stack([
		Gxx, Gxy, Gx + xGxx, Gx + yGxy, yGxx, xGxy, -BzxBxx - BzxByy, -BzxBx - xBzxBxx - xBzxByy,
		-BzxBy - yBzxBxx - yBzxByy, Gtx, Gxy, Gyy, Gy + xGxy, Gy + yGyy, yGxy, xGyy, -BzyBxx - BzyByy,
		-BzyBx - xBzyBxx - xBzyByy, -BzyBy - yBzyBxx - yBzyByy, Gty, Gx + xGxx, Gy + xGxy,
		G + 2 * xGx + xxGxx, G + xGx + xyGxy + yGy, xyGxx + yGx, xGy + xxGxy,
		-BzBxx - BzByy - xBzxBxx - xBzxByy, -BzBx - xBzBxx - xBzByy - xBzxBx - xxBzxBxx - xxBzxByy,
		-BzBy - xBzxBy - xyBzxBxx - xyBzxByy - yBzBxx - yBzByy, Ht + xGtx, Gx + yGxy, Gy + yGyy,
		G + xGx + xyGxy + yGy, G + 2 * yGy + yyGyy, yGx + yyGxy, xGy + xyGyy,
		-BzBxx - BzByy - yBzyBxx - yBzyByy, -BzBx - xBzBxx - xBzByy - xyBzyBxx - xyBzyByy - yBzyBx,
		-BzBy - yBzBxx - yBzByy - yBzyBy - yyBzyBxx - yyBzyByy, Ht + yGty, yGxx, yGxy, xyGxx + yGx,
		yGx + yyGxy, yyGxx, xyGxy, -yBzxBxx - yBzxByy, -xyBzxBxx - xyBzxByy - yBzxBx,
		-yBzxBy - yyBzxBxx - yyBzxByy, yGtx, xGxy, xGyy, xGy + xxGxy, xGy + xyGyy, xyGxy, xxGyy,
		-xBzyBxx - xBzyByy, -xBzyBx - xxBzyBxx - xxBzyByy, -xBzyBy - xyBzyBxx - xyBzyByy, xGty,
		-BzxBxx - BzxByy, -BzyBxx - BzyByy, -BzBxx - BzByy - xBzxBxx - xBzxByy,
		-BzBxx - BzByy - yBzyBxx - yBzyByy, -yBzxBxx - yBzxByy, -xBzyBxx - xBzyByy,
		BxxBxx + 2 * BxxByy + ByyByy, BxBxx + BxByy + xBxxBxx + 2 * xBxxByy + xByyByy,
		ByBxx + ByByy + yBxxBxx + 2 * yBxxByy + yByyByy, -BztBxx - BztByy, -BzxBx - xBzxBxx - xBzxByy,
		-BzyBx - xBzyBxx - xBzyByy, -BzBx - xBzBxx - xBzByy - xBzxBx - xxBzxBxx - xxBzxByy,
		-BzBx - xBzBxx - xBzByy - xyBzyBxx - xyBzyByy - yBzyBx,
		-xyBzxBxx - xyBzxByy - yBzxBx, -xBzyBx - xxBzyBxx - xxBzyByy,
		BxBxx + BxByy + xBxxBxx + 2 * xBxxByy + xByyByy,
		BxBx + 2 * xBxBxx + 2 * xBxByy + xxBxxBxx + 2 * xxBxxByy + xxByyByy,
		BxBy + xByBxx + xByByy + xyBxxBxx + 2 * xyBxxByy + xyByyByy + yBxBxx + yBxByy,
		-BztBx - xBztBxx - xBztByy, -BzxBy - yBzxBxx - yBzxByy, -BzyBy - yBzyBxx - yBzyByy,
		-BzBy - xBzxBy - xyBzxBxx - xyBzxByy - yBzBxx - yBzByy,
		-BzBy - yBzBxx - yBzByy - yBzyBy - yyBzyBxx - yyBzyByy,
		-yBzxBy - yyBzxBxx - yyBzxByy, -xBzyBy - xyBzyBxx - xyBzyByy,
		ByBxx + ByByy + yBxxBxx + 2 * yBxxByy + yByyByy,
		BxBy + xByBxx + xByByy + xyBxxBxx + 2 * xyBxxByy + xyByyByy + yBxBxx + yBxByy,
		ByBy + 2 * yByBxx + 2 * yByByy + yyBxxBxx + 2 * yyBxxByy + yyByyByy, -BztBy - yBztBxx - yBztByy,
		Gtx, Gty, Ht + xGtx, Ht + yGty, yGtx, xGty, -BztBxx - BztByy, -BztBx - xBztBxx - xBztByy,
		-BztBy - yBztBxx - yBztByy, Gtt
	], axis=0)

	return A


def pointingflux(dx,bx,by,bz,vx,vy,vz):
	"""Energy flux of magnetic and velocity fields (SI units)."""
	gauss2tesla = 1e-4
	km2m = 1e3
	mu0 = 4e-7*np.pi
	area_px = (dx*km2m)**2

	Bx = bx*gauss2tesla
	By = by*gauss2tesla
	Bz = bz*gauss2tesla

	Vx = vx*km2m
	Vy = vy*km2m
	Vz = vz*km2m

	Sn = area_px*((Bx**2 + By**2)*Vz)/mu0
	St = area_px*(-Bz*(Vx*Bx + Vy*By))/mu0
	Ss = Sn + St

	return Sn, St, Ss, np.sum(Sn), np.sum(St), np.sum(Ss)


def neutralline(bz,gaussian=False,fwhm_arc=3,scale=0.058,mag_qs=150):
	"""Compute neutral line mask for Bz field."""
	bz_positive = convolve((bz>mag_qs).astype('float32'),np.ones((9,9)),mode='same',method='fft')
	bz_negative = convolve((bz<-mag_qs).astype('float32'),np.ones((9,9)),mode='same',method='fft')

	nl = ((bz_positive>0.5) + (bz_negative>0.5))>1.5

	if gaussian:
		fwhm = fwhm_arc/scale
		std1 = fwhm / (2 * np.sqrt(2 * np.log(2)))
		nl = gaussian_filter(nl, std1, order=0)

	return nl


def v_perp(bx,by,bz,vx,vy,vz):
	"""Compute velocity component perpendicular to magnetic field."""
	factor = -1*(bx*vx + by*vy + bz*vz)/(bx**2 + by**2 + bz**2)
	vpx = vx + factor*bx
	vpy = vy + factor*by
	vpz = vz + factor*bz
	return vpx, vpy, vpz


# Helper: compute velocity gradients once when both divergence and vorticity are needed
def _velocity_gradients(vx, vy):
	"""
	Compute gradients for velocity components.

	Parameters
	----------
	vx, vy : 2D arrays
		Velocity components on a regular grid.

	Returns
	-------
	du_x, du_y, dv_x, dv_y : 2D arrays
		Partial derivatives of vx and vy.
	"""
	du_y, du_x = np.gradient(vx)
	dv_y, dv_x = np.gradient(vy)
	return du_x, du_y, dv_x, dv_y
