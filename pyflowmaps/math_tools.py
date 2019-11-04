
import numpy as np
from scipy.fftpack import fft2, ifft2,ifftshift
from astropy.convolution import convolve, Box2DKernel

__all__ = ['fivepoint','qfit2','crossD','divergence','vorticity','fft_differentiation', 'fft_poisson', 'smooth']
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


	if (dim < 2) and (dim > 4):
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
		a1 = .5*a1-cc[0,1,:,:]-cc[1,1:,:]-cc[2,1,:,:]
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


	if (dim < 2) and (dim > 4):
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
	du_y,du_x = np.gradient(vx)
	dv_y,dv_x = np.gradient(vy)
	div = du_x+dv_y

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
	du_y,du_x = np.gradient(vx)
	dv_y,dv_x = np.gradient(vy)
	vort = dv_x-du_y

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

	min = int((w-1)/2)
	max_x = int(nx - (w+1)/2)
	max_y = int(ny - (w+1)/2)

	kernel = Box2DKernel(w)
	smooth_0 = convolve(im2d,kernel)
	smooth_f = im2d.copy()
	smooth_f[min:max_y+1,min:max_x+1] = smooth_0[min:max_y+1,min:max_x+1]

	return smooth_f
