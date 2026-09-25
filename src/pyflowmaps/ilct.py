import warnings
import numpy as np
from .math_tools import fft_differentiation, fft_poisson, smooth
from collections import namedtuple,OrderedDict





def pyilct(velocityField,BField_comp,pix_size,interval,psi_opt=False,phi_opt=False,
		   threshold=10.0,check=False, mask=False,verbose=False):

	velocityField = np.asarray(velocityField, dtype=float)
	BField_comp = np.asarray(BField_comp, dtype=float)

	if velocityField.ndim != 3 or velocityField.shape[0] != 2:
		raise ValueError('velocityField must have shape (2, ny, nx), got {}'.format(velocityField.shape))
	if BField_comp.ndim != 3 or BField_comp.shape[0] != 4:
		raise ValueError('BField_comp must have shape (4, ny, nx), got {}'.format(BField_comp.shape))
	if velocityField.shape[1:] != BField_comp.shape[1:]:
		raise ValueError('Spatial dimensions of velocityField {} and BField_comp {} must match.'.format(
			velocityField.shape[1:], BField_comp.shape[1:]))

	if pix_size <= 0:
		raise ValueError('pix_size must be positive, got {}'.format(pix_size))
	if interval == 0:
		raise ValueError('interval must be non-zero, got {}'.format(interval))

	if threshold < 0:
		raise ValueError('threshold must be non-negative, got {}'.format(threshold))

	vel = velocityField.copy()
	mag = BField_comp.copy()

	structure = OrderedDict()

	u = vel[0,:,:]
	v = vel[1,:,:]

	Bx = mag[0,:,:]
	By = mag[1,:,:]
	Bz = mag[2,:,:]

	dBz_dt = mag[3,:,:]/interval


	# Find phi
	#=========

	phi = fft_poisson(dBz_dt,pix_size)

	dphi_dx, dphi_dy = fft_differentiation(phi,pix_size)


	uxBz = u*Bz
	uyBz = v*Bz

	duyBz_dx, _ = fft_differentiation(uyBz,pix_size)
	_ , duxBz_dy = fft_differentiation(uxBz,pix_size)

	curlulBz = duyBz_dx - duxBz_dy

	psi = fft_poisson(-curlulBz,pix_size)

	dpsi_dx, dpsi_dy = fft_differentiation(psi, pix_size)


	if threshold <= 1.0:
		nzthr = threshold*np.nanmax(np.abs(Bz))
	else:
		nzthr = threshold

	hiBz = np.where(np.abs(Bz) > nzthr)
	n_hiBz = len(hiBz[0])

	loBz = np.where(np.abs(Bz) <= nzthr)

	zeros = np.where(np.abs(Bz) < 1e-4)
	n_zeros = len(zeros[0])


	if n_zeros != 0:
		Bz[zeros] = np.nan 
		if verbose:
			print('Some zeros present, ILCT are ajusting them...')
	else:
		if verbose:
			print('No zeros present in the calculations')


	Fx = -dphi_dx + dpsi_dy
	Fy = -dphi_dy - dpsi_dx

	magb2 = Bx**2 + By**2 + Bz**2

	# Guard against division by zero in magb2 and Bz
	magb2_safe = magb2.copy()
	magb2_safe[magb2_safe == 0] = np.nan

	if n_hiBz == 0:
		raise ValueError('ILCT message: There is not magnetic data above the threshold limit = {0}. Stopping the code...'.format(threshold))
	else:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', RuntimeWarning)
			vz = -(Bx*Fx + By*Fy)/magb2_safe
			vx = (Fx + Bx*vz)/Bz
			vy = (Fy + By*vz)/Bz



	if check == True:
		# ~ induct_check = np.zeros((4,ny,nx),dtype=np.float64)

		induction = OrderedDict()

		x_sol = vx*Bz - vz*Bx
		y_sol = vy*Bz - vz*By

		dx_sol_dx , _ = fft_differentiation(x_sol,pix_size)
		_ , dy_sol_dy = fft_differentiation(y_sol,pix_size)

		div_vec = dx_sol_dx + dy_sol_dy

		dFx_dx , _ = fft_differentiation(Fx,pix_size)
		_ , dFy_dy = fft_differentiation(Fy,pix_size)

		div_scalars = dFx_dx + dFy_dy

		smp = 5

		duBz_dx , _ = fft_differentiation(smooth(Bz*u,smp),pix_size)
		_ , dvBz_dy = fft_differentiation(smooth(Bz*v,smp),pix_size)

		div_vel_LCT_Bz = duBz_dx + dvBz_dy

		if mask == True:
			dBz_dt[loBz] = 0.0
			div_vel_LCT_Bz[loBz] = 0.0
			div_vec[loBz] = 0.0
			div_scalars[loBz] = 0.0

		induction['dBzdt'] = dBz_dt
		induction['divLCT'] = -div_vel_LCT_Bz
		induction['divVec'] = -div_vec
		induction['divSca'] = -div_scalars
		InductionStructure = namedtuple('InductionStructure',sorted(induction))
		InductionCheck = InductionStructure(**induction)

		structure['InductionCheck'] = InductionCheck


	if mask == True:
		vx[loBz] = 0.0
		vy[loBz] = 0.0
		vz[loBz] = 0.0

		
	structure['vx'] = vx
	structure['vy'] = vy
	structure['vz'] = vz

	if psi_opt == True:
		structure['psi'] = psi
	if phi_opt == True:
		structure['phi'] = phi

	ILCTStructure = namedtuple('ILCTStructure',sorted(structure))

	return ILCTStructure(**structure)
