import numpy as np
from math_tools import fft_differentiation, fft_poisson, smooth
from collections import namedtuple,OrderedDict
import sys




def pyilct(velocityField,BField_comp,pix_size,interval,psi_opt=False,phi_opt=False,
		   threshold=100.0,check=False, mask=False):

	vel = velocityField.copy()
	mag = BField_comp.copy()

	structure = OrderedDict()
	ny, nx, n = mag.shape

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

	# ~ d2phi_dx2,dphi_dxdy = fft_differentiation(dphi_dx,pix_size.value)
	# ~ d2phi_dydx, d2phi_dy2 = fft_differentiation(dphi_dy,pix_size)
	# ~ laplace_phi = d2phi_dx2 + d2phi_dy2

	uxBz = u*Bz
	uyBz = v*Bz

	duyBz_dx, _ = fft_differentiation(uyBz,pix_size)
	_ , duxBz_dy = fft_differentiation(uxBz,pix_size)

	curlulBz = duyBz_dx - duxBz_dy

	psi = fft_poisson(-curlulBz,pix_size)

	dpsi_dx, dpsi_dy = fft_differentiation(psi, pix_size)

	# ~ d2psi_dx2,dpsi_dxdy = fft_differentiation(dpsi_dx,pix_size.value)
	# ~ d2psi_dydx, d2psi_dy2 = fft_differentiation(dpsi_dy,pix_size)
	# ~ laplace_psi = d2psi_dx2 + d2psi_dy2

	
	if threshold <= 1.0:
		nzthr = threshold*np.max(np.abs(Bz))
	else:
		nzthr = threshold

	hiBz = np.where(np.abs(Bz) > nzthr)
	n_hiBz = len(hiBz[0])

	loBz = np.where(np.abs(Bz) <= nzthr)
	n_loBz = len(loBz[0])

	zeros = np.where(np.abs(Bz) < 1e-4)
	n_zeros = len(zeros[0])


	if n_zeros != 0:
		Bzmiss = 1e10*np.max(np.abs(Bz))
		Bz[zeros] = Bzmiss
		print('Some zeros present, ILCT are ajusting them...')
	else:
		print('No zeros present in the calculations')

	# ~ vx = np.zeros([ny,nx],dtype=np.float64)*u.cm/u.s
	# ~ vy = np.zeros([ny,nx],dtype=np.float64)*u.cm/u.s
	# ~ vz = np.zeros([ny,nx],dtype=np.float64)*u.cm/u.s

	Fx = -dphi_dx + dpsi_dy
	Fy = -dphi_dy - dpsi_dx

	magb2 = Bx**2 + By**2 + Bz**2


	if n_hiBz == 0:
		raise ValueError('ILCT message: There is not magnetic data above the threshold limit = {0}. Stopping the code...'.format(threshold))
	else:
		vz = -(Bx*Fx + By*Fy)/magb2
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

	# ~ vel = np.zeros([3,ny,nx],dtype=np.float64)
	# ~ vel[0,:,:] = vx[:,:]
	# ~ vel[1,:,:] = vy[:,:]
	# ~ vel[2,:,:] = vz[:,:]

	structure['vx'] = vx
	structure['vy'] = vy
	structure['vz'] = vz

	if psi_opt == True:
		structure['psi'] = psi
	if phi_opt == True:
		structure['phi'] = phi

	ILCTStructure = namedtuple('ILCTStructure',sorted(structure))

	return ILCTStructure(**structure)
