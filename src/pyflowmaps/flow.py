#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np
from .flowmaker import pyflowmaker
from .math_tools import divergence, odiff
from collections import OrderedDict, namedtuple
from .ilct import pyilct
from .dave4vm import pydave4vm
# Physical/solar constants used in conversions
# Value of kilometers per arcsec (approx. plate scale at 1 AU)
KM_PER_ARCSEC = 725.0
# Mass-flux scale height (November 1989, ApJ, 344, 494)
MASS_FLUX_SCALE_HEIGHT = 150.0

__all__ = ['flowLCT','flowILCT','flowDAVE']
__authors__ = ["Jose Iván Campos Rozo, Santiago Vargas Domínguez"]
__email__ = "hypnus1803@gmail.com"

def flowLCT(mc, fwhm_arcsec, scale, cadence, reb=1,**kwargs):
    """
    Spanish: Programa para generar los mapas de flujos vectoriales vx, vy, vz  en km/s.
    English: Script to generate the vector flow maps vx, vy, vz in km/s.
    Parameters
    ----------
        mc :
            A `~numpy array` of shape ``(nt, ny,nx)`` where ``nt`` is the number of
            layers or images to calculate the flow maps.

        fwhm_arcsec :
            Window size for apodization (arcsec)
        scale :
            Size of pixel (Instument information)
        cadence:
            temporal sampling interval in the time series (seconds between image)

    kwargs flowmaker :
	----------------
	    method : {'square' | 'absolute' | 'corr'}
	    interpolation : {'fivepoint' | 'qfit2' | 'crossD'}
	    window : {'gaussian' | 'boxcar'}
		rebine: the rebinning factor if it is wished.
		lag: he lag between the images to be compared (number of images)
		mu : float, optional
		    Cosine of the heliocentric angle (0 < mu <= 1). When mu < 1 the
		    Gaussian window becomes elliptical to compensate for foreshortening.
		theta : float, optional
		    Position angle (radians, counterclockwise from +x) of the radial
		    direction. Used only when mu < 1.

    Output
    -------
        The function returns the velocity maps for vx, vy, and vz, in km/s. For
        solar images cases, vz = h_m * div(vx,vy). For mor information check
        divergence function inside math_tools.py
    """
    
    structure = OrderedDict()

    fwhm = fwhm_arcsec / scale

    kmperasec = KM_PER_ARCSEC
    h_m = MASS_FLUX_SCALE_HEIGHT

    delta_t = cadence  # time-lag in seconds
    factor = scale * kmperasec / delta_t
    v_limit = 2*reb +reb  # cota maxima velocidad en pixels.

    # ************************************************************
    #               Applying LCT
    # ************************************************************


    HorizontalVelocities = pyflowmaker(mc, fwhm, reb=reb, **kwargs)
    vx = HorizontalVelocities.vx
    vy = HorizontalVelocities.vy

    vx = vx.clip(-v_limit, v_limit)
    vy = vy.clip(-v_limit, v_limit)

    vx_kps = vx.copy() 
    vy_kps = vy.copy()

    vx_kps = vx_kps - np.mean(vx_kps)
    vy_kps = vy_kps - np.mean(vy_kps)


    div = divergence(vx_kps, vy_kps)

    vz_kps = h_m * div

    structure['vx'] = vx_kps*factor
    structure['vy'] = vy_kps*factor
    structure['vz'] = vz_kps

    FlowStructure = namedtuple('FlowStructure', sorted(structure))

    return FlowStructure(**structure)


def flowILCT(vels,BField_comp,fwhm_arcsec, scale,interval,threshold=10,**kwargs):
    """
    Spanish: Programa para generar los mapas de flujos vectoriales vx, vy, vz en km/s.
    English: Script to generate the vector flow maps vx, vy, vz  in km/s.
    Parameters
    ----------
        vels :
            A `~numpy array` of shape ``(2, ny,nx)`` with the horizontal velocity components
            [km/s] obtained from LOS magnetogram data using LCT. [vx,vy].
        BField_comp :
            A `~numpy array` of shape ``(4,ny,nx)`` with B-field components [Gauss]
            whit the time-centered B field (Bx_c,By_c,Bz_c) are stored in BField_comp[0:3,:,:],
            and BField_comp[3,:,:] is the change in Bz such that dBz = Bz_f - Bz_i.
        fwhm_arcsec :
            Window size for apodization (arcsec)
        scale :
            Size of pixel (Instument information)
        interval :
            Time between the first and last image
        threshold :
            According with the literature, the order for the noise magnetic field is in the range(5-15) Gauss.


    kwargs ILCT
    -----------
        psi_opt=False
        phi_opt=False
        check=False
        mask=False


    Results
    -------
        The function returns the velocity maps for vx, vy, and vz in km/s.
    """

    structure = OrderedDict()

    #fwhm = fwhm_arcsec / scale


    kmperasec = KM_PER_ARCSEC  # Value of kilometers per arcsec'

    #factor = scale * kmperasec / delta_t

    # ************************************************************
    #               Applying ILCT
    # ************************************************************


    pix_size = kmperasec * scale * (1.e5)

    vx_cps = vels[0,:,:] * 1.e5  # vx in cm/s
    vy_cps = vels[1,:,:] * 1.e5  # vy in cm/s

    HorizontalField = np.array([vx_cps, vy_cps])

    ILCTField = pyilct(HorizontalField, BField_comp, pix_size, interval,
                       threshold=threshold, **kwargs)

    vels = OrderedDict()
    optional = OrderedDict()

    for i in ILCTField._fields:
        if (i != 'vx') and (i != 'vy') and (i != 'vz'):
            optional[i] = getattr(ILCTField, i)
        else:
            vels[i] = getattr(ILCTField, i)

    vx_kps = vels['vx'] / 1e5
    vy_kps = vels['vy'] / 1e5
    vz_kps = vels['vz'] / 1e5

    vx_kps[np.abs(BField_comp[2,:,:]) < threshold] = np.nan
    vy_kps[np.abs(BField_comp[2,:,:]) < threshold] = np.nan
    vz_kps[np.abs(BField_comp[2,:,:]) < threshold] = np.nan

    vx_kps = vx_kps - np.nanmean(vx_kps)
    vy_kps = vy_kps - np.nanmean(vy_kps)
    vz_kps = vz_kps - np.nanmean(vz_kps)

    vx_kps[np.abs(BField_comp[2, :, :]) < threshold] = 0.0
    vy_kps[np.abs(BField_comp[2, :, :]) < threshold] = 0.0
    vz_kps[np.abs(BField_comp[2, :, :]) < threshold] = 0.0

    structure['vx'] = vx_kps
    structure['vy'] = vy_kps
    structure['vz'] = vz_kps

    if len(optional) > 0:
        structure['OptionalILCT'] = optional

    FlowStructure = namedtuple('FlowStructure', sorted(structure))

    return FlowStructure(**structure)


def flowDAVE(bx_init, by_init, bz_init,
             bx_final, by_final, bz_final,
             pixel_size_arcsec, date, delta_time, window_arcsec):
    """
    Spanish: Genera mapas de flujos horizontales/verticales usando DAVE4VM.
    English: Generate vx, vy, vz in km/s using DAVE4VM between two vector magnetograms.

    Parameters
    ----------
        bx_init, by_init, bz_init :
            Initial B-field components (Gauss), arrays of shape (ny, nx).
        bx_final, by_final, bz_final :
            Final B-field components (Gauss), arrays of shape (ny, nx).
        pixel_size_arcsec :
            Pixel size in arcseconds.
        date :
            Time string (e.g., '2013-06-13') used for Sun-Earth distance if available.
        delta_time :
            Time between final and initial images (seconds).
        window_arcsec :
            DAVE window side length in arcseconds.

    Results
    -------
        Returns a structure with vx, vy, vz in km/s.

        Notes
        -----
        - Pixel scale (km/pixel) is computed via Sun-Earth distance if SunPy is available; otherwise,
            a constant scale `KM_PER_ARCSEC * pixel_size_arcsec` is used.
        - Output velocities are in km/s; optional derivative fields (UX,VX,WX,UY,VY,WY) are local coefficients.
    """

    # Compute pixel size in km using Sun–Earth distance if possible; else constant scale
    try:
        from sunpy.coordinates.sun import earth_distance
        distance_km = earth_distance(date).to('km').value
        pixel_size_km = (np.radians(pixel_size_arcsec / 3600.0) * distance_km)
    except Exception:
        pixel_size_km = KM_PER_ARCSEC * pixel_size_arcsec

    dx = pixel_size_km
    dy = pixel_size_km.math

    window_size = int(window_arcsec / pixel_size_arcsec)
    if window_size % 2 == 0:
        window_size += 1

    bzt = (bz_final - bz_init) / delta_time
    bx = (bx_final + bx_init) / 2.0
    by = (by_final + by_init) / 2.0
    bz = (bz_final + bz_init) / 2.0

    bxx, bxy = odiff(bx)
    byx, byy = odiff(by)
    bzx, bzy = odiff(bz)

    magvm = {
        'bzt': bzt,
        'bx': bx,
        'bxx': bxx / dx,
        'bxy': bxy / dy,
        'by': by,
        'byx': byx / dx,
        'byy': byy / dy,
        'bz': bz,
        'bzx': bzx / dx,
        'bzy': bzy / dy,
        'dx': dx,
        'dy': dy,
        'dt': delta_time,
    }

    vel4vm = pydave4vm(magvm, window_size, dx, dy)

    structure = OrderedDict()
    structure['vx'] = vel4vm['U0']
    structure['vy'] = vel4vm['V0']
    structure['vz'] = vel4vm['W0']
    structure['OptionalDAVE'] = {
        'UX': vel4vm['UX'], 'VX': vel4vm['VX'], 'WX': vel4vm['WX'],
        'UY': vel4vm['UY'], 'VY': vel4vm['VY'], 'WY': vel4vm['WY']
    }

    FlowStructure = namedtuple('FlowStructure', sorted(structure))
    return FlowStructure(**structure)
