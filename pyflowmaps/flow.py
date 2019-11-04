#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np
from flowmaker import pyflowmaker
from math_tools import divergence
from collections import OrderedDict, namedtuple
from ilct import pyilct

__all__ = ['flowLCT','flowILCT']
__authors__ = ["Jose Ivan Campos Rozo, Santiago Vargas Dominguez"]
__email__ = "hypnus1803@gmail.com"

def flowLCT(mc, fwhm_arcsec, scale, cadence,verbose=False, **kwargs):
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
	    method : {'square' | 'absolute' | 'cross'}
	    interpolation : {'fivepoint' | 'qfit2' | 'crossD'}
	    window : {'gaussian' | 'boxcar'}
		rebine: the rebinning factor if it is wished.
		lag: he lag between the images to be compared (number of images)

    Output
    -------
        The function returns the velocity maps for vx, vy, and vz, in km/s. For
        solar images cases, vz = h_m * div(vx,vy). For mor information check
        divergence function inside math_tools.py
    """
    reb = 1
    structure = OrderedDict()
    
    fwhm = fwhm_arcsec / scale
    
    kmperasec = 725  # Value of kilometers per arcsec'
    h_m = 150  # input('mass-flux scale-heigth (November 1989, ApJ,344,494):')
    
    delta_t = cadence  # time-lag in seconds
    factor = scale * kmperasec / delta_t
    v_limit = 2*reb +reb  # cota maxima velocidad en pixels.
    
    # ************************************************************
    #               Applying LCT
    # ************************************************************
    
 
    HorizontalVelocities = pyflowmaker(mc, fwhm, **kwargs)
    vx = HorizontalVelocities.vx
    vy = HorizontalVelocities.vy
    
    vx = vx.clip(-v_limit, v_limit)
    vy = vy.clip(-v_limit, v_limit)
    
    vx_kps = vx * factor  # vx in km/s
    vy_kps = vy * factor  # vy in km/s
    
    vx_kps = vx_kps - np.mean(vx_kps)
    vy_kps = vy_kps - np.mean(vy_kps)
    
    div = divergence(vx_kps, vy_kps)
    
    vz_kps = h_m * div
    
    structure['vx'] = vx_kps
    structure['vy'] = vy_kps
    structure['vz'] = vz_kps

    FlowStructure = namedtuple('FlowStructure', sorted(structure))
    
    return FlowStructure(**structure)


def flowILCT(vels,BField_comp, scale,interval,threshold=10,
             verbose=False, **kwargs):
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


    
    kmperasec = 725  # Value of kilometers per arcsec'
    
    
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
    

    structure['vx'] = vx_kps
    structure['vy'] = vy_kps
    structure['vz'] = vz_kps
    
    if len(optional) > 0:
        structure['OptionalILCT'] = optional
    
    FlowStructure = namedtuple('FlowStructure', sorted(structure))
    
    return FlowStructure(**structure)