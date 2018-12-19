"""Generate random deployment errors  and apply them to an OPD.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        deployment_errors = get_deployment_errors(out_dir=out_dir)
        deployment_errors = reduce_deployment_errors(deployment_errors, reduction_factor=0.2, out_dir=out_dir)
        ote, segment_tilts = apply_deployment_errors(ote, deployment_errors, out_dir=out_dir)
        ote = remove_piston_tip_tilt(ote, out_dir=out_dir)
"""

import os
import time
import yaml

from astropy.io import fits
import numpy as np


def get_deployment_errors(save=True, out_dir=None):
    """Randomly generate a expected deployment tolerances.

    Parameters:
    -----------
    save : bool, optional
        Denotes whether to save out the deployment error dictionary as a yaml
    out_dir : str, optional
        Directory in which to store the saved yaml file

    Returns:
    --------
    deployment_errors : dict
        Dictionary containing lists of generated deployment errors

    Notes:
    ------
    Deployment tolerances taken from JWST WFS&C Commissioning and Operations Plan (OTE-24):
    D36168 / 2299462 Rev C Page 10
    """

    deployment_errors = {
        'sm_piston': np.random.normal(loc=0, scale=2500/5),  # microns
        'sm_tilt': np.random.normal(loc=0, scale=1300/5, size=2),  # microradians
        'sm_decenter': np.random.normal(loc=0, scale=2500/5, size=2),  # microns
        'pm_piston': np.random.normal(loc=0, scale=1500/5, size=18),  # microns
        'pm_tilt': np.random.normal(loc=0, scale=1100/5, size=(18, 2)),  # microradians
        'pm_decenter': np.random.normal(loc=0, scale=1300/5, size=(18, 2)),  # microns
        'pm_roc': np.random.normal(loc=0, scale=151/5, size=(18)),  # microns
        'pm_clocking': np.random.normal(loc=0, scale=1200/5, size=(18)),  # microradians
        'global_pm_piston': np.random.normal(loc=0, scale=700/5),  # microns
        'global_pm_tilt': np.random.normal(loc=0, scale=190/5, size=2),  # microradians
        'global_pm_decenter': np.random.normal(loc=0, scale=200/5, size=2),  # microns
        'global_pm_clocking': np.random.normal(loc=0, scale=150/5),  # microradians
    }

    # Save the deployments dictionary to a yaml file that can be opened
    if out_dir is not None and save:
        save_file = os.path.join(out_dir, 'deployment_errors_{}.yaml'.format(time.strftime("%Y%m%d_%H%M%S")))
        with open(save_file, 'w') as f:
            yaml.dump(deployment_errors, f, default_flow_style=False)
        print('Saved deployment errors to {}'.format(save_file))
    elif save:
        raise IOError('Cannot save deployment errors to yaml; no out_dir provided')

    return deployment_errors

def reduce_deployment_errors(deployment_errors, reduction_factor=0.2, save=True, out_dir=None):
    """Reduce an existing dictionary of deployment errors by a given factor.

    Parameters:
    -----------
    deployment_errors : dict
        Dictionary containing lists of deployment errors to reduce
    reduction_factor : optional
        Factor by which to reduce the input deployment errors. Default is 0.2.
    save : bool, optional
        Denotes whether to save out the deployment error dictionary as a yaml
    out_dir : str, optional
        Directory in which to store the saved yaml file

    Returns:
    --------
    deployment_errors : dict
        Dictionary containing lists of reduced deployment errors
    """
    for key, value in deployment_errors.items():
        deployment_errors[key] = value * reduction_factor

    # Save the deployments dictionary to a yaml file that can be opened
    if out_dir is not None and save:
        save_file = os.path.join(out_dir, 'deployment_errors_reduced_{}.yaml'.format(time.strftime("%Y%m%d_%H%M%S")))
        with open(save_file, 'w') as f:
            yaml.dump(deployment_errors, f, default_flow_style=False)
        print('Saved reduced ({}%) deployment errors to {}'.format(reduction_factor * 100, save_file))
    elif save:
        raise IOError('Cannot save deployment errors to yaml; no out_dir provided')

    return deployment_errors


def apply_deployment_errors(ote, deployment_errors, save=True, out_dir=None):
    """Generate an OPD with the defined deployment errors.

    Parameters:
    -----------
    ote : webbpsf.opds.OTE_Linear_Model_WSS object
        Adjustable OTE object upon which to apply the deployment errors
    deployment_errors : dict
        Dictionary containing lists of deployment errors to reduce
    save : bool, optional
        Denotes whether to save out the deployment error dictionary as a yaml
    out_dir : str, optional
        Directory in which to store the saved yaml file

    Returns:
    --------
    ote : webbpsf.opds.OTE_Linear_Model_WSS object
        Adjustable OTE object upon which deployment errors were applied
    segment_tilts : np.ndarray
        List of X & Y tilts for each segment, in microns, to be used later to
        calculate the segment PSF displacement in pixels
    """
    ote.reset()
    ote.remove_piston_tip_tilt = False  # Reset the piston/tip/tilt

    # Add SM moves
    ote.move_sm_local(
        piston=deployment_errors['sm_piston'], xtilt=deployment_errors['sm_tilt'][0],
        ytilt=deployment_errors['sm_tilt'][1], xtrans=deployment_errors['sm_decenter'][0],
        ytrans=deployment_errors['sm_decenter'][1])

    # Add PMSA (segment) moves
    for i, seg in enumerate(ote.segnames[0:18]):
        ote.move_seg_local(
            seg,
            piston=deployment_errors['pm_piston'][i] + deployment_errors['global_pm_piston'],
            xtilt=deployment_errors['pm_tilt'][i][0] + deployment_errors['global_pm_tilt'][0],
            ytilt=deployment_errors['pm_tilt'][i][1] + deployment_errors['global_pm_tilt'][1],
            xtrans=deployment_errors['pm_decenter'][i][0] + deployment_errors['global_pm_decenter'][0],
            ytrans=deployment_errors['pm_decenter'][i][1] + deployment_errors['global_pm_decenter'][1],
            roc=deployment_errors['pm_roc'][i],
            clocking=deployment_errors['pm_clocking'][i] + deployment_errors['global_pm_clocking']
        )

    # Save segment tilt to a list
    segment_tilts = np.copy(ote.segment_state[:, :2])

    if out_dir is not None and save:
        save_file = os.path.join(out_dir, 'OPD_withtilt_{}.fits'.format(time.strftime("%Y%m%d_%H%M%S")))
        hdu = fits.PrimaryHDU(ote.opd, header=ote.opd_header)
        hdu.writeto(save_file)
        print('Saved OPD to {}'.format(save_file))
    elif save:
        raise IOError('Cannot save deployment errors to yaml; no out_dir provided')

    return ote, segment_tilts


def remove_piston_tip_tilt(ote, save=True, out_dir=None):
    """Remove the piston/tip/tilt from the OPD.

    Parameters:
    -----------
    ote : webbpsf.opds.OTE_Linear_Model_WSS object
        Adjustable OTE object from which to remove piston/tip/tilt\
    save : bool, optional
        Denotes whether to save out the deployment error dictionary as a yaml
    out_dir : str, optional
        Directory in which to store the saved yaml file

    Returns:
    --------
    ote : webbpsf.opds.OTE_Linear_Model_WSS object
        Adjustable OTE object without piston/tip/tilt
    """
    ote.remove_piston_tip_tilt = True

    # Also manually zero out the piston/tip/tilt
    for segment in range(18):
        ote.segment_state[segment, :3] = 0
    ote.update_opd()

    if out_dir is not None and save:
        save_file = os.path.join(out_dir, 'OPD_notilt_{}.fits'.format(time.strftime("%Y%m%d_%H%M%S")))
        hdu = fits.PrimaryHDU(ote.opd, header=ote.opd_header)
        hdu.writeto(save_file)
        print('Saved OPD with piston/tip/tilt removed to {}'.format(save_file))
    elif save:
        raise IOError('Cannot save deployment errors to yaml; no out_dir provided')

    return ote
