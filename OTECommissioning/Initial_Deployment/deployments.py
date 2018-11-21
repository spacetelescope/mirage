"""Generate random deployment errors and apply them to an OPD.

Use
---
    ::
        deployment_errors = get_deployment_errors(reduction_factor=0.2)
        ote = apply_deployment_errors(ote, deployment_errors)
        ote = remove_piston_tip_tilt(ote)
"""

import numpy as np


def get_deployment_errors(reduction_factor=None):
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

    if reduction_factor is not None:
        for key, value in deployment_errors.items():
            deployment_errors[key] = value * reduction_factor

    return deployment_errors


def apply_deployment_errors(ote, deployment_errors):
    """Generate an OPD with the defined deployment errors.
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

    return ote, segment_tilts


def remove_piston_tip_tilt(ote):
    """Remove the piston/tip/tilt from the OPD."""
    ote.remove_piston_tip_tilt = True

    # Also manually zero out the piston/tip/tilt
    for segment in range(18):
        ote.segment_state[segment, :3] = 0
    ote.update_opd()

    return ote
