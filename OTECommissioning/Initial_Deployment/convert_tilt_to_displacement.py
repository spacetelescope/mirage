"""Functions used to convert a segment tilt in segment control coordinates
into a NIRCam SW pixel displacement in detector coordinates.

Use
---
    ::
        calc_location_after_tilt(segment, xtilt, ytilt)

"""
import numpy as np

# The following values are averages calculated from WebbPSF-generated
# PSFs for each segment. The data for each segment can be found here:
# /user/lchambers/OTECommSims/tilt_conversion_parameters.yaml
TILT_TO_PIXEL_SLOPE = 13.4
TILT_TO_PIXEL_INTERCEPT = 0


def calc_location_after_tilt(segment, xtilt, ytilt):
    """Calculate the location of a segment PSF after a tilt

    Parameters
    ----------
    segment : str
        Segment name, i.e 'C6'
    xtilt : float
        The amount of x tilt applied to the given segment
    ytilt : float
        The amount of x tilt applied to the given segment

    Returns
    -------
    x_displacement
        The shift of the segment PSF in NIRCam SW x pixels
    y_displacement
        The shift of the segment PSF in NIRCam SW y pixels
    """
    tilt_onto_x, tilt_onto_y = _convert_control_to_global(segment, xtilt, ytilt)
    # print('Tilt projected onto V2/V3:', tilt_onto_x, tilt_onto_y)

    x_displacement = tilt_onto_x * TILT_TO_PIXEL_SLOPE
    y_displacement = tilt_onto_y * TILT_TO_PIXEL_SLOPE

    return -x_displacement, -y_displacement


def _convert_control_to_global(segment, xtilt, ytilt):
    """Convert vectors coordinates in the local segment control
    coordinate system to NIRCam detector X and Y coordinates.
    At least proportionally."""
    control_xaxis_rotations = {
        'A1': 180, 'A2': 120, 'A3': 60, 'A4': 0,'A5': -60,
        'A6': -120, 'B1': 0, 'C1': 60, 'B2': -60, 'C2': 0,
        'B3': -120, 'C3': -60, 'B4': -180, 'C4': -120,
        'B5': -240, 'C5': -180, 'B6': -300, 'C6': -240
    }

    x_rot = control_xaxis_rotations[segment[:2]]  # degrees
    x_rot_rad = x_rot * np.pi / 180  # radians
    # print('Rotating by {} deg ({} rad) to account for segment {}'.format(x_rot, x_rot_rad, segment))

    y_det = (xtilt * np.cos(x_rot_rad)) - (ytilt * np.sin(x_rot_rad))
    x_det = (xtilt * np.sin(x_rot_rad)) + (ytilt * np.cos(x_rot_rad))

    return x_det, y_det
