"""Test the functions provided by mirage.yaml.yaml_generator

Authors
-------
    - Lauren Chambers

Use
---
    >>> pytest test_yaml_generator.py
"""
import os

import numpy as np
import pytest

from mirage.yaml.yaml_generator import SimInput

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Reset the MIRAGE_DATA env variable to be a real location so yaml_generator
# doesn't croak
os.environ["MIRAGE_DATA"] = __location__
os.environ["CRDS_PATH"] = os.path.join(__location__, "temp")


def test_get_psf_path():
    """Test that the get_psf_path method of the yaml_generator.SimInput class
    is working as expected.
    """
    # Make an instance of the SimInput class
    input_xml = os.path.join(__location__, 'test_data', 'NIRCam', '1144-OTE-10.xml')
    pointing_file = os.path.join(__location__, 'test_data', 'NIRCam', '1144-OTE-10.pointing')
    yam = SimInput(input_xml, pointing_file, offline=True)

    # Fake the activity IDs and instrument
    n_activities = 101
    act_ids = [np.base_repr(i, 36).zfill(2)  for i in range(n_activities)]
    yam.info = {}
    yam.info['act_id'] = act_ids
    yam.info['Instrument'] = ['NIRCam'] * n_activities

    # Test for a default path
    paths_out = yam.get_psf_path()
    assert len(paths_out) == n_activities,\
        'Default PSF path not properly provided.'
    assert 'nircam/gridded_psf_library' in paths_out[0],\
        'Default PSF path not properly provided.'
    np.testing.assert_array_equal(paths_out, paths_out[0],
                                  err_msg='Default PSF path not properly defined.')

    # Test for a single path
    yam.psf_paths = __location__
    paths_out = yam.get_psf_path()
    assert paths_out == [__location__] * n_activities,\
        'Single PSF path not properly assigned.'

    # Test for a list of paths of incorrect length
    with pytest.raises(ValueError) as e:
        yam.psf_paths = [__location__] * 99
        yam.get_psf_path()
        assert 'Please provide the psf_paths in the form of a list of strings ' \
               'with a length equal to the number of activities in the APT ' \
               'program (101), not equal to 99.' in e, \
            'Failed to reject psf_path of incorrect length.'

    # Test for a list of paths of correct length
    list_101_paths = [__location__] * 50 + [os.path.dirname(__location__)] * 51
    yam.psf_paths = list_101_paths
    paths_out = yam.get_psf_path()
    assert paths_out == sorted(list_101_paths),\
        'List of PSF paths not properly assigned.'

    # Test for a completely invalid path
    with pytest.raises(TypeError) as e:
        yam.psf_paths = 3.054
        yam.get_psf_path()
        assert 'Please provide the psf_paths in the form of a list or ' \
               'string, not float' in e, \
            'Failed to reject psf_path of incorrect type.'


def test_reffile_crds():
    """Make this a copy of test_reffile_crds_full_name (including reffile_overrides)
    but use 'crds'
    """
    catalogs = {'nircam': {'sw': 'sw_catalog.cat', 'lw': 'lw_cartalog.cat'},
                'niriss': 'nis_catalog.cat'}

    # Make an instance of the SimInput class
    input_xml = os.path.join(__location__, 'test_data/misc/12345/12345_nircam_imaging_prime_niriss_wfss_parallel.xml')
    pointing_file = os.path.join(__location__, 'test_data/misc/12345/12345_nircam_imaging_prime_niriss_wfss_parallel.pointing')
    temp_output_dir = os.path.join(__location__, "temp")

    yam = SimInput(input_xml, pointing_file, verbose=True, catalogs=catalogs,
                   offline=True, output_dir=temp_output_dir,
                   simdata_output_dir=temp_output_dir,
                   datatype='raw', reffile_defaults='crds')
    yam.use_linearized_darks = True
    yam.create_inputs()

    # Check that all reference files are simply the string 'crds'
    expected_value = 'crds'
    assert np.all(np.array(yam.info['superbias']) == expected_value)
    assert np.all(np.array(yam.info['linearity']) == expected_value)
    assert np.all(np.array(yam.info['saturation']) == expected_value)
    assert np.all(np.array(yam.info['gain']) == expected_value)
    assert np.all(np.array(yam.info['astrometric']) == expected_value)
    assert np.all(np.array(yam.info['pixelAreaMap']) == expected_value)
    assert np.all(np.array(yam.info['badpixmask']) == expected_value)
    assert np.all(np.array(yam.info['ipc']) == expected_value)

# Clean up
    os.system('rm -r {}'.format(temp_output_dir))


def test_reffile_crds_full_name():
    """Test that the correct values for reference files are found
    """
    reffile_overrides = {'nircam': {'superbias': {'nrcb5': {'bright1': 'my_reffiles/my_superbias_for_b5.fits',
                                                            'shallow4': 'my_reffiles/my_superbias_for_b5.fits'
                                                            },
                                                  'nrcb4': {'shallow2': 'my_reffiles/my_superbias_for_b4.fits'}
                                                  },
                                    'linearity': {'nrcb5': 'my_reffiles/my_linearity_for_b5.fits',
                                                  'nrcb4': 'my_reffiles/my_linearity_for_b4.fits'},
                                    'saturation': {'nrcb5': 'my_reffiles/my_saturation_for_b5.fits',
                                                   'nrcb4': 'my_reffiles/my_saturation_for_b4.fits'},
                                    'gain': {'nrcb5': 'my_reffiles/my_gain_for_b5.fits',
                                             'nrcb4': 'my_reffiles/my_gain_for_b4.fits'},
                                    'distortion': {'nrcb5': {'f322w2': {'NRC_IMAGE': 'my_reffiles/my_distortion_for_b5.asdf'}},
                                                   'nrcb4': {'f150w': {'NRC_IMAGE': 'my_reffiles/my_distortion_for_b4.asdf'}}},
                                    'area': {'nrcb5': {'f322w2': {'clear': {'nrc_image': 'my_reffiles/my_pam_for_b5.fits'}}},
                                             'nrcb4': {'f150w': {'clear': {'nrc_image': 'my_reffiles/my_pam_for_b4.fits'}}}},
                                    'badpixmask': {'nrcb5': 'my_reffiles/my_bpm_for_b5.fits',
                                                   'nrcb4': 'my_reffiles/my_bpm_for_b4.fits'},
                                    },
                         'niriss': {'superbias': {'nisrapid': 'my_niriss_supebias.fits'},
                                    'linearity': 'my_niriss_linearity,fits',
                                    'saturation': 'my_niriss_saturation.fits',
                                    'gain': 'my_niriss_gain.fits',
                                    'distortion': {'F115W': {'nis_image': 'my_niriss_disotrtion.asdf'}},
                                    'area': {'clear': {'f115w': {'nis_image': 'my_niriss_area.fits'}}},
                                    'badpixmask': 'my_niriss_badpixmask.fits'
                                    }
                         }

    catalogs = {'nircam': {'sw': 'sw_catalog.cat', 'lw': 'lw_cartalog.cat'},
                'niriss': 'nis_catalog.cat'}

    # Make an instance of the SimInput class
    input_xml = os.path.join(__location__, 'test_data/misc/12345/12345_nircam_imaging_prime_niriss_wfss_parallel.xml')
    pointing_file = os.path.join(__location__, 'test_data/misc/12345/12345_nircam_imaging_prime_niriss_wfss_parallel.pointing')
    temp_output_dir = os.path.join(__location__, "temp")

    yam = SimInput(input_xml, pointing_file, verbose=True, catalogs=catalogs,
                   offline=True, output_dir=temp_output_dir,
                   simdata_output_dir=temp_output_dir,
                   datatype='raw', reffile_defaults='crds_full_name',
                   reffile_overrides=reffile_overrides)
    yam.use_linearized_darks = True
    yam.create_inputs()

    # Make into numpy arrays for easier matching
    sw_filternames = np.array(yam.info['ShortFilter'])
    lw_filternames = np.array(yam.info['LongFilter'])
    sw_pupilnames = np.array(yam.info['ShortPupil'])
    lw_pupilnames = np.array(yam.info['LongPupil'])
    filternames = np.array(yam.info['FilterWheel'])
    pupilnames = np.array(yam.info['PupilWheel'])
    detectors = np.array(yam.info['detector'])
    instruments = np.array(yam.info['Instrument'])
    read_patterns = np.array(yam.info['ReadoutPattern'])

    match_nrc_sw_superbias = np.where((read_patterns == 'SHALLOW2') & (detectors == 'B4'))[0]
    match_nrc_lw_superbias = np.where((read_patterns == 'BRIGHT1') & (detectors == 'B5'))[0]
    match_nrc_sw_common = np.where(detectors == 'B4')[0]
    match_nrc_lw_common = np.where(detectors == 'B5')[0]
    match_nrc_sw_distortion_area = np.where((sw_filternames == 'F150W') & (detectors == 'B4') & (sw_pupilnames == 'CLEAR'))[0]
    match_nrc_lw_distortion_area = np.where((lw_filternames == 'F322W2') & (detectors == 'B5') & (lw_pupilnames == 'CLEAR'))[0]

    match_nis_superbias = np.where(read_patterns == 'NISRAPID')[0]
    match_nis_common = np.where(instruments == 'NIRISS')[0]
    match_nis_distortion_area = np.where((pupilnames == 'F115W') & (instruments == 'NIRISS') & (filternames == 'CLEAR'))[0]

    # Check that reference files that are covered by reffile_overrides
    # are equal to the override values.
    for index in match_nrc_sw_superbias:
        assert yam.info['superbias'][index] == 'my_reffiles/my_superbias_for_b4.fits'
    for index in match_nrc_lw_superbias:
        assert yam.info['superbias'][index] == 'my_reffiles/my_superbias_for_b5.fits'
    for index in match_nrc_sw_common:
        assert yam.info['linearity'][index] == 'my_reffiles/my_linearity_for_b4.fits'
        assert yam.info['saturation'][index] == 'my_reffiles/my_saturation_for_b4.fits'
        assert yam.info['gain'][index] == 'my_reffiles/my_gain_for_b4.fits'
        assert yam.info['badpixmask'][index] == 'my_reffiles/my_bpm_for_b4.fits'
    for index in match_nrc_lw_common:
        assert yam.info['linearity'][index] == 'my_reffiles/my_linearity_for_b5.fits'
        assert yam.info['saturation'][index] == 'my_reffiles/my_saturation_for_b5.fits'
        assert yam.info['gain'][index] == 'my_reffiles/my_gain_for_b5.fits'
        assert yam.info['badpixmask'][index] == 'my_reffiles/my_bpm_for_b5.fits'
    for index in match_nrc_sw_distortion_area:
        assert yam.info['astrometric'][index] == 'my_reffiles/my_distortion_for_b4.asdf'
        assert yam.info['pixelAreaMap'][index] == 'my_reffiles/my_pam_for_b4.fits'
    for index in match_nrc_lw_distortion_area:
        assert yam.info['astrometric'][index] == 'my_reffiles/my_distortion_for_b5.asdf'
        assert yam.info['pixelAreaMap'][index] == 'my_reffiles/my_pam_for_b5.fits'

    for index in match_nis_superbias:
        assert yam.info['superbias'][index] == 'my_niriss_supebias.fits'
    for index in match_nis_common:
        assert yam.info['linearity'][index] == 'my_niriss_linearity,fits'
        assert yam.info['saturation'][index] == 'my_niriss_saturation.fits'
        assert yam.info['gain'][index] == 'my_niriss_gain.fits'
        assert yam.info['badpixmask'][index] == 'my_niriss_badpixmask.fits'
    for info in match_nis_distortion_area:
        assert yam.info['astrometric'][index] == 'my_niriss_disotrtion.asdf'
        assert yam.info['pixelAreaMap'][index] == 'my_niriss_area.fits'

    # Check that reference files covered by reffile_overrides contain
    # the CRDS_PATH, which here has been set to the temp directory
    match_nrc_sw_defaults = np.where((instruments == 'NIRCAM') & (sw_filternames == 'F070W') & (detectors != 'B5'))[0]
    match_nrc_lw_defaults = np.where((instruments == 'NIRCAM') & (lw_filternames == 'F335M') & (detectors == 'B5'))[0]

    for index in list(match_nrc_sw_defaults) + list(match_nrc_lw_defaults):
        assert temp_output_dir in yam.info['superbias'][index]
        assert temp_output_dir in yam.info['astrometric'][index]
        assert temp_output_dir in yam.info['pixelAreaMap'][index]
        assert temp_output_dir in yam.info['ipc'][index]

    # Clean up
    os.system('rm -r {}'.format(temp_output_dir))
