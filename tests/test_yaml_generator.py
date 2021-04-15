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
# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')

if not ON_GITHUB:
    orig_mirage_data = os.environ['MIRAGE_DATA']

os.environ["MIRAGE_DATA"] = __location__
os.environ["CRDS_PATH"] = os.path.join(__location__, "temp")


def call_with_assorted_user_inputs(catalogs, crs, dates, backgrounds, roll_angles):
    """Test various user-inputs to yaml_generator, including dates, backgrounds,
    roll angles, cosmic rays, and catalogs
    """
    # Make an instance of the SimInput class
    input_xml = os.path.join(__location__, 'test_data/misc/54321/54321_niriss_wfss_prime_nircam_imaging_parallel_mult_targets.xml')
    pointing_file = os.path.join(__location__, 'test_data/misc/54321/54321_niriss_wfss_prime_nircam_imaging_parallel_mult_targets.pointing')
    temp_output_dir = os.path.join(__location__, "temp")

    yam = SimInput(input_xml, pointing_file, verbose=True, catalogs=catalogs,
                   offline=True, output_dir=temp_output_dir,
                   simdata_output_dir=temp_output_dir,
                   datatype='raw', reffile_defaults='crds',
                   cosmic_rays=crs, dates=dates, background=backgrounds,
                   roll_angle=roll_angles)
    yam.use_linearized_darks = True
    yam.create_inputs()
    return yam.info


def test_user_inputs_basic():
    """Check that use-input dates, roll angles, cosmic ray info,
    backgrounds, and catalogs are propagated correctly in cases
    with the more basic user inputs.
    """
    cats = {'TARG1': {'point_source': 'ptsrc1.cat',
                      'galaxy': 'galaxy1.cat',
                      'extended': 'ex1.cat',
                      'moving_pointsource': 'mt_ptsrc1.cat',
                      'moving_sersic': 'mt_gal_1.cat',
                      'moving_extended': 'mt_ext_1.cat',
                      'moving_target_to_track': 'mt_track_1.cat'
                      },
            'TARG2': {'point_source': 'ptsrc2.cat',
                      'galaxy': 'galaxy2.cat',
                      'extended': 'ex2.cat',
                      'moving_pointsource': 'mt_ptsrc2.cat',
                      'moving_sersic': 'mt_gal_2.cat',
                      'moving_extended': 'mt_ext_2.cat',
                      'moving_target_to_track': 'mt_track_2.cat'
                      }
            }
    cr = {'library': 'FLARE', 'scale': 44.0}
    date = '2019-05-25'
    background = {'001': 'high', '002': 'medium', '003': 22.3}
    roll_angle = 34.5

    tab = call_with_assorted_user_inputs(cats, cr, date, background, roll_angle)

    nrc_obs1 = ((np.array(tab['obs_num']) == '001') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs1 = ((np.array(tab['obs_num']) == '001') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs2 = ((np.array(tab['obs_num']) == '002') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs2 = ((np.array(tab['obs_num']) == '002') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs3 = ((np.array(tab['obs_num']) == '003') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs3 = ((np.array(tab['obs_num']) == '003') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs = np.vstack([nrc_obs1, nrc_obs2, nrc_obs3])
    nis_obs = np.vstack([nis_obs1, nis_obs2, nis_obs3])

    # Check dates
    date_to_match = '{} 00:00:00'.format(date)
    date_match = [str(entry) == date_to_match for entry in tab['Date']]

    assert all(date_match) is True

    # Check roll angle
    pav3_match = [entry == str(roll_angle) for entry in tab['PAV3']]
    assert all(pav3_match) is True

    # Check cosmic rays
    cr_lib_match = [entry == cr['library'] for entry in tab['CosmicRayLibrary']]
    assert all(cr_lib_match) is True

    cr_scale_match = [entry == str(cr['scale']) for entry in tab['CosmicRayScale']]
    assert all(cr_scale_match) is True

    # Check backgrounds
    for i, key in enumerate(background):
        back01_nis_match = [entry == str(background[key]) for entry in np.array(tab['BackgroundRate'])[nis_obs[i, :]]]
        assert all(back01_nis_match) is True
        back01_nrcsw_match = [entry == str(background[key]) for entry in np.array(tab['sw_bkgd'])[nrc_obs[i, :]]]
        assert all(back01_nrcsw_match) is True
        back01_nrclw_match = [entry == str(background[key]) for entry in np.array(tab['lw_bkgd'])[nrc_obs[i, :]]]
        assert all(back01_nrclw_match) is True

    # Check catalogs
    cwd = os.getcwd()

    for i, key in enumerate(['TARG1', 'TARG2', 'TARG1']):
        # NIRCAM
        cat01_nrc_ptsrc_match = [entry == os.path.join(cwd, cats[key]['point_source']) for entry in np.array(tab['sw_ptsrc'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_ptsrc_match) is True
        cat01_nrc_galaxy_match = [entry == os.path.join(cwd, cats[key]['galaxy']) for entry in np.array(tab['sw_galcat'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_galaxy_match) is True
        cat01_nrc_ext_match = [entry == os.path.join(cwd, cats[key]['extended']) for entry in np.array(tab['sw_ext'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_ext_match) is True
        cat01_nrc_mt_ptsrc_match = [entry == os.path.join(cwd, cats[key]['moving_pointsource']) for entry in np.array(tab['sw_movptsrc'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_ptsrc_match) is True
        cat01_nrc_mt_gal_match = [entry == os.path.join(cwd, cats[key]['moving_sersic']) for entry in np.array(tab['sw_movgal'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_gal_match) is True
        cat01_nrc_mt_ext_match = [entry == os.path.join(cwd, cats[key]['moving_extended']) for entry in np.array(tab['sw_movext'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_ext_match) is True
        cat01_nrc_mt_track_match = [entry == os.path.join(cwd, cats[key]['moving_target_to_track']) for entry in np.array(tab['sw_solarsys'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_track_match) is True

        # NIRISS
        cat01_nis_ptsrc_match = [entry == os.path.join(cwd, cats[key]['point_source']) for entry in np.array(tab['PointSourceCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_ptsrc_match) is True
        cat01_nis_galaxy_match = [entry == os.path.join(cwd, cats[key]['galaxy']) for entry in np.array(tab['GalaxyCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_galaxy_match) is True
        cat01_nis_ext_match = [entry == os.path.join(cwd, cats[key]['extended']) for entry in np.array(tab['ExtendedCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_ext_match) is True
        cat01_nis_mt_ptsrc_match = [entry == os.path.join(cwd, cats[key]['moving_pointsource']) for entry in np.array(tab['MovingTargetList'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_ptsrc_match) is True
        cat01_nis_mt_gal_match = [entry == os.path.join(cwd, cats[key]['moving_sersic']) for entry in np.array(tab['MovingTargetSersic'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_gal_match) is True
        cat01_nis_mt_ext_match = [entry == os.path.join(cwd, cats[key]['moving_extended']) for entry in np.array(tab['MovingTargetExtended'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_ext_match) is True
        cat01_nis_mt_track_match = [entry == os.path.join(cwd, cats[key]['moving_target_to_track']) for entry in np.array(tab['MovingTargetToTrack'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_track_match) is True

    # Clean up
    temp_output_dir = os.path.join(__location__, "temp")
    os.system('rm -r {}'.format(temp_output_dir))


def test_user_inputs_complex():
    """Check that use-input dates, roll angles, cosmic ray info,
    backgrounds, and catalogs are propagated correctly in cases
    with the more complex user inputs.
    """
    cats = {'TARG1': {'nircam': {'point_source': 'ptsrc_nrc_1.cat',
                                 'galaxy': 'galaxy_nrc_1.cat',
                                 'extended': 'ex_nrc_1.cat',
                                 'moving_pointsource': 'mt_ptsrc_nrc_1.cat',
                                 'moving_sersic': 'mt_gal_nrc_1.cat',
                                 'moving_extended': 'mt_ext_nrc_1.cat',
                                 'moving_target_to_track': 'mt_track_nrc_1.cat'
                                 },
                      'niriss': {'point_source': 'ptsrc_nis_1.cat',
                                 'galaxy': 'galaxy_nis_1.cat',
                                 'extended': 'ex_nis_1.cat',
                                 'moving_pointsource': 'mt_ptsrc_nis_1.cat',
                                 'moving_sersic': 'mt_gal_nis_1.cat',
                                 'moving_extended': 'mt_ext_nis_1.cat',
                                 'moving_target_to_track': 'mt_track_nis_1.cat'
                                 }
                      },
            'TARG2': {'nircam': {'point_source': 'ptsrc_nrc_2.cat',
                                 'galaxy': 'galaxy_nrc_2.cat',
                                 'extended': 'ex_nrc_2.cat',
                                 'moving_pointsource': 'mt_ptsrc_nrc_2.cat',
                                 'moving_sersic': 'mt_gal_nrc_2.cat',
                                 'moving_extended': 'mt_ext_nrc_2.cat',
                                 'moving_target_to_track': 'mt_track_nrc_2.cat'
                                 },
                      'niriss': {'point_source': 'ptsrc_nis_2.cat',
                                 'galaxy': 'galaxy_nis_2.cat',
                                 'extended': 'ex_nis_2.cat',
                                 'moving_pointsource': 'mt_ptsrc_nis_2.cat',
                                 'moving_sersic': 'mt_gal_nis_2.cat',
                                 'moving_extended': 'mt_ext_nis_2.cat',
                                 'moving_target_to_track': 'mt_track_nis_2.cat'
                                 }
                      },
            }

    cr = {'001': {'library': 'FLARE', 'scale': 1.2},
          '002': {'library': 'SUNMIN', 'scale': 5.5},
          '003': {'library': 'SUNMAX', 'scale': 4.4}}

    date = {'001': '2019-05-25', '002': '2019-11-15T12:12:12.120000', '003': '2020-10-14T19:20:21'}

    background = {'001': {'nircam': {'sw': 0.2, 'lw': 0.3}, 'niriss': 0.4},
                  '002': {'nircam': {'sw': 'medium', 'lw': 'high'}, 'niriss': 'low'},
                  '003': {'nircam': {'sw': 0.75, 'lw': 'high'}, 'niriss': 0.2}}

    roll_angle = {'001': 34.5, '002': 154.5, '003': 3.14}

    tab = call_with_assorted_user_inputs(cats, cr, date, background, roll_angle)

    nrc_obs1 = ((np.array(tab['obs_num']) == '001') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs1 = ((np.array(tab['obs_num']) == '001') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs2 = ((np.array(tab['obs_num']) == '002') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs2 = ((np.array(tab['obs_num']) == '002') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs3 = ((np.array(tab['obs_num']) == '003') & (np.array(tab['Instrument']) == 'NIRCAM'))
    nis_obs3 = ((np.array(tab['obs_num']) == '003') & (np.array(tab['Instrument']) == 'NIRISS'))
    nrc_obs = np.vstack([nrc_obs1, nrc_obs2, nrc_obs3])
    nis_obs = np.vstack([nis_obs1, nis_obs2, nis_obs3])

    # Check dates
    for i, key in enumerate(date):
        if 'T' not in date[key]:
            date_to_match = '{} 00:00:00'.format(date[key])
        else:
            date_to_match = date[key].replace('T', ' ')

        date_match = [str(entry) == date_to_match for entry in np.array(tab['Date'])[nis_obs[i, :]]]
        assert all(date_match) is True
        date_match = [str(entry) == date_to_match for entry in np.array(tab['Date'])[nrc_obs[i, :]]]
        assert all(date_match) is True

    # Check roll angle
    for i, key in enumerate(roll_angle):
        pav3_match = [entry == str(roll_angle[key]) for entry in np.array(tab['PAV3'])[nis_obs[i, :]]]
        assert all(pav3_match) is True
        pav3_match = [entry == str(roll_angle[key]) for entry in np.array(tab['PAV3'])[nrc_obs[i, :]]]
        assert all(pav3_match) is True

    # Check cosmic rays
    for i, key in enumerate(cr):
        cr_lib_match = [entry == cr[key]['library'] for entry in np.array(tab['CosmicRayLibrary'])[nis_obs[i, :]]]
        assert all(cr_lib_match) is True
        cr_lib_match = [entry == cr[key]['library'] for entry in np.array(tab['CosmicRayLibrary'])[nrc_obs[i, :]]]
        assert all(cr_lib_match) is True

        cr_scale_match = [entry == str(cr[key]['scale']) for entry in np.array(tab['CosmicRayScale'])[nis_obs[i, :]]]
        assert all(cr_scale_match) is True
        cr_scale_match = [entry == str(cr[key]['scale']) for entry in np.array(tab['CosmicRayScale'])[nrc_obs[i, :]]]
        assert all(cr_scale_match) is True

    # Check backgrounds
    for i, key in enumerate(background):
        back01_nis_match = [entry == str(background[key]['niriss']) for entry in np.array(tab['BackgroundRate'])[nis_obs[i, :]]]
        assert all(back01_nis_match) is True
        back01_nrcsw_match = [entry == str(background[key]['nircam']['sw']) for entry in np.array(tab['sw_bkgd'])[nrc_obs[i, :]]]
        assert all(back01_nrcsw_match) is True
        back01_nrclw_match = [entry == str(background[key]['nircam']['lw']) for entry in np.array(tab['lw_bkgd'])[nrc_obs[i, :]]]
        assert all(back01_nrclw_match) is True

    # Check catalogs
    cwd = os.getcwd()

    for i, key in enumerate(['TARG1', 'TARG2', 'TARG1']):
        # NIRCAM
        cat01_nrc_ptsrc_match = [entry == os.path.join(cwd, cats[key]['nircam']['point_source']) for entry in np.array(tab['sw_ptsrc'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_ptsrc_match) is True
        cat01_nrc_galaxy_match = [entry == os.path.join(cwd, cats[key]['nircam']['galaxy']) for entry in np.array(tab['sw_galcat'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_galaxy_match) is True
        cat01_nrc_ext_match = [entry == os.path.join(cwd, cats[key]['nircam']['extended']) for entry in np.array(tab['sw_ext'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_ext_match) is True
        cat01_nrc_mt_ptsrc_match = [entry == os.path.join(cwd, cats[key]['nircam']['moving_pointsource']) for entry in np.array(tab['sw_movptsrc'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_ptsrc_match) is True
        cat01_nrc_mt_gal_match = [entry == os.path.join(cwd, cats[key]['nircam']['moving_sersic']) for entry in np.array(tab['sw_movgal'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_gal_match) is True
        cat01_nrc_mt_ext_match = [entry == os.path.join(cwd, cats[key]['nircam']['moving_extended']) for entry in np.array(tab['sw_movext'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_ext_match) is True
        cat01_nrc_mt_track_match = [entry == os.path.join(cwd, cats[key]['nircam']['moving_target_to_track']) for entry in np.array(tab['sw_solarsys'])[nrc_obs[i, :]]]
        assert all(cat01_nrc_mt_track_match) is True

        # NIRISS
        cat01_nis_ptsrc_match = [entry == os.path.join(cwd, cats[key]['niriss']['point_source']) for entry in np.array(tab['PointSourceCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_ptsrc_match) is True
        cat01_nis_galaxy_match = [entry == os.path.join(cwd, cats[key]['niriss']['galaxy']) for entry in np.array(tab['GalaxyCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_galaxy_match) is True
        cat01_nis_ext_match = [entry == os.path.join(cwd, cats[key]['niriss']['extended']) for entry in np.array(tab['ExtendedCatalog'])[nis_obs[i, :]]]
        assert all(cat01_nis_ext_match) is True
        cat01_nis_mt_ptsrc_match = [entry == os.path.join(cwd, cats[key]['niriss']['moving_pointsource']) for entry in np.array(tab['MovingTargetList'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_ptsrc_match) is True
        cat01_nis_mt_gal_match = [entry == os.path.join(cwd, cats[key]['niriss']['moving_sersic']) for entry in np.array(tab['MovingTargetSersic'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_gal_match) is True
        cat01_nis_mt_ext_match = [entry == os.path.join(cwd, cats[key]['niriss']['moving_extended']) for entry in np.array(tab['MovingTargetExtended'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_ext_match) is True
        cat01_nis_mt_track_match = [entry == os.path.join(cwd, cats[key]['niriss']['moving_target_to_track']) for entry in np.array(tab['MovingTargetToTrack'])[nis_obs[i, :]]]
        assert all(cat01_nis_mt_track_match) is True

    # Clean up
    temp_output_dir = os.path.join(__location__, "temp")
    os.system('rm -r {}'.format(temp_output_dir))


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

    # Put the entry_number entry of yam.apt_xml_dict into yam.info.
    # Typically this is done when running create_inputs() in the
    # yaml_generator. This step is skipped here to save time and
    # avoid having to deal with lots of output files.
    yam.info['entry_number'] = yam.apt_xml_dict['entry_number']

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
               'program (103), not equal to 99.' in e, \
            'Failed to reject psf_path of incorrect length.'

    # Test for a list of paths of correct length
    list_paths = [__location__] * 9 + [os.path.dirname(__location__)] * 9
    yam.psf_paths = list_paths
    paths_out = yam.get_psf_path()

    assert paths_out == list_paths,\
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
    catalogs = None

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
                                    'transmission': {'nrcb5': {'f322w2': {'clear': 'my_reffiles/my_transmission_for_b5.fits'},
                                                               'f444w': {'clear': 'my_reffiles/my_transmission_for_b5.fits'},
                                                               'f335m': {'clear': 'my_reffiles/my_transmission_for_b5.fits'},
                                                               'f300m': {'clear': 'my_reffiles/my_transmission_for_b5.fits'}},
                                                     'nrcb1': {'f150w': {'clear': 'my_reffiles/my_transmission_for_b1.fits'},
                                                               'f070w': {'clear': 'my_reffiles/my_transmission_for_b1.fits'},
                                                               'f150w2': {'clear': 'my_reffiles/my_transmission_for_b1.fits'},
                                                               'f187n': {'clear': 'my_reffiles/my_transmission_for_b1.fits'}},
                                                     'nrcb2': {'f150w': {'clear': 'my_reffiles/my_transmission_for_b2.fits'},
                                                               'f070w': {'clear': 'my_reffiles/my_transmission_for_b2.fits'},
                                                               'f150w2': {'clear': 'my_reffiles/my_transmission_for_b2.fits'},
                                                               'f187n': {'clear': 'my_reffiles/my_transmission_for_b2.fits'}},
                                                     'nrcb3': {'f150w': {'clear': 'my_reffiles/my_transmission_for_b3.fits'},
                                                               'f070w': {'clear': 'my_reffiles/my_transmission_for_b3.fits'},
                                                               'f150w2': {'clear': 'my_reffiles/my_transmission_for_b3.fits'},
                                                               'f187n': {'clear': 'my_reffiles/my_transmission_for_b3.fits'}},
                                                     'nrcb4': {'f150w': {'clear': 'my_reffiles/my_transmission_for_b4.fits'},
                                                               'f070w': {'clear': 'my_reffiles/my_transmission_for_b4.fits'},
                                                               'f150w2': {'clear': 'my_reffiles/my_transmission_for_b4.fits'},
                                                               'f187n': {'clear': 'my_reffiles/my_transmission_for_b4.fits'}},
                                                    },
                                    'badpixmask': {'nrcb5': 'my_reffiles/my_bpm_for_b5.fits',
                                                   'nrcb4': 'my_reffiles/my_bpm_for_b4.fits'},
                                    'pixelflat': {'nrcb5': {'f322w2': {'clear': 'my_reffiles/my_flatfield_for_b5.fits'}}}
                                    },
                         'niriss': {'superbias': {'nisrapid': 'my_niriss_supebias.fits'},
                                    'linearity': 'my_niriss_linearity,fits',
                                    'saturation': 'my_niriss_saturation.fits',
                                    'gain': 'my_niriss_gain.fits',
                                    'distortion': {'F115W': {'nis_image': 'my_niriss_disotrtion.asdf'}},
                                    'area': {'clear': {'f115w': {'nis_image': 'my_niriss_area.fits'}}},
                                    'transmission': {'clear': {'f115w': 'my_niriss_transmission.fits'},
                                                     'gr150c': {'f115w': 'my_niriss_gr_transmission.fits'}
                                                     },
                                    'badpixmask': 'my_niriss_badpixmask.fits',
                                    'pixelflat': {'clear': {'f115w': 'my_niriss_flatfield.fits'}}
                                    }
                         }

    catalogs = None

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
    match_nrc_lw_flat = np.where((lw_filternames == 'F322W2') & (detectors == 'B5') & (lw_pupilnames == 'CLEAR'))[0]

    match_nis_superbias = np.where(read_patterns == 'NISRAPID')[0]
    match_nis_common = np.where(instruments == 'NIRISS')[0]
    match_nis_distortion_area = np.where((pupilnames == 'F115W') & (instruments == 'NIRISS') & (filternames == 'CLEAR'))[0]
    match_nis_flat = np.where((pupilnames == 'F115W') & (instruments == 'NIRISS') & (filternames == 'CLEAR'))[0]

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
        assert yam.info['transmission'][index] == 'my_reffiles/my_transmission_for_b4.fits'
    for index in match_nrc_lw_distortion_area:
        assert yam.info['astrometric'][index] == 'my_reffiles/my_distortion_for_b5.asdf'
        assert yam.info['pixelAreaMap'][index] == 'my_reffiles/my_pam_for_b5.fits'
        assert yam.info['transmission'][index] == 'my_reffiles/my_transmission_for_b5.fits'
    for index in match_nrc_lw_flat:
        assert yam.info['pixelflat'][index] == 'my_reffiles/my_flatfield_for_b5.fits'

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
        assert yam.info['transmission'][index] == 'my_niriss_transmission.fits'
    for info in match_nis_flat:
        assert yam.info['pixelflat'][index] == 'my_niriss_flatfield.fits'

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


# Return environment variable to original value. This is helpful when
# calling many tests at once, some of which need the real value.
if not ON_GITHUB:
    os.environ['MIRAGE_DATA'] = orig_mirage_data
