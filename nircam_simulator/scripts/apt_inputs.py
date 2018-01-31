# ! /usr/bin/env python

'''
Given APT output files, read in data relevant to the data simulator,
organize, and create input files for the simulator.

Inputs:

xml file - Name of xml file exported from APT.
pointing file - Name of associated pointing file exported from APT.
siaf - Name of csv version of SIAF.

Optional inputs:

epoch_list - Name of ascii file which lists observation labels and
             associated starting observation time, as well as telescope
             roll angle (PAV3). If you wish to have observations back-
             to-back, give them all the same starting time.

Outputs:

output_csv - Ascii table containing all relevant information needed
             to construct a ramp simulator input file for each
             observation/dither/detector combination.

Dependencies:

argparse, lxml, astropy, numpy, collections

JWST Calibration pipeline (only for the use of set_telescope_pointing.py)
in order to translate PAV3 values to local roll angles for each detector.

rotations.py - Colin Cox's module of WCS-related functions


HISTORY:

July 2017 - V0: Initial version. Bryan Hilbert

'''
import os
import sys
import re
import collections
import argparse
from lxml import etree
from astropy.table import Table, Column
from astropy.io import ascii
import numpy as np
import yaml
import pprint
from . import rotations
from . import set_telescope_pointing_separated as set_telescope_pointing

class AptInput:
    """Summary

    Attributes:
        exposure_tab (TYPE): Description
        input_xml (str): Description
        observation_table (str): Description
        obstab (TYPE): Description
        output_csv (TYPE): Description
        pointing_file (str): Description
        siaf (str): Description
    """

    def __init__(self):
        self.input_xml = ''  # e.g. 'GOODSS_ditheredDatasetTest.xml'
        self.output_csv = None  # e.g. 'GOODSS_ditheredDatasetTest.csv'
        self.pointing_file = '' #  e.g. 'GOODSS_ditheredDatasetTest.pointing'
        self.siaf = ''
        self.observation_table = ''

    def create_input_table(self):
        # Expand paths to full paths
        self.input_xml = os.path.abspath(self.input_xml)
        self.pointing_file = os.path.abspath(self.pointing_file)
        self.siaf = os.path.abspath(self.siaf)
        if self.output_csv is not None:
            self.output_csv = os.path.abspath(self.output_csv)
        if self.observation_table is not None:
            self.observation_table = os.path.abspath(self.observation_table)

        main_dir = os.path.split(self.input_xml)[0]

        # Read in xml file
        tab = self.read_xml(self.input_xml)

        # ascii.write(Table(tab), 'as_read_in.csv', format='csv', overwrite=True)

        # Expand the dictionary for multiple dithers. Expand such that there
        # is one entry in each list for each exposure, rather than one entry
        # for each set of dithers
        xmltab = self.expand_for_dithers(tab)

        # ascii.write(Table(xmltab), 'expand_for_dithers.csv', format='csv', overwrite=True)

        # read in the pointing file and produce dictionary
        pointing_tab = self.get_pointing_info(self.pointing_file, xmltab['ProposalID'][0])

        # combine the dictionaries
        obstab = self.combine_dicts(xmltab, pointing_tab)

        # ascii.write(Table(obstab), 'add_pointing_info.csv', format='csv', overwrite=True)

        # add epoch information
        # obstab = self.add_epochs(obstab)

        # add epoch and catalog information
        obstab = self.add_observation_info(obstab)

        # Expand for detectors. Create one entry in each list for each
        # detector, rather than a single entry for 'ALL' or 'BSALL'
        self.exposure_tab = self.expand_for_detectors(obstab)

        detectors_file = os.path.join(main_dir, 'expand_for_detectors.csv')
        ascii.write(Table(self.exposure_tab), detectors_file, format='csv', overwrite=True)
        print('Wrote exposure table to {}'.format(detectors_file))

        # Calculate the correct V2, V3 and RA, Dec for each exposure/detector
        self.ra_dec_update()

        # Output to a csv file.
        if self.output_csv is None:
            indir, infile = os.path.split(self.input_xml)
            self.output_csv = os.path.join(indir, 'Observation_table_for_' + infile + '.csv')
        ascii.write(Table(self.exposure_tab), self.output_csv, format='csv', overwrite=True)
        print('Final csv exposure list written to {}'.format(self.output_csv))

    def read_xml(self, infile):
        """Read in the .xml file from APT. Can read templates for NircamImaging,
        NircamEngineeringImaging, and NircamWfss modes.

        Arguments
        =========
        infile (str):
            Path to input .xml file

        Returns
        =======
        dict:
            Dictionary with extracted observation parameters

        Raises
        ======
        ValueError:
            If an .xml file is provided that includes an APT template that is not
            supported
        """

        # Open XML file, get element tree of the APT proposal
        with open(infile) as f:
            tree = etree.parse(f)

        # Define the needed namespaces
        apt = '{http://www.stsci.edu/JWST/APT}'
        ncei = "{http://www.stsci.edu/JWST/APT/Template/NircamEngineeringImaging}"
        nci = "{http://www.stsci.edu/JWST/APT/Template/NircamImaging}"
        ncwfss = "{http://www.stsci.edu/JWST/APT/Template/NircamWfss}"
        wfscc = "{http://www.stsci.edu/JWST/APT/Template/WfscCommissioning}"
        wfscga = "{http://www.stsci.edu/JWST/APT/Template/WfscGlobalAlignment}"

        # Set up dictionary of observation parameters to be populated
        ProposalParams_keys = ['PI_Name', 'Proposal_category', 'ProposalID',
                               'Science_category', 'Title']
        ObsParams_keys = ['Module', 'Subarray',
                          'PrimaryDitherType', 'PrimaryDithers', 'SubpixelPositions',
                          'SubpixelDitherType', 'CoordinatedParallel',
                          'ObservationID', 'TileNumber', 'APTTemplate']
        FilterParams_keys = ['ShortFilter', 'LongFilter', 'ShortPupil', 'LongPupil',
                             'ReadoutPattern', 'Groups', 'Integrations']
        OtherParams_keys = ['Mode', 'Grism']

        APTObservationParams_keys = ProposalParams_keys + ObsParams_keys + \
                                    FilterParams_keys + OtherParams_keys

        APTObservationParams = {}
        for key in APTObservationParams_keys:
            APTObservationParams[key] = []

        # Get high-level information: proposal info - - - - - - - - - - - - - -

        # Set default values
        propid_default = 42424
        proptitle_default = 'Looking for my towel'
        scicat_default = 'Planets and Planet Formation'
        piname_default = 'D.N. Adams'
        propcat_default = 'GO'

        # Get just the element with the proposal information
        proposal_info = tree.find(apt + 'ProposalInformation')

        # Title
        try:
            prop_title = proposal_info.find(apt + 'Title').text
        except:
            prop_title = proptitle_default

        # Proposal ID
        try:
            prop_id = proposal_info.find(apt + 'ProposalID').text
        except:
            prop_id = propid_default

        # Proposal Category
        try:
            prop_category = proposal_info.find(apt + 'ProposalCategory')[0]
            prop_category = etree.QName(prop_category).localname
        except:
            prop_category = propcat_default

        # Science Category
        try:
            science_category = proposal_info.find(apt + 'ScientificCategory').text
        except:
            science_category = scicat_default

        # Principal Investigator Name
        try:
            pi_firstname = proposal_info.find('.//' + apt + 'FirstName').text
            pi_lastname = proposal_info.find('.//' + apt + 'LastName').text
            pi_name = ' '.join([pi_firstname, pi_lastname])
        except:
            pi_name = piname_default

        # Get parameters for each observation  - - - - - - - - - - - - - - - -

        # Find all observations (but use only those that use NIRCam or are WFSC)
        observation_data = tree.find(apt + 'DataRequests')
        obs_results = observation_data.findall('.//' + apt + 'Observation')

        observations = []
        i_observations = []
        obs_indices = range(len(obs_results))
        for o, i_obs in zip(obs_results, obs_indices):
            if o.find(apt + 'Instrument').text in ['NIRCAM', 'WFSC']:
                observations.append(o)
                i_observations.append(i_obs)

        # Get parameters out!
        for i_obs, obs in zip(i_observations, observations):

            # Determine what template is used for the observation
            template = obs.find(apt + 'Template')[0]
            template_name = etree.QName(template).localname

            # Are all the templates in the XML file something that we can handle?
            known_APT_templates = ['NircamImaging', 'NircamWfss', 'WfscCommissioning',
                                   'NircamEngineeringImaging', 'WfscGlobalAlignment']
            if template_name not in known_APT_templates:
                # If not, turn back now.
                raise ValueError('No protocol written to read {} template.'.format(template_name))

            obs_tuple_list = []

            # Get observation label
            label_ele = obs.find(apt + 'Label')
            if label_ele is not None:
                label = label_ele.text
                if (' (' in label) and (')' in label):
                    label = re.split(r' \(|\)', label)[0]

            else:
                label = 'None'

            # Get coordinated parallel (?)
            coordparallel = obs.find(apt + 'CoordinatedParallel').text

            # Determine pointing offset?
            offset = obs.find('.//' + apt + 'Offset')
            try:
                offset_x = offset.get('Xvalue')
                offset_y = offset.get('Yvalue')
            except AttributeError:
                offset_x, offset_y = 0, 0

            if (offset_x != 0) or (offset_y != 0):
                print('* * * OFFSET OF ({}, {}) IN OBS {} NOT APPLIED ***'.format(offset_x,
                                                                                  offset_y,
                                                                                  i_obs + 1))

            # If template is NircamImaging or NircamEngineeringImaging
            if template_name in ['NircamImaging', 'NircamEngineeringImaging']:
                # Set namespace
                if template_name == 'NircamImaging':
                    ns = nci
                elif template_name == 'NircamEngineeringImaging':
                    ns = ncei

                # Set parameters that are constant for all imaging obs
                typeflag = template_name
                grismval = 'N/A'
                short_pupil = 'CLEAR'

                # Find observation-specific parameters
                mod = template.find(ns + 'Module').text
                subarr = template.find(ns + 'Subarray').text
                pdithtype = template.find(ns + 'PrimaryDitherType').text

                # Determine if there is an aperture override
                override = obs.find('.//' + apt + 'FiducialPointOverride')
                if override is not None:
                    mod = override.text
                    if 'FULL' not in mod:
                        config = ascii.read('../config/NIRCam_subarray_definitions.list')
                        try:
                            i_sub = list(config['AperName']).index(mod)
                        except ValueError:
                            i_sub = i_sub = [mod in name for name in np.array(config['AperName'])]
                            i_sub = np.where(i_sub)[0]
                            if len(i_sub) > 1:
                                raise ValueError('Unable to match \
                                    FiducialPointOverride {} to valid \
                                    aperture.'.format(mod))

                        subarr = config['Name'][i_sub][0]
                        print('Aperture override: subarray {}'.format(subarr[0]))

                try:
                    pdither = template.find(ns + 'PrimaryDithers').text
                except:
                    pdither = '1'

                sdithtype = template.find(ns + 'SubpixelDitherType').text

                try:
                    sdither = template.find(ns + 'SubpixelPositions').text
                except:
                    try:
                        stemp = template.find(ns + 'CoordinatedParallelSubpixelPositions').text
                        sdither = np.int(stemp[0])
                    except:
                        sdither = '1'

                # Find filter parameters for all filter configurations within obs
                filter_configs = template.findall('.//' + ns + 'FilterConfig')

                for filt in filter_configs:
                    sfilt = filt.find(ns + 'ShortFilter').text
                    lfilt = filt.find(ns + 'LongFilter').text
                    rpatt = filt.find(ns + 'ReadoutPattern').text
                    grps = filt.find(ns + 'Groups').text
                    ints = filt.find(ns + 'Integrations').text

                    # Separate pupil and filter in case of filter that is
                    # mounted in the pupil wheel
                    if ' + ' in sfilt:
                        split_ind = sfilt.find(' + ')
                        short_pupil = sfilt[0:split_ind]
                        sfilt = sfilt[split_ind + 1:]
                    else:
                        short_pupil = 'CLEAR'

                    if ' + ' in lfilt:
                        p = lfilt.find(' + ')
                        long_pupil = lfilt[0:p]
                        lfilt = lfilt[p + 1:]
                    else:
                        long_pupil = 'CLEAR'

                    # Add all parameters to dictionary
                    tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                  science_category, typeflag, mod, subarr, pdithtype,
                                  pdither, sdithtype, sdither, sfilt, lfilt,
                                  rpatt, grps, ints, short_pupil,
                                  long_pupil, grismval, coordparallel,
                                  i_obs + 1, 1, template_name)
                    APTObservationParams = self.add_exposure(APTObservationParams, tup_to_add)
                    obs_tuple_list.append(tup_to_add)

            # If template is WFSC Commissioning
            if template_name in ['WfscCommissioning']:
                # Set namespace
                if template_name == 'WfscCommissioning':
                    ns = wfscc
                # elif template_name == 'NircamEngineeringImaging':
                #     ns = ncei

                # Set parameters that are constant for all WFSC obs
                typeflag = template_name
                grismval = 'N/A'
                short_pupil = 'CLEAR'
                subarr = 'FULL'
                pdithtype = 'NONE'
                pdither = '1'
                sdithtype = 'STANDARD'
                sdither = '1'

                # Find observation-specific parameters
                mod = template.find(ns + 'Module').text
                num_WFCgroups = int(template.find(ns + 'ExpectedWfcGroups').text)

                # Determine if there is an aperture override
                override = obs.find('.//' + apt + 'FiducialPointOverride')
                if override is not None:
                    mod = override.text
                    if 'FULL' not in mod:
                        config = ascii.read('../config/NIRCam_subarray_definitions.list')
                        try:
                            i_sub = list(config['AperName']).index(mod)
                        except ValueError:
                            i_sub = i_sub = [mod in name for name in  np.array(config['AperName'])]
                            i_sub = np.where(i_sub)[0]
                            if len(i_sub) > 1:
                                raise ValueError('Unable to match \
                                    FiducialPointOverride {} to valid \
                                    aperture.'.format(mod))

                        subarr = config['Name'][i_sub][0]
                        print('Aperture override: subarray {}'.format(subarr[0]))

                # Find filter parameters for all filter configurations within obs
                filter_configs = template.findall('.//' + ns + 'FilterConfig')

                for filt in filter_configs:
                    sfilt = filt.find(ns + 'ShortFilter').text
                    lfilt = filt.find(ns + 'LongFilter').text
                    rpatt = filt.find(ns + 'ReadoutPattern').text
                    grps = filt.find(ns + 'Groups').text
                    ints = filt.find(ns + 'Integrations').text

                    # Separate pupil and filter in case of filter that is
                    # mounted in the pupil wheel
                    if ' + ' in sfilt:
                        split_ind = sfilt.find(' + ')
                        short_pupil = sfilt[0:split_ind]
                        sfilt = sfilt[split_ind + 1:]
                    else:
                        short_pupil = 'CLEAR'

                    if ' + ' in lfilt:
                        p = lfilt.find(' + ')
                        long_pupil = lfilt[0:p]
                        lfilt = lfilt[p + 1:]
                    else:
                        long_pupil = 'CLEAR'

                    # Repeat for the number of expected WFSC groups + 1
                    for j in range(num_WFCgroups + 1):
                        # Add all parameters to dictionary
                        tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                      science_category, typeflag, mod, subarr, pdithtype,
                                      pdither, sdithtype, sdither, sfilt, lfilt,
                                      rpatt, grps, ints, short_pupil,
                                      long_pupil, grismval, coordparallel,
                                      i_obs + 1, j + 1, template_name)

                        APTObservationParams = self.add_exposure(APTObservationParams, tup_to_add)
                        obs_tuple_list.append(tup_to_add)

            # If template is WFSC Global Alignment
            if template_name in ['WfscGlobalAlignment']:
                ns = wfscga

                # Set parameters that are constant for all WFSC obs
                typeflag = template_name
                grismval = 'N/A'
                short_pupil = 'CLEAR'
                subarr = 'FULL'
                pdither = '1'
                pdithtype = 'NONE'
                sdithtype = 'STANDARD'
                sdither = '1'

                # Determine the Global Alignment Iteration Type
                GA_iteration = obs.find('.//' + wfscga + 'GaIteration').text

                if GA_iteration == 'ADJUST1':
                    n_exp = 3
                elif GA_iteration == 'ADJUST2':
                    n_exp = 6  # technically 5, but 3 is repeated?
                elif GA_iteration == 'BSCORRECT':
                    # Technically has 2 dithers, but that doesn't seem to be incorporated...
                    n_exp = 2
                elif GA_iteration == 'CORRECT+ADJUST':
                    n_exp = 6  # technically 5, but 3 is repeated?
                elif GA_iteration == 'CORRECT':
                    n_exp = 3

                # Find observation-specific parameters
                mod = template.find(ns + 'Module').text
                # num_WFCgroups = int(template.find(ns + 'ExpectedWfcGroups').text)

                # Determine if there is an aperture override
                override = obs.find('.//' + apt + 'FiducialPointOverride')
                if override is not None:
                    mod = override.text
                    if 'FULL' not in mod:
                        config = ascii.read('../config/NIRCam_subarray_definitions.list')
                        try:
                            i_sub = list(config['AperName']).index(mod)
                        except ValueError:
                            i_sub = i_sub = [mod in name for name in  np.array(config['AperName'])]
                            i_sub = np.where(i_sub)[0]
                            if len(i_sub) > 1:
                                raise ValueError('Unable to match \
                                    FiducialPointOverride {} to valid \
                                    aperture.'.format(mod))

                        subarr = config['Name'][i_sub][0]
                        print('Aperture override: subarray {}'.format(subarr[0]))

                # Find filter parameters for all filter configurations within obs
                ga_nircam_configs = template.findall('.//' + ns + 'NircamParameters')

                for conf in ga_nircam_configs:
                    sfilt = conf.find(ns + 'ShortFilter').text
                    lfilt = conf.find(ns + 'LongFilter').text
                    rpatt = conf.find(ns + 'ReadoutPattern').text
                    grps = conf.find(ns + 'Groups').text
                    ints = conf.find(ns + 'Integrations').text

                    # Separate pupil and filter in case of filter that is
                    # mounted in the pupil wheel
                    if ' + ' in sfilt:
                        split_ind = sfilt.find(' + ')
                        short_pupil = sfilt[0:split_ind]
                        sfilt = sfilt[split_ind + 1:]
                    else:
                        short_pupil = 'CLEAR'

                    if ' + ' in lfilt:
                        p = lfilt.find(' + ')
                        long_pupil = lfilt[0:p]
                        lfilt = lfilt[p + 1:]
                    else:
                        long_pupil = 'CLEAR'

                # Repeat for the number of exposures + 1
                for j in range(n_exp + 1):
                    # Add all parameters to dictionary
                    tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                  science_category, typeflag, mod, subarr, pdithtype,
                                  pdither, sdithtype, sdither, sfilt, lfilt,
                                  rpatt, grps, ints, short_pupil,
                                  long_pupil, grismval, coordparallel,
                                  i_obs + 1, j + 1, template_name)

                    APTObservationParams = self.add_exposure(APTObservationParams, tup_to_add)
                    obs_tuple_list.append(tup_to_add)

            # If template is WFSS
            if template_name == 'NircamWfss':
                # Set namespace
                ns = ncwfss

                mod = template.find(ns + 'Module').text
                subarr = template.find(ns + 'Subarray').text
                grismval = template.find(ns + 'Grism').text
                if grismval == 'BOTH':
                    grismval = ['GRISMR', 'GRISMC']
                else:
                    grismval = [grismval]
                # pdithtype = template.find(ns + 'PrimaryDitherType').text
                # pdither = template.find(ns + 'PrimaryDithers').text
                # sdither = template.find(ns + 'SubpixelPositions').text
                # sdithtype = template.find(ns + 'SubpixelPositions').text
                explist = template.find(ns + 'ExposureList')
                expseqs = explist.findall(ns + 'ExposureSequences')

                # if BOTH was specified for the grism,
                # then we need to repeat the sequence of
                # grism/direct/grism/direct/outoffield for each grism
                for gnum in range(len(grismval)):
                    for expseq in expseqs:
                        # sequence = grism, direct, grism, direct, outoffield
                        # if grism == both, sequence is done for grismr,
                        # then repeated for grismc
                        grismvalue = grismval[gnum]
                        # need to switch the order of the grism and direct
                        # exposures in order for them to be chronological
                        grismexp = expseq.find(ns + 'GrismExposure')
                        typeflag = 'WFSS'
                        sfilt = grismexp.find(ns + 'ShortFilter').text
                        lfilt = grismexp.find(ns + 'LongFilter').text
                        rpatt = grismexp.find(ns + 'ReadoutPattern').text
                        groups = grismexp.find(ns + 'Groups').text
                        integrations = grismexp.find(ns + 'Integrations').text

                        pdithtype = template.find(ns + 'PrimaryDitherType').text
                        pdither = template.find(ns + 'PrimaryDithers').text
                        sdither = template.find(ns + 'SubpixelPositions').text
                        sdithtype = template.find(ns + 'SubpixelPositions').text

                        # separate pupil and filter in case of filter
                        # that is mounted in the pupil wheel
                        if ' + ' in sfilt:
                            p = sfilt.find(' + ')
                            short_pupil = sfilt[0:p]
                            sfilt = sfilt[p + 1:]
                        else:
                            short_pupil = 'CLEAR'

                        long_pupil = grismvalue
                        tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                      science_category, typeflag, mod, subarr,
                                      pdithtype, pdither, sdithtype,
                                      sdither, sfilt, lfilt, rpatt, groups,
                                      integrations, short_pupil, long_pupil,
                                      grismvalue, coordparallel,
                                      i_obs + 1, 1, template_name)

                        APTObservationParams = self.add_exposure(APTObservationParams, tup_to_add)
                        obs_tuple_list.append(tup_to_add)

                        directexp = expseq.find(ns + 'DiExposure')
                        typeflag = template_name
                        pdither = '1'  # direct image has no dithers
                        sdither = '1'  # direct image has no dithers
                        sdithtype = '1'  # direct image has no dithers
                        grismvalue = 'N/A'
                        sfilt = directexp.find(ns + 'ShortFilter').text
                        lfilt = directexp.find(ns + 'LongFilter').text
                        rpatt = directexp.find(ns + 'ReadoutPattern').text
                        grps = directexp.find(ns + 'Groups').text
                        ints = directexp.find(ns + 'Integrations').text

                        # separate pupil and filter in case of filter
                        # that is mounted in the pupil wheel
                        if ' + ' in sfilt:
                            p = sfilt.find(' + ')
                            short_pupil = sfilt[0:p]
                            sfilt = sfilt[p + 1:]
                        else:
                            short_pupil = 'CLEAR'

                        if ' + ' in lfilt:
                            p = lfilt.find(' + ')
                            long_pupil = lfilt[0:p]
                            lfilt = lfilt[p + 1:]
                        else:
                            long_pupil = 'CLEAR'

                        direct_tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                             science_category, typeflag, mod, subarr, pdithtype,
                                             pdither, sdithtype, sdither, sfilt, lfilt,
                                             rpatt, grps, ints, short_pupil, long_pupil,
                                             grismvalue, coordparallel,
                                             i_obs + 1, 1, template_name)
                        APTObservationParams = self.add_exposure(APTObservationParams, direct_tup_to_add)
                        obs_tuple_list.append(tup_to_add)

                    # Now we need to add the two out-of-field exposures, which are
                    # not present in the APT file (but are in the associated pointing
                    # file from APT. We can just
                    # duplicate the entries for the direct images taken immediately
                    # prior. BUT, will there ever be a case where there is no preceding
                    # direct image?
                    APTObservationParams = self.add_exposure(APTObservationParams, direct_tup_to_add)
                    APTObservationParams = self.add_exposure(APTObservationParams, direct_tup_to_add)
                    obs_tuple_list.append(tup_to_add)
                    obs_tuple_list.append(tup_to_add)

            # Now we need to look for mosaic details, if any
            mostile = obs.findall('.//' + apt + 'MosaicTiles')
            n_tiles = len(mostile)

            if n_tiles > 1:
                for i in range(n_tiles - 1):
                    for tup in obs_tuple_list:
                        APTObservationParams = self.add_exposure(APTObservationParams, tup)

            # If WFSC, look at expected groups rather than mosaic tiles:
            if n_tiles == 0 and template_name in ['WfscCommissioning']:
                if num_WFCgroups:
                    n_tiles = num_WFCgroups + 1
            if n_tiles == 0 and template_name in ['WfscGlobalAlignment']:
                n_tiles = n_exp + 1

            # Get observation name, if there is one
            if label == 'None':
                label = ''
            else:
                label = '({})'.format(label)

            print("Found {} mosaic tile(s) for observation {} {}".format(n_tiles, i_obs + 1, label))

        return APTObservationParams


    def add_exposure(self, dictionary, tup):
        # add an exposure to the dictionary
        dictionary['PI_Name'].append(tup[0])
        dictionary['ProposalID'].append(tup[1])
        dictionary['Title'].append(tup[2])
        dictionary['Proposal_category'].append(tup[3])
        dictionary['Science_category'].append(tup[4])
        dictionary['Mode'].append(tup[5])
        dictionary['Module'].append(tup[6])
        dictionary['Subarray'].append(tup[7])
        dictionary['PrimaryDitherType'].append(tup[8])
        dictionary['PrimaryDithers'].append(tup[9])
        dictionary['SubpixelDitherType'].append(tup[10])
        dictionary['SubpixelPositions'].append(tup[11])
        dictionary['ShortFilter'].append(tup[12])
        dictionary['LongFilter'].append(tup[13])
        dictionary['ReadoutPattern'].append(tup[14])
        dictionary['Groups'].append(tup[15])
        dictionary['Integrations'].append(tup[16])
        dictionary['ShortPupil'].append(tup[17])
        dictionary['LongPupil'].append(tup[18])
        dictionary['Grism'].append(tup[19])
        dictionary['CoordinatedParallel'].append(tup[20])
        dictionary['ObservationID'].append(tup[21])
        dictionary['TileNumber'].append(tup[22])
        dictionary['APTTemplate'].append(tup[23])
        return dictionary

    def extract_value(self, line):
        # extract text from xml line
        gt = line.find('>')
        lt = line.find('<', gt)
        return line[gt + 1:lt]


    def expand_for_dithers(self, indict):
        # Expand a given dictionary to create one entry
        # for each dither
        # define the dictionary to hold the expanded entries

        # in here we should also reset the primary and subpixel dither
        # numbers to 1, to avoid confusion.

        expanded = {}
        for key in indict:
            expanded[key] = []

        # loop over entries in dict and duplicate by the
        # number of dither positions
        # keys = np.array(indict.keys())
        keys = indict.keys()
        for i in range(len(indict['PrimaryDithers'])):
            # entry = np.array([item[i] for item in dict.values()])
            arr = np.array([item[i] for item in indict.values()])
            entry = dict(zip(keys, arr))

            # subpix = entry[keys == 'SubpixelPositions']
            subpix = entry['SubpixelPositions']
            if subpix == '0':
                subpix = [[1]]
            if subpix == '4-Point':
                subpix = [[4]]
            if subpix == '9-Point':
                subpix = [[9]]
            # in WFSS, SubpixelPositions will be either '4-Point' or '9-Point'
            # primary = entry[keys == 'PrimaryDithers']
            primary = entry['PrimaryDithers']
            if primary == '0':
                primary = [1]
            reps = np.int(subpix[0][0]) * np.int(primary[0])
            for key in keys:
                for j in range(reps):
                    expanded[key].append(indict[key][i])
        return expanded


    def base36encode(self, integer):
        chars, encoded = '0123456789abcdefghijklmnopqrstuvwxyz', ''

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            encoded = chars[remainder] + encoded

        return encoded.zfill(2)


    def get_pointing_info(self, file, propid):
        # read in information from APT's pointing file
        tar = []
        tile = []
        exp = []
        dith = []
        aperture = []
        targ1 = []
        targ2 = []
        ra = []
        dec = []
        basex = []
        basey = []
        dithx = []
        dithy = []
        v2 = []
        v3 = []
        idlx = []
        idly = []
        level = []
        type_str = []
        expar = []
        dkpar = []
        ddist = []
        observation_number = []
        visit_number = []
        visit_id = []
        visit_grp = []
        activity_id = []
        observation_label = []
        observation_id = []
        seq_id = []

        act_counter = 1
        with open(file) as f:
            for line in f:
                if len(line) > 1:
                    elements = line.split()

                    # look for lines that give visit/observation numbers
                    if line[0:2] == '* ':
                        paren = line.rfind('(')
                        if paren == -1:
                            obslabel = line[2:]
                            obslabel = obslabel.strip()
                        else:
                            obslabel = line[2:paren-1]
                            obslabel = obslabel.strip()
                        if (' (' in obslabel) and (')' in obslabel):
                            obslabel = re.split(r' \(|\)', obslabel)[0]

                        if 'FGS' in obslabel:
                            skip = True
                        else:
                            skip = False

                    if line[0:2] == '**':
                        v = elements[2]
                        obsnum, visitnum = v.split(':')
                        obsnum = str(obsnum).zfill(3)
                        visitnum = str(visitnum).zfill(3)
                        if skip == True:
                            print('Skipping observation {} ({})'.format(obsnum, obslabel))

                    try:
                        # skip the line at the beginning of each
                        # visit that gives info on the target,
                        # but is not actually an observation
                        # These lines have 'Exp' values of 0,
                        # while observations have a value of 1
                        # (that I've seen so far)
                        #
                        # Also, skip non-NIRCam lines. Check for NRC in aperture name
                        if ((np.int(elements[1]) > 0) & ('NRC' in elements[4])):
                            if skip:
                                act_counter += 1
                                continue
                            act = self.base36encode(act_counter)
                            activity_id.append(act)
                            observation_label.append(obslabel)
                            observation_number.append(obsnum)
                            visit_number.append(visitnum)
                            vid = str(propid) + visitnum + obsnum
                            visit_id.append(vid)
                            vgrp = '01'
                            visit_grp.append(vgrp)
                            seq = '1'
                            seq_id.append(seq)
                            tar.append(np.int(elements[0]))
                            tile.append(np.int(elements[1]))
                            exnum = str(elements[2]).zfill(5)
                            exp.append(exnum)
                            dith.append(np.int(elements[3]))

                            ap = elements[4]
                            if ('GRISMR_WFSS' in elements[4]):
                                ap = ap.replace('GRISMR_WFSS', 'FULL')
                            elif ('GRISMC_WFSS' in elements[4]):
                                ap = ap.replace('GRISMC_WFSS', 'FULL')

                            aperture.append(ap)
                            targ1.append(np.int(elements[5]))
                            targ2.append(elements[6])
                            ra.append(elements[7])
                            dec.append(elements[8])
                            basex.append(elements[9])
                            basey.append(elements[10])
                            dithx.append(np.float(elements[11]))
                            dithy.append(np.float(elements[12]))
                            v2.append(np.float(elements[13]))
                            v3.append(np.float(elements[14]))
                            idlx.append(np.float(elements[15]))
                            idly.append(np.float(elements[16]))
                            level.append(elements[17])
                            type_str.append(elements[18])
                            expar.append(np.int(elements[19]))
                            dkpar.append(np.int(elements[20]))
                            ddist.append(np.float(elements[21]))
                            observation_id.append('V' + vid + 'P00000000' + vgrp + seq + act)
                            act_counter  += 1

                    except:
                        pass

        pointing = {'exposure': exp, 'dither': dith, 'aperture': aperture,
                    'targ1': targ1, 'targ2': targ2, 'ra': ra, 'dec': dec,
                    'basex': basex, 'basey': basey, 'dithx': dithx,
                    'dithy': dithy, 'v2': v2, 'v3': v3, 'idlx': idlx,
                    'idly': idly, 'obs_label': observation_label,
                    'obs_num': observation_number, 'visit_num': visit_number,
                    'act_id': activity_id, 'visit_id': visit_id, 'visit_group': visit_grp,
                    'sequence_id': seq_id, 'observation_id': observation_id}
        return pointing


    def combine_dicts(self, dict1, dict2):
        # Now combine the dictionaries from the xml file and the pointing file
        combined = dict1.copy()
        combined.update(dict2)
        return combined


    def expand_for_detectors(self, dict):
        # Expand dictionary to have one line per detector, rather than the
        # one line per module that is in the input
        finaltab = {}
        for key in dict:
            finaltab[key] = []
            # print(key, len(dict[key]))
        finaltab['detector'] = []

        for i in range(len(dict['PrimaryDithers'])):
            # Determine module of the observation
            module = dict['Module'][i]
            if module == 'ALL':
                detectors = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
            elif module == 'A':
                detectors = ['A1', 'A2', 'A3', 'A4', 'A5']
            elif module == 'B':
                detectors = ['B1', 'B2', 'base36encode', 'B4', 'B5']
            elif 'A3' in module:
                detectors = ['A3']
            elif 'B4' in module:
                detectors = ['B4']
            else:
                raise ValueError('Unknown module {}'.format(module))

            for key in dict:
                finaltab[key].extend(([dict[key][i]] * len(detectors)))
            finaltab['detector'].extend(detectors)

        return finaltab


    def ra_dec_update(self):
        # given the v2, v3 values in each entry, calculate RA, Dec

        # read in siaf
        distortionTable = ascii.read(self.siaf, header_start=1)

        aperture_ra = []
        aperture_dec = []
        for i in range(len(self.exposure_tab['Module'])):

            # first find detector
            # need ra, dec and v2, v3 pairs from entry
            # to calculate ra, dec at each detector's reference location
            config = ascii.read('/Users/lchambers/TEL/nircam_simulator/nircam_simulator/config/NIRCam_subarray_definitions.list')
            detector = 'NRC' + self.exposure_tab['detector'][i]
            sub = self.exposure_tab['Subarray'][i]
            aperture = detector + '_' + sub
            if aperture in config['AperName']:
                pass
            else:
                aperture = [apername for apername, name in \
                            np.array(config['AperName', 'Name']) if \
                            (sub in apername) or (sub in name)]
                if len(aperture) > 1 or len(aperture) == 0:
                    raise ValueError('Cannot combine detector {} and subarray {}\
                        into valid aperture name.'.format(detector, sub))
                else:
                    aperture = aperture[0]

            pointing_ra = np.float(self.exposure_tab['ra'][i])
            pointing_dec = np.float(self.exposure_tab['dec'][i])
            pointing_v2 = np.float(self.exposure_tab['v2'][i])
            pointing_v3 = np.float(self.exposure_tab['v3'][i])
            pav3 = np.float(self.exposure_tab['pav3'][i])

            # print('pointing_v2, pointing_v3: ',pointing_v2, pointing_v3)
            # print('pointing ra & dec: ', pointing_ra, pointing_dec)

            # calculate local roll angle
            local_roll = set_telescope_pointing.compute_local_roll(pav3, pointing_ra,
                                                                   pointing_dec,
                                                                   pointing_v2,
                                                                   pointing_v3)
            # print('local_roll: ', local_roll)
            # create attitude_matrix
            attitude_matrix = rotations.attitude(pointing_v2, pointing_v3,
                                                 pointing_ra, pointing_dec, local_roll)
            # print(attitude_matrix)

            # find v2, v3 of the reference location for the detector
            aperture_v2, aperture_v3 = self.ref_location(distortionTable, aperture)

            # print('aperture_v2, aperture_v3: ',aperture_v2[0], aperture_v3[0])

            # print('aperture_v2, aperture_v3: ',aperture_v2[0], aperture_v3[0])

            # calculate RA, Dec of reference location for the detector
            ra, dec = rotations.pointing(attitude_matrix, aperture_v2, aperture_v3)
            # if dec < 0:
            #    print("ra, dec: {}, {}".format(ra, dec))
            #    print('attitude matrix {}'.format(attitude_matrix))
            #    print('aperture v2, v3: {}, {}'.format(aperture_v2.data, aperture_v3.data))
            #    print(pointing_v2, pointing_v3, pointing_ra, pointing_dec, local_roll, pav3)
            #    stop

            # print('ra & dec: ', ra, dec, '\n')

            aperture_ra.append(ra)
            aperture_dec.append(dec)
        self.exposure_tab['ra_ref'] = aperture_ra
        self.exposure_tab['dec_ref'] = aperture_dec

    def ref_location(self, siaf, det):
        # find v2, v3 of detector reference location
        match = siaf['AperName'] == det
        if np.any(match) == False:
            raise ValueError("Aperture name {} not found in input CSV file {}.".
                             format(det, self.siaf))
        v2 = siaf[match]['V2Ref']
        v3 = siaf[match]['V3Ref']
        return v2, v3

    def add_observation_info(self, intab):
        # Add information about each observation. Catalog names,
        # dates, PAV3 values, etc.

        with open(self.observation_table, 'r') as infile:
            self.obstab = yaml.load(infile)

        onames = []
        onums = []
        for key1 in self.obstab:
            onames.append(self.obstab[key1]['Name'])
            onums.append(key1)
        onames = np.array(onames)

        obs_start = []
        obs_pav3 = []
        obs_sw_ptsrc = []
        obs_sw_galcat = []
        obs_sw_ext = []
        obs_sw_extscl = []
        obs_sw_extcent = []
        obs_sw_movptsrc = []
        obs_sw_movgal = []
        obs_sw_movext = []
        obs_sw_movconv = []
        obs_sw_solarsys = []
        obs_sw_bkgd = []
        obs_lw_ptsrc = []
        obs_lw_galcat = []
        obs_lw_ext = []
        obs_lw_extscl = []
        obs_lw_extcent = []
        obs_lw_movptsrc = []
        obs_lw_movgal = []
        obs_lw_movext = []
        obs_lw_movconv = []
        obs_lw_solarsys = []
        obs_lw_bkgd = []

        # This will allow the filter values
        # in the observation table to override
        # what comes from APT. This is useful for
        # WFSS observations, where you want to create
        # 'direct images' in various filters to hand
        # to the disperser software in order to create
        # simulated dispersed data. In that case, you
        # would need to create 'direct images' using
        # several filters inside the filter listed in APT
        # in order to get broadband colors to disperse.
        # If filter names are present in the observation
        # yaml file, then they will be used. If not, the
        # filters from APT will be kept
        obs_sw_filt = []
        obs_lw_filt = []

        for obs in intab['obs_label']:
            match = np.where(obs == onames)[0]
            if len(match) == 0:
                print("No valid epoch line found for observation {}".format(obs))
                print(type(obs))
                print(onames, obs)
                sys.exit()
            else:
                # print('Matching {} from xml with {} from observation listfile'.format(obs, onames[match[0]]))
                # obslist = self.obstab['Observation{}'.format(match[0] + 1)]
                obslist = self.obstab[onums[match[0]]]
                obs_start.append(obslist['Date'].strftime('%Y-%m-%d'))
                obs_pav3.append(obslist['PAV3'])
                obs_sw_ptsrc.append(obslist['SW']['PointSourceCatalog'])
                obs_sw_galcat.append(obslist['SW']['GalaxyCatalog'])
                obs_sw_ext.append(obslist['SW']['ExtendedCatalog'])
                obs_sw_extscl.append(obslist['SW']['ExtendedScale'])
                obs_sw_extcent.append(obslist['SW']['ExtendedCenter'])
                obs_sw_movptsrc.append(obslist['SW']['MovingTargetList'])
                obs_sw_movgal.append(obslist['SW']['MovingTargetSersic'])
                obs_sw_movext.append(obslist['SW']['MovingTargetExtended'])
                obs_sw_movconv.append(obslist['SW']['MovingTargetConvolveExtended'])
                obs_sw_solarsys.append(obslist['SW']['MovingTargetToTrack'])
                obs_sw_bkgd.append(obslist['SW']['BackgroundRate'])
                obs_lw_ptsrc.append(obslist['LW']['PointSourceCatalog'])
                obs_lw_galcat.append(obslist['LW']['GalaxyCatalog'])
                obs_lw_ext.append(obslist['LW']['ExtendedCatalog'])
                obs_lw_extscl.append(obslist['LW']['ExtendedScale'])
                obs_lw_extcent.append(obslist['LW']['ExtendedCenter'])
                obs_lw_movptsrc.append(obslist['LW']['MovingTargetList'])
                obs_lw_movgal.append(obslist['LW']['MovingTargetSersic'])
                obs_lw_movext.append(obslist['LW']['MovingTargetExtended'])
                obs_lw_movconv.append(obslist['LW']['MovingTargetConvolveExtended'])
                obs_lw_solarsys.append(obslist['LW']['MovingTargetToTrack'])
                obs_lw_bkgd.append(obslist['LW']['BackgroundRate'])

                # Override filters if given
                try:
                    obs_sw_filt.append(obslist['SW']['Filter'])
                except:
                    pass
                try:
                    obs_lw_filt.append(obslist['LW']['Filter'])
                except:
                    pass

        intab['epoch_start_date'] = obs_start
        intab['pav3'] = obs_pav3
        intab['sw_ptsrc'] = obs_sw_ptsrc
        intab['sw_galcat'] = obs_sw_galcat
        intab['sw_ext'] = obs_sw_ext
        intab['sw_extscl'] = obs_sw_extscl
        intab['sw_extcent'] = obs_sw_extcent
        intab['sw_movptsrc'] = obs_sw_movptsrc
        intab['sw_movgal'] = obs_sw_movgal
        intab['sw_movext'] = obs_sw_movext
        intab['sw_movconv'] = obs_sw_movconv
        intab['sw_solarsys'] = obs_sw_solarsys
        intab['sw_bkgd'] = obs_sw_bkgd
        intab['lw_ptsrc'] = obs_lw_ptsrc
        intab['lw_galcat'] = obs_lw_galcat
        intab['lw_ext'] = obs_lw_ext
        intab['lw_extscl'] = obs_lw_extscl
        intab['lw_extcent'] = obs_lw_extcent
        intab['lw_movptsrc'] = obs_lw_movptsrc
        intab['lw_movgal'] = obs_lw_movgal
        intab['lw_movext'] = obs_lw_movext
        intab['lw_movconv'] = obs_lw_movconv
        intab['lw_solarsys'] = obs_lw_solarsys
        intab['lw_bkgd'] = obs_lw_bkgd

        # Here we override the filters read from APT
        # if they are given in the observation yaml file
        if len(obs_sw_filt) > 0:
            intab['ShortFilter'] = obs_sw_filt
        if len(obs_lw_filt) > 0:
            intab['LongFilter'] = obs_lw_filt
        return intab


    def add_epochs(self, intab):
        # add information on the epoch of each observation
        # if the user entered a list of epochs, read that in
        default_date = '2020-10-14'

        if self.epoch_list is not None:
            epochs = ascii.read(self.epoch_list, header_start=0, data_start=1)
        else:
            epochs = Table()
            epochs['observation'] = intab['obs_label']
            epochs['date'] = ['2018-10-14'] * len(intab['obs_label'])
            epochs['pav3'] = [0.] * len(intab['obs_label'])

        # insert epoch info for each observation into dictionary
        epoch_start = []
        epoch_pav3 = []
        for obs in intab['obs_label']:
            match = obs == epochs['observation'].data
            if np.sum(match) == 0:
                print("No valid epoch line found for observation {}".format(obs))
                print(epochs['observation'].data)
                epoch_start.append(default_date)
                epoch_pav3.append(0.)
            else:
                epoch_start.append(epochs['date'][match].data[0])
                epoch_pav3.append(epochs['pav3'][match].data[0])
        intab['epoch_start_date'] = epoch_start
        intab['pav3'] = epoch_pav3
        return intab

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("input_xml", help='XML file from APT describing the observations.')
        parser.add_argument("pointing_file", help='Pointing file from APT describing observations.')
        parser.add_argument("siaf", help='Name of CSV version of SIAF')
        parser.add_argument("--output_csv", help="Name of output CSV file containing list of observations.", default=None)
        parser.add_argument("--observation_table", help='Ascii file containing a list of observations, times, and roll angles, catalogs', default=None)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: apt_inputs.py NIRCam_obs.xml NIRCam_obs.pointing SIAF_March2017.csv'

    input = AptInput()
    parser = input.add_options(usage = usagestring)
    args = parser.parse_args(namespace=input)
    input.create_input_table()
