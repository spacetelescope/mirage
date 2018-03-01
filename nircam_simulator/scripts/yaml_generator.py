# ! /usr/bin/env python

'''
Function to produce yaml files that can be used as input for
the ramp simulator

Inputs:

xml file - Name of xml file exported from APT.
pointing file - Name of associated pointing file exported from APT.
siaf - Name of csv version of SIAF.

Optional inputs:

output_dir - Directory into which the output yaml files are written

simdata_output_dir - Directory to place in the output_directory field of the yaml files.
                     This is the directory where the simulated ramps will be saved.

table_file - Ascii table containing observation info. This is the
             output from apt_inputs.py Use this if you are
             not providing xml and pointing files from APT.

datatype - Specifies the type of output data to save. Can be "raw", in which
           case the raw (uncalibrated) file is saved, "linear", where the
           linearized file is saved and ready to be run through the jump
           detection and ramp-fitting steps of the pipeline, or
           "linear, raw", where both versions are saved.

use_nonstsci_names - set to True to override the use of the standard
                     STScI naming convention for output files

subarray_def_file - Ascii file containing NIRCam subarray definitions

readpatt_def_file - Ascii file containing NIRCam readout pattern definitions

point_source - point source catalog file. Can be a single file, or a list of
               catalogs. If it is a list, each filename is expected to contain
               the filter name for which it is to be used. Catalogs and filters
               will then be matched up in the output yaml files.

galaxyListFile - galaxy (sersic) source catalog file. Can be a single name, or
                 a list of names. Behavior is identical to point_source above.

extended - extended source catalog file. Behavior is identical to point_source
            above.

convolveExtended - Set to True to convolve extended sources with NIRCam PSF

movingTarg - Moving (point source) target catalog (sources moving through fov) names.
             Behavior is the same as point_sources above.

movingTargSersic - Moving galaxy (sersic) target catalog (sources moving through fov)
                   names. Behavior is the same as point_sources above.

movingTargExtended - Moving extended source target catalog (sources moving through fov)
                     names. Behavior is the same as point_sources above.

movingTargToTrack - Catalog of non-sidereal targets for non-sidereal tracking observations.
                    Behavior is the same as point_sources above.

bkgdrate - Uniform background rate (e-/s) to add to observation.

epoch_list - Ascii table file containing epoch start times and telescope roll angles
             to use for each observation.


Dependencies:
argparse, astropy, numpy, glob, copy

apt_inputs.py - Functions for reading and parsing xml and pointing files from APT.

HISTORY:

July 2017 - V0: Initial version. Bryan Hilbert
Feb 2018 - V1: Updates to accomodate multiple filter pairs per 
               observation. Launren Chambers

'''

import sys
import os
import pkg_resources
from glob import glob
from copy import deepcopy
import argparse
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.io import ascii
from . import apt_inputs


class SimInput:
    def __init__(self):
        # Set the NIRCAM_SIM_DATA environment variable
        # if it's not already
        stsci_local = '/ifs/jwst/wit/nircam/nircam_simulator_data/'
        self.datadir = os.environ.get('NIRCAM_SIM_DATA')
        if self.datadir is None:
            if os.path.exists(stsci_local):
                self.datadir = stsci_local
            else:
                print("WARNING: NIRCAM_SIM_DATA environment variable")
                print("has not been set, and it appears you do not")
                print("have access to the STScI disks. You must set")
                print("the environment variable to point at the directory")
                print("containing the data files needed by the simulator.")
                sys.exit()

        self.info = {}
        self.input_xml = None
        self.pointing_file = None
        self.siaf = None
        self.reffile_setup()
        self.datatype = 'linear'
        self.output_dir = './'
        self.table_file = None
        self.use_nonstsci_names = False
        self.subarray_def_file = 'config'
        self.readpatt_def_file = 'config'
        self.crosstalk = 'config'
        self.filtpupil_pairs = 'config'
        self.fluxcal = 'config'
        self.dq_init_config = 'config'
        self.saturation_config = 'config'
        self.superbias_config = 'config'
        self.refpix_config = 'config'
        self.linearity_config = 'config'
        self.filter_throughput = 'config'
        self.observation_table = None
        self.use_JWST_pipeline = True
        self.use_linearized_darks = False
        self.simdata_output_dir = './'
        self.psfpath = '/ifs/jwst/wit/nircam/nircam_simulator_data/webbpsf_library'
        self.psfbasename = 'nircam'
        self.psfpixfrac = 0.25
        self.psfwfe = 'predicted'
        self.psfwfegroup = 0
        self.resets_bet_ints = 1 # NIRCam should be 1
        self.tracking = 'sidereal'

        # Prepare to find files listed as 'config'
        self.modpath = pkg_resources.resource_filename('nircam_simulator', '')
        self.configfiles = {}
        self.configfiles['subarray_def_file'] = 'NIRCam_subarray_definitions.list'
        self.configfiles['fluxcal'] = 'NIRCam_zeropoints.list'
        self.configfiles['filtpupil_pairs'] = 'nircam_filter_pupil_pairings.list'
        self.configfiles['readpatt_def_file'] = 'nircam_read_pattern_definitions.list'
        self.configfiles['crosstalk'] = 'xtalk20150303g0.errorcut.txt'
        self.configfiles['dq_init_config'] = 'dq_init.cfg'
        self.configfiles['saturation_config'] = 'saturation.cfg'
        self.configfiles['superbias_config'] = 'superbias.cfg'
        self.configfiles['refpix_config'] = 'refpix.cfg'
        self.configfiles['linearity_config'] = 'linearity.cfg'
        self.configfiles['filter_throughput'] = 'placeholder.txt'

    def create_inputs(self):
        # Use full paths for inputs
        self.path_defs()

        if ((self.input_xml is not None) &
            (self.pointing_file is not None) &
            (self.siaf is not None) &
            (self.observation_table is not None)):

            print('Using {}, \n      {}, \n      {}, and \n      {} \n      to generate observation table.\n'.
                  format(self.observation_table, self.input_xml, self.pointing_file, self.siaf))

            # Define directories and paths
            indir, infile = os.path.split(self.input_xml)
            final_file = os.path.join(self.output_dir,
                                      'Observation_table_for_' + infile + \
                                      '_with_yaml_parameters.csv')

            # Read XML file and make observation table
            apt = apt_inputs.AptInput()
            apt.input_xml = self.input_xml
            apt.pointing_file = self.pointing_file
            apt.siaf = self.siaf
            apt.observation_table = self.observation_table
            apt.create_input_table()
            self.info = apt.exposure_tab
            # pprint.pprint(self.info)

            # Add start time info to each element
            self.make_start_times()

            # Add a list of output yaml names to the dictionary
            self.make_output_names()

            # Add source catalogs
            # self.add_catalogs()

        elif self.table_file is not None:
            print('Reading table file: {}'.format(self.table_file))
            info = ascii.read(self.table_file)
            self.info = self.table_to_dict(info)
            final_file = self.table_file + '_with_yaml_parameters.csv'

        else:
            print("WARNING. You must include either an ascii table file of observations.")
            print("or xml and pointing files from APT plus an ascii siaf table.")
            sys.exit()

        # For each element in the lists below, use the detector name to
        # find the appropriate reference files. Create lists, and add
        # these lists to the dictionary
        darks = []
        lindarks = []
        superbias = []
        linearity = []
        saturation = []
        gain = []
        astrometric = []
        ipc = []
        pam = []

        for det in self.info['detector']:
            darks.append(self.get_dark(det))
            lindarks.append(self.get_lindark(det))
            superbias.append(self.get_reffile(self.superbias_list, det))
            linearity.append(self.get_reffile(self.linearity_list, det))
            saturation.append(self.get_reffile(self.saturation_list, det))
            gain.append(self.get_reffile(self.gain_list, det))
            astrometric.append(self.get_reffile(self.astrometric_list, det))
            ipc.append(self.get_reffile(self.ipc_list, det))
            pam.append(self.get_reffile(self.pam_list, det))
        self.info['dark'] = darks
        # If linearized darks are to be used, set the darks to None
        if self.use_linearized_darks:
            self.info['dark'] = [None] * len(darks)
            self.info['lindark'] = lindarks
        else:
            self.info['dark'] = darks
            self.info['lindark'] = [None] * len(lindarks)
        self.info['superbias'] = superbias
        self.info['linearity'] = linearity
        self.info['saturation'] = saturation
        self.info['gain'] = gain
        self.info['astrometric'] = astrometric
        self.info['ipc'] = ipc
        self.info['pixelAreaMap'] = pam

        # Add setting describing whether JWST pipeline will be used
        self.info['use_JWST_pipeline'] = [self.use_JWST_pipeline] * len(darks)

        # add background rate to the table
        # self.info['bkgdrate'] = np.array([self.bkgdrate]*len(self.info['Mode']))

        # grism entries
        grism_source_image = ['False'] * len(self.info['Mode'])
        grism_input_only = ['False'] * len(self.info['Mode'])
        for i in range(len(self.info['Mode'])):
            if self.info['Mode'][i] == 'WFSS':
                grism_source_image[i] = 'True'
                grism_input_only[i] = 'True'
                # SW detectors shouldn't be wfss
                if self.info['detector'][i][-1] != '5':
                    self.info['Mode'][i] = 'imaging'
        self.info['grism_source_image'] = grism_source_image
        self.info['grism_input_only'] = grism_input_only

        # level-3 associated keywords that are not present in APT file.
        # not quite sure how to populate these
        self.info['visit_group'] = ['01'] * len(self.info['Mode'])
        # self.info['sequence_id'] = ['1'] * len(self.info['Mode'])
        seq = []
        for par in self.info['CoordinatedParallel']:
            if par.lower() == 'true':
                seq.append('2')
            if par.lower() == 'false':
                seq.append('1')
        self.info['sequence_id'] = seq
        self.info['obs_template'] = ['NIRCam Imaging'] * len(self.info['Mode'])

        # write out the updated table, including yaml filenames
        # start times, and reference files
        print('Updated observation table file saved to {}'.format(final_file))
        ascii.write(Table(self.info), final_file, format='csv', overwrite=True)

        # Now go through the lists one element at a time
        # and create a yaml file for each.
        yamls = []
        for i in range(len(self.info['detector'])):
            file_dict = {}
            for key in self.info:
                file_dict[key] = self.info[key][i]

            # break dither number into numbers for primary
            # and subpixel dithers
            tot_dith = np.int(file_dict['dither'])
            primarytot = np.int(file_dict['PrimaryDithers'])
            try:
                subpixtot = np.int(file_dict['SubpixelPositions'])
            except:
                subpixtot = np.int(file_dict['SubpixelPositions'][0])
            primary_dither = np.ceil(1. * tot_dith / subpixtot)
            subpix_dither = tot_dith - (primary_dither * primarytot * subpixtot - subpixtot)
            file_dict['primary_dither_num'] = primary_dither
            file_dict['subpix_dither_num'] = subpix_dither
            file_dict['siaf'] = self.siaf
            fname = self.write_yaml(file_dict)
            yamls.append(fname)
            
        # Write out summary of all written yaml files
        yaml_path = os.path.join(self.output_dir, 'V*.yaml')
        #yamls = glob(yaml_path)

        filenames = [y.split('/')[-1] for y in yamls]
        mosaic_numbers = sorted(list(set([f.split('_')[0] for f in filenames])))
        obs_ids = sorted(list(set([m[9:12] for m in mosaic_numbers])))

        print('\n')
        i_mod = 0
        for obs in obs_ids:
            n_visits = len(list(set([m[6:9] for m in mosaic_numbers if m[9:12] == obs])))
            n_tiles = len(list(set([m[-2:] for m in mosaic_numbers if m[9:12] == obs])))
            module = self.info['Module'][i_mod]

            if module in ['A', 'B']:
                n_det = 5
                module = ' ' + module
            if module == 'ALL':
                n_det = 10
                module = 's A and B'
            if 'A3' in module:
                n_det = 1
                module = ' A3'
            if 'B4' in module:
                n_det = 1
                module = ' B4'

            i_mod += n_tiles * n_det

            print('Observation {}: \n   {} visit(s) \n   {} exposure(s)\n   {} detector(s) in module{}'.format(obs, n_visits, n_tiles, n_det, module))
        print('\n{} exposures total.'.format(len(mosaic_numbers)))
        print('{} output files written to: {}'.format(len(yamls), self.output_dir))

    def path_defs(self):
        """Expand input files to have full paths"""
        self.input_xml = os.path.abspath(self.input_xml)
        self.pointing_file = os.path.abspath(self.pointing_file)
        self.siaf = os.path.abspath(self.siaf)
        self.output_dir = os.path.abspath(self.output_dir)
        self.simdata_output_dir = os.path.abspath(self.simdata_output_dir)
        if self.table_file is not None:
            self.table_file = os.path.abspath(self.table_file)
        self.subarray_def_file = self.set_config(self.subarray_def_file, 'subarray_def_file')
        self.readpatt_def_file = self.set_config(self.readpatt_def_file, 'readpatt_def_file')
        self.filtpupil_pairs = self.set_config(self.filtpupil_pairs, 'filtpupil_pairs')
        self.fluxcal = self.set_config(self.fluxcal, 'fluxcal')
        self.filter_throughput = self.set_config(self.filter_throughput, 'filter_throughput')
        self.dq_init_config = self.set_config(self.dq_init_config, 'dq_init_config')
        self.refpix_config = self.set_config(self.refpix_config, 'refpix_config')
        self.saturation_config = self.set_config(self.saturation_config, 'saturation_config')
        self.superbias_config = self.set_config(self.superbias_config, 'superbias_config')
        self.linearity_config = self.set_config(self.linearity_config, 'linearity_config')
        if self.observation_table is not None:
            self.observation_table = os.path.abspath(self.observation_table)
        if self.crosstalk not in [None, 'config']:
            self.crosstalk = os.path.abspath(self.crosstalk)
        elif self.crosstalk == 'config':
            self.crosstalk = os.path.join(self.modpath, 'config', self.configfiles['crosstalk'])

    def set_config(self, file, prop):
        """
        If a given file is listed as 'config'
        then set it in the yaml output as being in
        the config subdirectory.

        Parameters:
        -----------
        file -- Name of the input file
        prop -- String. Type of file that file is.

        Returns:
        --------
        Full path name to the input file
        """
        if file.lower() not in ['config']:
            file = os.path.abspath(file)
        elif file.lower() == 'config':
            file = os.path.join(self.modpath, 'config', self.configfiles[prop])
        return file

    def add_catalogs(self):
        """
        Add list(s) of source catalogs to the table containing the
        observation information
        """
        self.info['point_source'] = [None] * len(self.info['Module'])
        self.info['galaxyListFile'] = [None] * len(self.info['Module'])
        self.info['extended'] = [None] * len(self.info['Module'])
        self.info['convolveExtended'] = [False] * len(self.info['Module'])
        self.info['movingTarg'] = [None] * len(self.info['Module'])
        self.info['movingTargSersic'] = [None] * len(self.info['Module'])
        self.info['movingTargExtended'] = [None] * len(self.info['Module'])
        self.info['movingTargToTrack'] = [None] * len(self.info['Module'])

        for i in range(len(self.info['ShortFilter'])):
            if np.int(self.info['detector'][i][-1]) < 5:
                filtkey = 'ShortFilter'
                pupilkey = 'ShortPupil'
            else:
                filtkey = 'LongFilter'
                pupilkey = 'LongPupil'
            filt = self.info[filtkey][i]
            pup = self.info[pupilkey][i]

            if self.point_source[0] is not None:
                # In here, we assume the user provided a catalog to go with each filter
                # so now we need to find the filter for each entry and generate a list that makes sense
                self.info['point_source'][i] = self.catalog_match(filt, pup, self.point_source, 'point source')
            if self.galaxyListFile[0] is not None:
                self.info['galaxyListFile'][i] = self.catalog_match(filt, pup, self.galaxyListFile, 'galaxy')
            if self.extended[0] is not None:
                self.info['extended'][i] = self.catalog_match(filt, pup, self.extended, 'extended')
            if self.movingTarg[0] is not None:
                self.info['movingTarg'][i] = self.catalog_match(filt, pup, self.movingTarg, 'moving point source target')
            if self.movingTargSersic[0] is not None:
                self.info['movingTargSersic'][i] = self.catalog_match(filt, pup, self.movingTargSersic, 'moving sersic target')
            if self.movingTargExtended[0] is not None:
                self.info['movingTargExtended'][i] = self.catalog_match(filt, pup, self.movingTargExtended, 'moving extended target')
            if self.movingTargToTrack[0] is not None:
                self.info['movingTargToTrack'][i] = self.catalog_match(filt, pup, self.movingTargToTrack, 'non-sidereal moving target')

        if self.convolveExtended == True:
            self.info['convolveExtended'] = [True] * len(self.info['Module'])


    def catalog_match(self, filter, pupil, catalog_list, cattype):
        """
        Given a filter and pupil value, along with a list of input
        catalogs, find the catalog names that contain each filter/
        pupil name.

        Parameters:
        -----------
        filter -- String name of a filter element
        pupil -- String name of a pupil element
        catalog_list -- List of catalog filenames
        cattype -- String containing the type of catalog in
                   the list.

        Returns:
        --------
        Name of catalog that contains the name of the
        input filter/pupil element
        """
        if pupil[0].upper() == 'F':
            match = [s for s in catalog_list if pupil.lower() in s.lower()]
            if len(match) == 0:
                self.no_catalog_match(pupil, cattype)
                return None
            elif len(match) > 1:
                self.multiple_catalog_match(pupil, cattype, match)
            return match[0]
        else:
            match = [s for s in catalog_list if filter.lower() in s.lower()]
            if len(match) == 0:
                self.no_catalog_match(filter, cattype)
                return None
            elif len(match) > 1:
                self.multiple_catalog_match(filter, cattype, match)
            return match[0]

    def no_catalog_match(self, filter, cattype):
        """
        Tell user if no catalog match was found.

        Parameters:
        -----------
        filter -- String name of filter element
        cattype -- String, type of catalog (e.g. pointsource)
        """
        print("WARNING: unable to find filter ({}) name".format(filter))
        print("in any of the given {} inputs".format(cattype))
        print("Using the first input for now. Make sure input catalog names have")
        print("the appropriate filter name in the filename to get matching to work.")

    def multiple_catalog_match(self, filter, cattype, matchlist):
        """
        Tell the user if more than one catalog matches the filter/pupil

        Parameters:
        filter -- String name of filter element
        cattype -- String, type of catalog (e.g. pointsource)
        matchlist -- List of matching catalog names
        """
        print("WARNING: multiple {} catalogs matched! Using the first.".format(cattype))
        print("Observation filter: {}".format(filter))
        print("Matched point source catalogs: {}".format(matchlist))

    def table_to_dict(self, tab):
        """
        Convert the ascii table of observations to a dictionary

        Parameters:
        -----------
        tab -- astropy Table containing observation information

        Returns:
        --------
        Dictionary of observation information
        """
        dict = {}
        for colname in tab.colnames:
            dict[colname] = tab[colname].data
        return dict

    def make_start_times(self):
        """
        Create exposure start times for each entry in
        the observation dictionary
        """
        date_obs = []
        time_obs = []
        expstart = []
        nframe = []
        nskip = []
        namp = []

        # choose arbitrary start time for each epoch
        epoch_base_time = '16:44:12'
        epoch_base_time0 = deepcopy(epoch_base_time)

        epoch_base_date = self.info['epoch_start_date'][0]
        base = Time(epoch_base_date + 'T' + epoch_base_time)
        base_date, base_time = base.iso.split()

        # Pick some arbirary overhead values
        act_overhead = 90  # seconds. (filter change)
        visit_overhead = 600  # seconds. (slew)

        # Get visit, activity_id info for first exposure
        actid = self.info['act_id'][0]
        visit = self.info['visit_num'][0]
        obsname = self.info['obs_label'][0]

        # Read in file containing subarray definitions
        subarray_def = self.get_subarray_defs()

        # Now read in readpattern definitions
        readpatt_def = self.get_readpattern_defs()

        for i in range(len(self.info['Module'])):
            # Get dither/visit
            # Files with the same activity_id should have the same start time
            # Overhead after a visit break should be large, smaller between
            # exposures within a visit
            next_actid = self.info['act_id'][i]
            next_visit = self.info['visit_num'][i]
            next_obsname = self.info['obs_label'][i]

            # Get the values of nframes, nskip, and namp
            readpatt = self.info['ReadoutPattern'][i]

            # Find the readpattern of the file
            readpatt = self.info['ReadoutPattern'][i]
            groups = np.int(self.info['Groups'][i])
            integrations = np.int(self.info['Integrations'][i])

            match2 = readpatt == readpatt_def['name']
            if np.sum(match2) == 0:
                print("WARNING!! Readout pattern {} not found in definition file.".format(readpatt))
                sys.exit()

            # Now get nframe and nskip so we know how many frames in a group
            fpg = np.int(readpatt_def['nframe'][match2][0])
            spg = np.int(readpatt_def['nskip'][match2][0])
            nframe.append(fpg)
            nskip.append(spg)

            # need to find number of amps used
            sub = self.info['Subarray'][i]
            det = 'NRC' + self.info['detector'][i]
            aperture = det + '_' + sub

            match = aperture == subarray_def['AperName']

            if np.sum(match) == 0:
                config = ascii.read(self.subarray_def_file)
                aperture = [apername for apername, name in \
                            np.array(config['AperName', 'Name']) if \
                            (sub in apername) or (sub in name)]

                match = aperture == subarray_def['AperName']

                if len(aperture) > 1 or len(aperture) == 0 or np.sum(match) == 0:
                    raise ValueError('Cannot combine detector {} and subarray {}\
                        into valid aperture name.'.format(det, sub))

            amp = subarray_def['num_amps'][match][0]
            namp.append(amp)

            # same activity ID
            if next_actid == actid:
                # in this case, the start time should remain the same
                date_obs.append(base_date)
                time_obs.append(base_time)
                expstart.append(base.mjd)
                # print(actid, visit, obsname, base_date, base_time)
                continue

            epoch_date = self.info['epoch_start_date'][i]
            epoch_time = deepcopy(epoch_base_time0)
            # new epoch - update the base time
            if epoch_date != epoch_base_date:
                epoch_base_date = deepcopy(epoch_date)
                base = Time(epoch_base_date + 'T' + epoch_base_time)
                base_date, base_time = base.iso.split()
                basereset = True
                date_obs.append(base_date)
                time_obs.append(base_time)
                expstart.append(base.mjd)
                actid = deepcopy(next_actid)
                visit = deepcopy(next_visit)
                obsname = deepcopy(next_obsname)
                continue

            # new visit
            if next_visit != visit:
                # visit break. Larger overhead
                overhead = visit_overhead
            elif ((next_actid > actid) & (next_visit == visit)):
                # same visit, new activity. Smaller overhead
                overhead = act_overhead
            else:
                # should never get in here
                print("Error. Fix me")
                sys.exit()

            # For cases where the base time needs to change
            # continue down here
            xs = subarray_def['xstart'][match][0]
            xe = subarray_def['xend'][match][0]
            ys = subarray_def['ystart'][match][0]
            ye = subarray_def['yend'][match][0]
            xd = xe - xs + 1
            yd = ye - ys + 1
            frametime = self.calcFrameTime(xd, yd, amp)

            # Estimate total exposure time
            exptime = ((fpg + spg) * groups + fpg) * integrations * frametime

            # Delta should include the exposure time, plus overhead
            delta = TimeDelta(exptime + overhead, format='sec')
            base += delta
            base_date, base_time = base.iso.split()

            # Add updated dates and times to the list
            date_obs.append(base_date)
            time_obs.append(base_time)
            expstart.append(base.mjd)

            # increment the activity ID and visit
            actid = deepcopy(next_actid)
            visit = deepcopy(next_visit)
            obsname = deepcopy(next_obsname)

        self.info['date_obs'] = date_obs
        self.info['time_obs'] = time_obs
        # self.info['expstart'] = expstart
        self.info['nframe'] = nframe
        self.info['nskip'] = nskip
        self.info['namp'] = namp

    def get_readpattern_defs(self):
        """Read in the readpattern definition file and return table"""
        tab = ascii.read(self.readpatt_def_file)
        return tab

    def get_subarray_defs(self):
        """Read in subarray definition file and return table"""
        sub = ascii.read(self.subarray_def_file)
        return sub

    def calcFrameTime(self, xd, yd, namp):
        """
        Calculate the exposure time of a single frame of the proposed output ramp
        based on the size of the croped dark current integration

        Parameters:
        -----------
        xd -- Integer x-coordinate size of the aperture to read out, in pixels
        yd -- Integer y-coordinate size of the aperture to read out, in pixels
        nanp -- Integer number of amplifiers used to read out the aperture

        Returns:
        --------
        The amount of time needed, in seconds, to read out the detector
        """
        return (xd / namp + 12.) * (yd + 1) * 10.00 * 1.e-6

    def make_output_names(self):
        """
        Create output yaml file names to go with all of the
        entries in the dictionary
        """
        onames = []
        fnames = []
        for i in range(len(self.info['Module'])):
            act = str(self.info['act_id'][i]).zfill(2)
            det = self.info['detector'][i]
            mode = self.info['Mode'][i]
            dither = str(self.info['dither'][i]).zfill(2)
            onames.append(os.path.abspath(os.path.join(self.output_dir, 'Act{}_{}_{}_Dither{}.yaml'.format(act, det, mode, dither))))
            fnames.append('Act{}_{}_{}_Dither{}_uncal.fits'.format(act, det, mode, dither))
        self.info['yamlfile'] = onames
        self.info['outputfits'] = fnames

    def get_dark(self, detector):
        """
        Return the name of a dark current file to use as input
        based on the detector being used

        Parameters:
        -----------
        detector -- string name of detector being used

        Returns:
        --------
        Name of a dark current file to use for this detector
        """
        files = self.dark_list[detector]
        rand_index = np.random.randint(0, len(files) - 1)
        return files[rand_index]

    def get_lindark(self, detector):
        """
        Return the name of a linearized dark current file to 
        use as input based on the detector being used

        Parameters:
        -----------
        detector -- string name of detector being used

        Returns:
        --------
        Name of a linearized dark current file to use for this detector
        """
        files = self.lindark_list[detector]
        rand_index = np.random.randint(0, len(files) - 1)
        return files[rand_index]

    def get_reffile(self, refs, detector):
        """
        Return the appropriate reference file for detector
        and given reference file dictionary. 

        Parameters:
        -----------
        refs -- dictionary in the form of:
             {'A1':'filenamea1.fits', 'A2':'filenamea2.fits'...}
        detector -- String name of detector

        Returns:
        --------
        Name of reference file appropriate for given detector
        """
        for key in refs:
            if detector in key:
                return refs[key]
        print("WARNING: no file found for detector {} in {}"
              .format(detector, refs))

    def write_yaml(self, input):
        """
        Create yaml file for a single exposure/detector

        Parameters:
        -----------
        input -- dictionary containing all needed exposure 
                 information for one exposure
        """

        # select the right filter
        if np.int(input['detector'][-1]) < 5:
            filtkey = 'ShortFilter'
            pupilkey = 'ShortPupil'
            catkey = 'sw'
        else:
            filtkey = 'LongFilter'
            pupilkey = 'LongPupil'
            catkey = 'lw'

        if self.use_nonstsci_names:
            outtf = False
            outfile = input['outputfits']
            yamlout = input['yamlfile']
        else:
            outtf = True
            outfile = input['observation_id'] + '_' + input['detector'] + '_' + input[filtkey] + '_uncal.fits'
            yamlout = input['observation_id'] + '_' + input['detector'] + '_' + input[filtkey] + '.yaml'

        yamlout = os.path.join(self.output_dir, yamlout)
        with open(yamlout, 'w') as f:
            f.write('Inst:\n')
            f.write('  instrument: {}          # Instrument name\n'.format('NIRCam'))
            f.write('  mode: {}                # Observation mode (e.g. imaging, WFSS)\n'.format(input['Mode']))
            f.write('  use_JWST_pipeline: {}   # Use pipeline in data transformations\n'.format(input['use_JWST_pipeline']))
            f.write('\n')
            f.write('Readout:\n')
            f.write('  readpatt: {}        # Readout pattern (RAPID, BRIGHT2, etc) overrides nframe, nskip unless it is not recognized\n'.format(input['ReadoutPattern']))
            f.write('  ngroup: {}              # Number of groups in integration\n'.format(input['Groups']))
            f.write('  nint: {}          # Number of integrations per exposure\n'.format(input['Integrations']))
            f.write('  resets_bet_ints: {} #Number of detector resets between integrations\n'.format(self.resets_bet_ints))

            apunder = input['aperture'].find('_')
            full_ap = 'NRC' + input['detector'] + '_' + input['aperture'][apunder + 1:]

            scripts_path = os.path.dirname(os.path.realpath(__file__))
            modpath = os.path.split(scripts_path)[0]
            subarray_def_file = os.path.join(modpath, 'config', 'NIRCam_subarray_definitions.list')
            config = ascii.read(subarray_def_file)

            if full_ap not in config['AperName']:
                full_ap = [apername for apername, name in \
                           np.array(config['AperName', 'Name']) if \
                           (full_ap in apername) or (full_ap in name)]
                if len(full_ap) > 1 or len(full_ap) == 0:
                    raise ValueError('Cannot match {} with valid aperture name.'
                                     .format(full_ap))
                else:
                    full_ap = full_ap[0]

            f.write('  array_name: {}    # Name of array (FULL, SUB160, SUB64P, etc) overrides subarray_bounds below\n'.format(full_ap))
            f.write('  filter: {}       # Filter of simulated data (F090W, F322W2, etc)\n'.format(input[filtkey]))
            f.write('  pupil: {}        # Pupil element for simulated data (CLEAR, GRISMC, etc)\n'.format(input[pupilkey]))
            f.write('\n')
            f.write('Reffiles:                                 # Set to None or leave blank if you wish to skip that step\n')
            f.write('  dark: {}   # Dark current integration used as the base\n'.format(input['dark']))
            f.write('  linearized_darkfile: {}   # Linearized dark ramp to use as input. Supercedes dark above\n'.format(input['lindark']))
            f.write('  badpixmask: None   # If linearized dark is used, populate output DQ extensions using this file\n')
            f.write('  superbias: {}     # Superbias file. Set to None or leave blank if not using\n'.format(input['superbias']))
            f.write('  linearity: {}    # linearity correction coefficients\n'.format(input['linearity']))
            f.write('  saturation: {}    # well depth reference files\n'.format(input['saturation']))
            f.write('  gain: {} # Gain map\n'.format(input['gain']))
            f.write('  pixelflat: None \n')
            f.write('  illumflat: None                               # Illumination flat field file\n')
            f.write('  astrometric: {}  # Astrometric distortion file (asdf)\n'.format(input['astrometric']))
            f.write('  distortion_coeffs: {}        # CSV file containing distortion coefficients\n'.format(input['siaf']))
            f.write('  ipc: {} # File containing IPC kernel to apply\n'.format(input['ipc']))
            f.write('  invertIPC: True       # Invert the IPC kernel before the convolution. True or False. Use True if the kernel is designed for the removal of IPC effects, like the JWST reference files are.\n')
            f.write('  occult: None                                    # Occulting spots correction image\n')
            f.write('  pixelAreaMap: {}      # Pixel area map for the detector. Used to introduce distortion into the output ramp.\n'.format(input['pixelAreaMap']))
            f.write('  subarray_defs: {} # File that contains a list of all possible subarray names and coordinates\n'.format(self.subarray_def_file))
            f.write('  readpattdefs: {}  # File that contains a list of all possible readout pattern names and associated NFRAME/NSKIP values\n'.format(self.readpatt_def_file))
            f.write('  crosstalk: {}   # File containing crosstalk coefficients\n'.format(self.crosstalk))
            f.write('  filtpupilcombo: {}   # File that lists the filter wheel element / pupil wheel element combinations. Used only in writing output file\n'.format(self.filtpupil_pairs))
            f.write('  flux_cal: {} # File that lists flux conversion factor and pivot wavelength for each filter. Only used when making direct image outputs to be fed into the grism disperser code.\n'.format(self.fluxcal))
            f.write('  filter_throughput: {} #File containing filter throughput curve\n'.format(self.filter_throughput))
            f.write('\n')
            f.write('nonlin:\n')
            f.write('  limit: 60000.0                           # Upper singal limit to which nonlinearity is applied (ADU)\n')
            f.write('  accuracy: 0.000001                        # Non-linearity accuracy threshold\n')
            f.write('  maxiter: 10                              # Maximum number of iterations to use when applying non-linearity\n')
            f.write('  robberto:  False                         # Use Massimo Robberto type non-linearity coefficients\n')
            f.write('\n')
            f.write('cosmicRay:\n')
            f.write('  path: /ifs/jwst/wit/witserv/data4/nrc/hilbert/simulated_data/cosmic_ray_library/               # Path to CR library\n')
            f.write('  library: SUNMIN    # Type of cosmic rayenvironment (SUNMAX, SUNMIN, FLARE)\n')
            f.write('  scale: 1.5     # Cosmic ray scaling factor\n')
            f.write('  suffix: IPC_NIRCam_{}    # Suffix of library file names\n'.format(input['detector']))
            f.write('  seed: {}                           # Seed for random number generator\n'.format(np.random.randint(1, 2**32-2)))
            f.write('\n')
            f.write('simSignals:\n')
            f.write('  pointsource: {}   #File containing a list of point sources to add (x, y locations and magnitudes)\n'.format(input['{}_ptsrc'.format(catkey)]))   #'point_source']))
            f.write('  psfpath: {}   #Path to PSF library\n'.format(self.psfpath))
            f.write('  psfbasename: {}      #Basename of the files in the psf library\n'.format(self.psfbasename))
            f.write('  psfpixfrac: {}       #Fraction of a pixel between entries in PSF library (e.g. 0.1 = files for PSF centered at 0.1 pixel intervals within pixel)\n'.format(self.psfpixfrac))
            f.write('  psfwfe: {}   #PSF WFE value (predicted or requirements)\n'.format(self.psfwfe))
            f.write('  psfwfegroup: {}      #WFE realization group (0 to 4)\n'.format(self.psfwfegroup))
            f.write('  galaxyListFile: {}    #File containing a list of positions/ellipticities/magnitudes of galaxies to simulate\n'.format(input['{}_galcat'.format(catkey)]))   #'galaxyListFile']))
            f.write('  extended: {}          #Extended emission count rate image file name\n'.format(input['{}_ext'.format(catkey)]))     #'extended']))
            f.write('  extendedscale: {}                          #Scaling factor for extended emission image\n'.format(input['{}_extscl'.format(catkey)]))
            f.write('  extendedCenter: {}                   #x, y pixel location at which to place the extended image if it is smaller than the output array size\n'.format(input['{}_extcent'.format(catkey)]))
            f.write('  PSFConvolveExtended: True #Convolve the extended image with the PSF before adding to the output image (True or False)\n')
            f.write('  movingTargetList: {}          #Name of file containing a list of point source moving targets (e.g. KBOs, asteroids) to add.\n'.format(input['{}_movptsrc'.format(catkey)]))   #'movingTarg']))
            f.write('  movingTargetSersic: {}  #ascii file containing a list of 2D sersic profiles to have moving through the field\n'.format(input['{}_movgal'.format(catkey)]))  #'movingTargSersic']))
            f.write('  movingTargetExtended: {}      #ascii file containing a list of stamp images to add as moving targets (planets, moons, etc)\n'.format(input['{}_movext'.format(catkey)]))  #'movingTargExtended']))
            f.write('  movingTargetConvolveExtended: {}       #convolve the extended moving targets with PSF before adding.\n'.format(input['{}_movconv'.format(catkey)]))
            f.write('  movingTargetToTrack: {} #File containing a single moving target which JWST will track during observation (e.g. a planet, moon, KBO, asteroid)	This file will only be used if mode is set to "moving_target" \n'.format(input['{}_solarsys'.format(catkey)]))  #'movingTargToTrack']))
            f.write('  zodiacal:  None                          #Zodiacal light count rate image file \n')
            f.write('  zodiscale:  1.0                            #Zodi scaling factor\n')
            f.write('  scattered:  None                          #Scattered light count rate image file\n')
            f.write('  scatteredscale: 1.0                        #Scattered light scaling factor\n')
            f.write("  bkgdrate: {}                         #Constant background count rate (ADU/sec/pixel) or 'high','medium','low' similar to what is used in the ETC\n".format(input['{}_bkgd'.format(catkey)]))  #'bkgdrate']))
            f.write('  poissonseed: {}                  #Random number generator seed for Poisson simulation)\n'.format(np.random.randint(1, 2**32-2)))
            f.write('  photonyield: True                         #Apply photon yield in simulation\n')
            f.write('  pymethod: True                            #Use double Poisson simulation for photon yield\n')
            f.write('\n')
            f.write('Telescope:\n')
            f.write('  ra: {}                      # RA of simulated pointing\n'.format(input['ra_ref']))
            f.write('  dec: {}                    # Dec of simulated pointing\n'.format(input['dec_ref']))
            f.write('  rotation: {}                    # y axis rotation (degrees E of N)\n'.format(input['pav3']))
            f.write('  tracking: {}   #Telescope tracking. Can be sidereal or non-sidereal\n'.format(self.tracking))
            f.write('\n')
            f.write('newRamp:\n')
            f.write('  dq_configfile: {}\n'.format(self.dq_init_config))
            f.write('  sat_configfile: {}\n'.format(self.saturation_config))
            f.write('  superbias_configfile: {}\n'.format(self.superbias_config))
            f.write('  refpix_configfile: {}\n'.format(self.refpix_config))
            f.write('  linear_configfile: {}\n'.format(self.linearity_config))
            f.write('\n')
            f.write('Output:\n')
            # f.write('  use_stsci_output_name: {} # Output filename should follow STScI naming conventions (True/False)\n'.format(outtf))
            f.write('  directory: {}  # Output directory\n'.format(self.simdata_output_dir))
            f.write('  file: {}   # Output filename\n'.format(outfile))
            f.write("  datatype: {} # Type of data to save. 'linear' for linearized ramp. 'raw' for raw ramp. 'linear, raw' for both\n".format(self.datatype))
            f.write('  format: DMS          # Output file format Options: DMS, SSR(not yet implemented)\n')
            f.write('  save_intermediates: False   # Save intermediate products separately (point source image, etc)\n')
            f.write('  grism_source_image: {}   # grism\n'.format(input['grism_source_image']))
            f.write('  unsigned: True   # Output unsigned integers? (0-65535 if true. -32768 to 32768 if false)\n')
            f.write('  dmsOrient: True    # Output in DMS orientation (vs. fitswriter orientation).\n')
            f.write('  program_number: {}    # Program Number\n'.format(input['ProposalID']))
            f.write('  title: {}   # Program title\n'.format(input['Title'].replace(':', ', ')))
            f.write('  PI_Name: {}  # Proposal PI Name\n'.format(input['PI_Name']))
            f.write('  Proposal_category: {}  # Proposal category\n'.format(input['Proposal_category']))
            f.write('  Science_category: {}  # Science category\n'.format(input['Science_category']))
            f.write("  observation_number: '{}'    # Observation Number\n".format(input['obs_num']))
            f.write('  observation_label: {}    # User-generated observation Label\n'.format(input['obs_label'].strip()))
            f.write("  visit_number: '{}'    # Visit Number\n".format(input['visit_num']))
            f.write("  visit_group: '{}'    # Visit Group\n".format(input['visit_group']))
            f.write("  visit_id: '{}'    # Visit ID\n".format(input['visit_id']))
            f.write("  sequence_id: '{}'    # Sequence ID\n".format(input['sequence_id']))
            f.write("  activity_id: '{}'    # Activity ID. Increment with each exposure.\n".format(input['act_id']))
            f.write("  exposure_number: '{}'    # Exposure Number\n".format(input['exposure']))
            f.write("  obs_id: '{}'   # Observation ID number\n".format(input['observation_id']))
            f.write("  date_obs: '{}'  # Date of observation\n".format(input['date_obs']))
            f.write("  time_obs: '{}'  # Time of observation\n".format(input['time_obs']))
            f.write("  obs_template: '{}'  # Observation template\n".format(input['obs_template']))
            f.write("  primary_dither_type: {}  # Primary dither pattern name\n".format(input['PrimaryDitherType']))
            f.write("  total_primary_dither_positions: {}  # Total number of primary dither positions\n".format(input['PrimaryDithers']))
            f.write("  primary_dither_position: {}  # Primary dither position number\n".format(np.int(input['primary_dither_num'])))
            f.write("  subpix_dither_type: {}  # Subpixel dither pattern name\n".format(input['SubpixelDitherType']))
            # For WFSS we need to strip out the '-Points' from
            # the number of subpixel positions entry
            try:
                dash = input['SubpixelPositions'].find('-')
                val = input['SubpixelPositions'][0:dash]
            except:
                val = input['SubpixelPositions']
            f.write("  total_subpix_dither_positions: {}  # Total number of subpixel dither positions\n".format(val))
            f.write("  subpix_dither_position: {}  # Subpixel dither position number\n".format(np.int(input['subpix_dither_num'])))
            f.write("  xoffset: {}  # Dither pointing offset in x (arcsec)\n".format(input['idlx']))
            f.write("  yoffset: {}  # Dither pointing offset in y (arcsec)\n".format(input['idly']))
        return yamlout
            
    def reffile_setup(self):
        """Create lists of reference files associate with each detector"""
        self.det_list = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
        sb_dir = os.path.join(self.datadir, 'reference_files/superbias')
        lin_dir = os.path.join(self.datadir, 'reference_files/linearity')
        gain_dir = os.path.join(self.datadir, 'reference_files/gain')
        sat_dir = os.path.join(self.datadir, 'reference_files/saturation')
        ipc_dir = os.path.join(self.datadir, 'reference_files/ipc')
        dist_dir = os.path.join(self.datadir, 'reference_files/distortion')
        pam_dir = os.path.join(self.datadir, 'reference_files/pam')
        rawdark_dir = os.path.join(self.datadir, 'darks/raw')
        lindark_dir = os.path.join(self.datadir, 'darks/linearized')
        self.superbias_list = {}
        self.linearity_list = {}
        self.gain_list = {}
        self.saturation_list = {}
        self.ipc_list = {}
        self.astrometric_list = {}
        self.pam_list = {}
        self.dark_list = {}
        self.lindark_list = {}
        for det in self.det_list:
            sbfiles = glob(os.path.join(sb_dir, '*fits'))
            self.superbias_list[det] = [d for d in sbfiles if 'NRC' + det in d][0]
            linfiles = glob(os.path.join(lin_dir, '*fits'))
            longdet = deepcopy(det)
            if '5' in det:
                longdet = det.replace('5', 'LONG')
            self.linearity_list[det] = [d for d in linfiles if 'NRC' + longdet in d][0]

            gainfiles = glob(os.path.join(gain_dir, '*fits'))
            self.gain_list[det] = [d for d in gainfiles if 'NRC' + det in d][0]

            satfiles = glob(os.path.join(sat_dir, '*fits'))
            self.saturation_list[det] = [d for d in satfiles if 'NRC' + det in d][0]

            ipcfiles = glob(os.path.join(ipc_dir, '*fits'))
            self.ipc_list[det] = [d for d in ipcfiles if 'NRC' + det in d][0]

            distfiles = glob(os.path.join(dist_dir, '*asdf'))
            self.astrometric_list[det] = [d for d in distfiles if 'NRC' + det in d][0]

            pamfiles = glob(os.path.join(pam_dir, '*fits'))
            self.pam_list[det] = [d for d in pamfiles if det in d][0]

            self.dark_list[det] = glob(os.path.join(rawdark_dir, det, '*.fits'))
            self.lindark_list[det] = glob(os.path.join(lindark_dir, det, '*.fits'))

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("--input_xml", help='XML file from APT describing the observations.')
        parser.add_argument("--pointing_file", help='Pointing file from APT describing observations.')
        parser.add_argument("--siaf", help='CSV version of SIAF. Needed only in conjunction with input_xml + pointing.')
        parser.add_argument("--datatype", help='Type of data to save. Can be "linear", "raw" or "linear, raw"', default="linear")
        parser.add_argument("--output_dir", help='Directory into which the yaml files are output', default='./')
        parser.add_argument("--table_file", help='Ascii table containing observation info. Use this or xml + pointing + siaf files.', default=None)
        parser.add_argument("--use_nonstsci_names", help="Use STScI naming convention for output files", action='store_true')
        parser.add_argument("--subarray_def_file", help="Ascii file containing subarray definitions", default='config')
        parser.add_argument("--readpatt_def_file", help='Ascii file containing readout pattern definitions', default='config')
        parser.add_argument("--crosstalk", help="Crosstalk coefficient file", default='config')
        parser.add_argument("--filtpupil_pairs", help="List of paired filter/pupil elements", default='config')
        parser.add_argument("--fluxcal", help="File with zeropoints per filter", default='config')
        parser.add_argument("--dq_init_config", help="DQ Initialization config file", default='config')
        parser.add_argument("--saturation_config", help="Saturation config file", default='config')
        parser.add_argument("--superbias_config", help="Superbias subtraction config file", default='config')
        parser.add_argument("--refpix_config", help="Refpix subtraction config file", default='config')
        parser.add_argument("--linearity_config", help="Linearity config file", default='config')
        parser.add_argument("--observation_table", help="Table file containing epoch start times, telescope roll angles, catalogs for each observation", default=None)
        parser.add_argument("--use_JWST_pipeline", help='True/False', action='store_true')
        parser.add_argument("--use_linearized_darks", help='True/False', action='store_true')
        parser.add_argument("--simdata_output_dir", help='Output directory for simulated exposure files', default='./')
        parser.add_argument("--psfpath", help='Directory containing PSF library',
                            default='/ifs/jwst/wit/nircam/nircam_simulator_data/webbpsf_library')
        parser.add_argument("--psfbasename", help="Basename of the files in the PSF library", default='nircam')
        parser.add_argument("--psfpixfrac", help="Subpixel centering resolution of files in PSF library", default=0.25)
        parser.add_argument("--psfwfe", help="Wavefront error value to use for PSFs", default='predicted')
        parser.add_argument("--psfwfegroup", help="Realization index number for PSF files", default=0)
        parser.add_argument("--resets_bet_ints", help="Number of detector resets between integrations", default=1)
        parser.add_argument("--tracking", help="Type of telescope tracking: 'sidereal' or 'non-sidereal'", default='sidereal')

        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: yaml_generator.py NIRCam_obs.xml NIRCam_obs.pointing'

    input = SimInput()
    parser = input.add_options(usage=usagestring)
    args = parser.parse_args(namespace=input)
    input.reffile_setup()
    input.create_inputs()
