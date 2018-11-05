# ! /usr/bin/env python

'''
Given APT output files, read in data relevant to the data simulator,
organize, and create input files for the simulator.

Inputs:

xml file - Name of xml file exported from APT.
pointing file - Name of associated pointing file exported from APT.

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


HISTORY:

July 2017 - V0: Initial version. Bryan Hilbert
Feb 2018  - V1: Updated to work for multiple filter pairs per observation
            Lauren Chambers
August 2018 - V2: Replaced manual Ra, Dec calculations with pysiaf functionality
October 2018 - Major modifications to read programs of all science instruments and parallels
               Johannes Sahlmann
'''
import copy
import os
import re
import argparse

from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
from pysiaf import rotations
import yaml

from . import read_apt_xml
from ..utils import siaf_interface


class AptInput:
    """Summary

    Attributes:
        exposure_tab (TYPE): Description
        input_xml (str): Description
        observation_list_file (str): Description
        obstab (TYPE): Description
        output_csv (TYPE): Description
        pointing_file (str): Description
    """

    def __init__(self, input_xml=None, pointing_file=None):

        self.input_xml = input_xml
        self.pointing_file = pointing_file

        self.output_csv = None
        self.observation_list_file = None

    def add_epochs(self, intab):
        """NOT CURRENTLY USED"""
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

    def add_observation_info(self, intab):
        """Add information about each observation.

        Catalog names, dates, PAV3 values, etc., which are retrieved from the observation list
        yaml file.

        Parameters
        ----------
        intab : obj
            astropy.table.Table containing exposure information

        Returns
        -------
        intab : obj
            Updated table with information from the observation list
            yaml file added.

        """
        with open(self.observation_list_file, 'r') as infile:
            self.obstab = yaml.load(infile)

        OBSERVATION_LIST_FIELDS = 'Date PAV3 Filter PointSourceCatalog GalaxyCatalog ' \
                                  'ExtendedCatalog ExtendedScale ExtendedCenter MovingTargetList ' \
                                  'MovingTargetSersic MovingTargetExtended ' \
                                  'MovingTargetConvolveExtended MovingTargetToTrack ' \
                                  'BackgroundRate DitherIndex'.split()

        nircam_mapping = {'ptsrc': 'PointSourceCatalog',
                          'galcat': 'GalaxyCatalog',
                          'ext': 'ExtendedCatalog',
                          'extscl': 'ExtendedScale',
                          'extcent': 'ExtendedCenter',
                          'movptsrc': 'MovingTargetList',
                          'movgal': 'MovingTargetSersic',
                          'movext': 'MovingTargetExtended',
                          'movconv': 'MovingTargetConvolveExtended',
                          'solarsys': 'MovingTargetToTrack',
                          'bkgd': 'BackgroundRate',
                          }

        unique_instrument_names = [name.lower() for name in np.unique(intab['Instrument'])]

        # initialize dictionary keys
        for key in OBSERVATION_LIST_FIELDS:
            intab[key] = []

        if 'nircam' in unique_instrument_names:
            for channel in ['SW', 'LW']:
                for name, item in nircam_mapping.items():
                    key = '{}_{}'.format(channel.lower(), name)
                    intab[key] = []

        # loop over entries in input dictionary
        for index, instrument in enumerate(intab['Instrument']):
            instrument = instrument.lower()

            # retrieve corresponding entry from observation list
            entry = get_entry(self.obstab, intab['entry_number'][index])

            if instrument == 'nircam':
                # keep the number of entries in the dictionary consistent
                for key in OBSERVATION_LIST_FIELDS:
                    if key == 'Date':
                        value = entry[key].strftime('%Y-%m-%d')
                    elif key in ['PAV3', 'Instrument']:
                        value = str(entry[key])
                    else:
                        value = str(None)

                    intab[key].append(value)

                for channel in ['SW', 'LW']:
                    for name, item in nircam_mapping.items():
                        key = '{}_{}'.format(channel.lower(), name)
                        if item in 'ExtendedScale ExtendedCenter MovingTargetConvolveExtended BackgroundRate'.split():
                            intab[key].append(entry['FilterConfig'][channel][item])
                        else:
                            intab[key].append(self.full_path(entry['FilterConfig'][channel][item]))

            else:
                for key in OBSERVATION_LIST_FIELDS:
                    if key == 'Date':
                        value = entry[key].strftime('%Y-%m-%d')
                    else:
                        value = str(entry[key])

                    intab[key].append(value)

                # keep the number of entries in the dictionary consistent
                if 'nircam' in unique_instrument_names:
                    for channel in ['SW', 'LW']:
                        for name, item in nircam_mapping.items():
                            key = '{}_{}'.format(channel.lower(), name)
                            intab[key].append(str(None))

        intab['epoch_start_date'] = intab['Date']

        return intab

    def base36encode(self, integer):
        """
        Translate a base 10 integer to base 36

        Parameters
        ----------
        integer : int
            a base 10 integer

        Returns
        -------
        integer : int
            The integer translated to base 36
        """
        chars, encoded = '0123456789abcdefghijklmnopqrstuvwxyz', ''

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            encoded = chars[remainder] + encoded

        return encoded.zfill(2)

    def combine_dicts(self, dict1, dict2):
        """Combine two dictionaries into a single dictionary.

        Parameters
        ----------
        dict1 : dict
            dictionary
        dict2 : dict
            dictionary

        Returns
        -------
        combined : dict
            Combined dictionary
        """
        combined = dict1.copy()
        combined.update(dict2)
        return combined

    def create_input_table(self, verbose=False):
        """

        Expansion for dithers is done upstream.

        Parameters
        ----------
        verbose

        Returns
        -------

        """
        # Expand paths to full paths
        # self.input_xml = os.path.abspath(self.input_xml)
        # self.pointing_file = os.path.abspath(self.pointing_file)
        if self.output_csv is not None:
            self.output_csv = os.path.abspath(self.output_csv)
        if self.observation_list_file is not None:
            self.observation_list_file = os.path.abspath(self.observation_list_file)

        # main_dir = os.path.split(self.input_xml)[0]

        # if APT.xml content has already been generated during observation list creation
        # (generate_observationlist.py) load it here

        if self.apt_xml_dict is None:
            raise RuntimeError('self.apt_xml_dict is not defined')
            # tab = self.apt_xml_dict
        # if self.apt_xml_dict is not None:
        # else:
        #     tab = self.apt_xml_dict
        # else:
        #     # Read in xml file
        #     readxml_obj = read_apt_xml.ReadAPTXML()
        #     tab = readxml_obj.read_xml(self.input_xml)


        # This affects on NIRCam exposures (right?)
        # If the number of dithers is set to '3TIGHT'
        # (currently only used in NIRCam)
        # remove 'TIGHT' from the entries and leave
        # only the number behind
        # tight = [True if 'TIGHT' in str(val) else False for val in tab['PrimaryDithers']]
        # if np.any(tight):
        #     tab = self.tight_dithers(tab)

        # if verbose:
        # for key in tab.keys():
        #     print('{:<25}: number of elements is {:>5}'.format(key, len(tab[key])))

        # xmltab = tab

        # Read in the pointing file and produce dictionary
        pointing_dictionary = self.get_pointing_info(self.pointing_file, propid=self.apt_xml_dict['ProposalID'][0])

        # Check that the .xml and .pointing files agree
        assert len(self.apt_xml_dict['ProposalID']) == len(pointing_dictionary['obs_num']),\
            'Inconsistent table size from XML file ({}) and pointing file ({}). Something was not processed correctly in apt_inputs.'.format(len(self.apt_xml_dict['ProposalID']), len(pointing_dictionary['obs_num']))

        # Combine the dictionaries
        observation_dictionary = self.combine_dicts(self.apt_xml_dict, pointing_dictionary)

        # Add epoch and catalog information
        observation_dictionary = self.add_observation_info(observation_dictionary)

        # if verbose:
        #     print('Summary of observation dictionary:')
        #     for key in observation_dictionary.keys():
        #         print('{:<25}: number of elements is {:>5}'.format(key, len(observation_dictionary[key])))

        # NIRCam case: Expand for detectors. Create one entry in each list for each
        # detector, rather than a single entry for 'ALL' or 'BSALL'

        # test if Module is always 'None', i.e. when NIRCam is not used
        if observation_dictionary['Module'].count('None') == len(observation_dictionary['Module']):
            # set detector key
            observation_dictionary['detector'] = []
            for i, instrument in enumerate(observation_dictionary['Instrument']):
                if instrument.lower() == 'niriss':
                    observation_dictionary['detector'].append('NIS')
                elif instrument.lower() == 'nirspec':
                    observation_dictionary['detector'].append('NRS')
                elif instrument.lower() == 'nirspec':
                    if 'NRS1' in observation_dictionary['aperture'][i]:
                        observation_dictionary['detector'].append('NRS1')
                    elif 'NRS2' in observation_dictionary['aperture'][i]:
                        observation_dictionary['detector'].append('NRS1')
                elif instrument.lower() == 'fgs':
                    if 'FGS1' in observation_dictionary['aperture'][i]:
                        observation_dictionary['detector'].append('G1')
                    elif 'FGS2' in observation_dictionary['aperture'][i]:
                        observation_dictionary['detector'].append('G2')
                elif instrument.lower() == 'miri':
                        observation_dictionary['detector'].append('MIR')
            self.exposure_tab = observation_dictionary
        else:
            self.exposure_tab = self.expand_for_detectors(observation_dictionary)

        if verbose:
            for key in self.exposure_tab.keys():
                print('{:>20} has {:>10} items'.format(key, len(self.exposure_tab[key])))

        detectors_file = os.path.join(self.output_dir, 'expand_for_detectors.csv')

        ascii.write(Table(self.exposure_tab), detectors_file, format='csv', overwrite=True)
        print('Wrote exposure table to {}'.format(detectors_file))

        # Calculate the correct V2, V3 and RA, Dec for each exposure/detector
        self.ra_dec_update()

        # Output to a csv file.
        if self.output_csv is None:
            indir, infile = os.path.split(self.input_xml)
            self.output_csv = os.path.join(self.output_dir, 'Observation_table_for_' + infile.split('.')[0] + '.csv')
        ascii.write(Table(self.exposure_tab), self.output_csv, format='csv', overwrite=True)
        print('csv exposure list written to {}'.format(self.output_csv))

    def expand_for_detectors(self, obstab):
        """
        Expand dictionary to have one entry per detector, rather than the
        one line per module that is in the input

        Parameters
        ----------
        obstab : dict
            dictionary containing one entry per module

        Returns
        -------
        finaltab : dict
            dictionary expanded to have one entry per detector
        """
        finaltab = {}
        for key in obstab:
            finaltab[key] = []
        finaltab['detector'] = []

        for index, instrument in enumerate(obstab['Instrument']):
            instrument = instrument.lower()
            if instrument == 'nircam':

                # Determine module of the observation
                module = obstab['Module'][index]
                if module == 'ALL':
                    detectors = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
                elif module == 'A':
                    detectors = ['A1', 'A2', 'A3', 'A4', 'A5']
                elif module == 'B':
                    detectors = ['B1', 'B2', 'B3', 'B4', 'B5']
                elif 'A3' in module:
                    detectors = ['A3']
                elif 'B4' in module:
                    detectors = ['B4']
                elif 'DHSPIL' in module:
                    if module[-1] == 'A':
                        detectors = ['A3']
                    elif module[-1] == 'B':
                        detectors = ['B4']
                    else:
                        ValueError('Unknown module {}'.format(module))
                else:
                    raise ValueError('Unknown module {}'.format(module))

                n_detectors = len(detectors)
                for key in obstab:
                    finaltab[key].extend(([obstab[key][index]] * n_detectors))
                finaltab['detector'].extend(detectors)

        return finaltab

    def extract_value(self, line):
        """Extract text from xml line

        Parameters
        ----------
        line : str
            Line from xml file

        Returns
        -------
        line : str
            Text between > and < in the input line
        """
        gt = line.find('>')
        lt = line.find('<', gt)
        return line[gt + 1:lt]

    def full_path(self, in_path):
        """
        If the input path is not None, expand
        any environment variables and make an
        absolute path. Return the updated path.

        Parameters
        ----------
        in_path : str
            Path to be expanded

        Returns
        -------
        in_path : str
            Expanded, absolute path
        """
        if in_path.lower() == 'none':
            return in_path
        else:
            return os.path.abspath(os.path.expandvars(in_path))

    def get_pointing_info(self, file, propid=0, verbose=False):
        """Read in information from APT's pointing file.

        Parameters
        ----------
        file : str
            Name of APT-exported pointing file to be read
        propid : int
            Proposal ID number (integer). This is used to
            create various ID fields

        Returns
        -------
        pointing : dict
            Dictionary of pointing-related information

        TODO
        ----
            extract useful information from header?
            check visit numbers
            set parallel proposal number correctly

        """
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

                #skip comments and new lines
                if (line[0] == '#') or (line in ['\n']) or ('=====' in line):
                    continue
                # extract proposal ID
                elif line.split()[0] == 'JWST':
                    propid_header = line.split()[7]
                    try:
                        propid = np.int(propid_header)
                    except ValueError:
                        #adopt value passed to function
                        pass
                    if verbose:
                        print('Extracted proposal ID {}'.format(propid))
                    continue


                elif (len(line) > 1):
                    elements = line.split()

                    # Look for lines that give visit/observation numbers
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

                    skip = False

                    if line[0:2] == '**':
                        v = elements[2]
                        obsnum, visitnum = v.split(':')
                        obsnum = str(obsnum).zfill(3)
                        visitnum = str(visitnum).zfill(3)
                        if (skip is True) and (verbose):
                            print('Skipping observation {} ({})'.format(obsnum, obslabel))

                    try:
                        # Skip the line at the beginning of each
                        # visit that gives info on the target,
                        # but is not actually an observation
                        # These lines have 'Exp' values of 0,
                        # while observations have a value of 1
                        # (that I've seen so far)

                        if ((np.int(elements[1]) > 0) & ('NRC' in elements[4]
                                                         or 'NIS' in elements[4]
                                                         or 'FGS' in elements[4]
                                                         or 'NRS' in elements[4]
                                                         or 'MIR' in elements[4])
                            ):
                            if (elements[18] == 'PARALLEL') and ('MIRI' in elements[4]):
                                skip = True

                            if skip:
                                act_counter += 1
                                continue
                            act = self.base36encode(act_counter)
                            activity_id.append(act)
                            observation_label.append(obslabel)
                            observation_number.append(obsnum)
                            visit_number.append(visitnum)
                            prop_5digit = "{0:05d}".format(int(propid))
                            vid = "{}{}{}".format(prop_5digit, obsnum, visitnum)
                            visit_id.append(vid)
                            # Visit group hard coded to 1. It's not clear how APT divides visits up into visit
                            # groups. For now just keep everything in a single visit group.
                            vgrp = '01'
                            visit_grp.append(vgrp)
                            # Parallel sequence id hard coded to 1 (Simulated instrument as prime rather than
                            # parallel) at the moment. Future improvements may allow the proper sequence
                            # number to be constructed.
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
                            if elements[18] == 'PARALLEL':
                                ddist.append(None)
                            else:
                                ddist.append(np.float(elements[21]))
                            # For the moment we assume that the instrument being simulated is not being
                            # run in parallel, so the parallel proposal number will be all zeros,
                            # as seen in the line below.
                            observation_id.append("V{}P{}{}{}{}".format(vid, '00000000', vgrp, seq, act))
                            act_counter += 1

                    except ValueError as e:
                        if verbose:
                            print('Skipping line:\n{}\nproducing error:\n{}'.format(line, e))
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

    def tight_dithers(self, input_dict):
        """
        In NIRCam, when the 'FULL' dither pattern is
        used, it is possible to set the number of primary
        dithers to '3TIGHT' rather than just a number
        (e.g. '3'). If the number of dithers is set to '3TIGHT'
        remove 'TIGHT' from the entries and leave
        only the number behind.

        Parameters
        ----------

        input_dict : dict
           Dictionary where each key points to a list containing
           observation details. For example, input_dict['PrimarDither']
           is a list of the number of primary dithers for all observations

        Returns
        -------

        input_dict : dict
            Updated dictionary where 'TIGHT' has been removed from
            PrimaryDither list
        """
        inlist = input_dict['PrimaryDithers']
        modlist = [v if 'TIGHT' not in v else v.strip('TIGHT') for v in inlist]
        input_dict['PrimaryDithers'] = modlist
        return input_dict

    def ra_dec_update(self, verbose=False):
        """Given the V2, V3 values for the reference pixels associated
        with detector apertures, calculate corresponding RA, Dec.
        """
        aperture_ra = []
        aperture_dec = []

        for i in range(len(self.exposure_tab['Module'])):
            siaf_instrument = self.exposure_tab["Instrument"][i]
            if siaf_instrument == 'NIRSPEC':
                siaf_instrument = 'NIRSpec'
            aperture_name = self.exposure_tab['aperture'][i]
            pointing_ra = np.float(self.exposure_tab['ra'][i])
            pointing_dec = np.float(self.exposure_tab['dec'][i])
            pointing_v2 = np.float(self.exposure_tab['v2'][i])
            pointing_v3 = np.float(self.exposure_tab['v3'][i])

            if 'pav3' in self.exposure_tab.keys():
                pav3 = np.float(self.exposure_tab['pav3'][i])
            else:
                pav3 = np.float(self.exposure_tab['PAV3'][i])

            telescope_roll = pav3

            aperture, local_roll, attitude_matrix, fullframesize, subarray_boundaries = \
                siaf_interface.get_siaf_information(
                siaf_instrument, aperture_name, pointing_ra, pointing_dec, telescope_roll,
                    v2_arcsec=pointing_v2, v3_arcsec=pointing_v3)

            # calculate RA, Dec of reference location for the detector
            ra, dec = rotations.pointing(attitude_matrix, aperture.V2Ref, aperture.V3Ref)
            aperture_ra.append(ra)
            aperture_dec.append(dec)

        self.exposure_tab['ra_ref'] = aperture_ra
        self.exposure_tab['dec_ref'] = aperture_dec

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("input_xml", help='XML file from APT describing the observations.')
        parser.add_argument("pointing_file", help='Pointing file from APT describing observations.')
        parser.add_argument("--output_csv", help="Name of output CSV file containing list of observations.", default=None)
        parser.add_argument("--observation_list_file", help='Ascii file containing a list of observations, times, and roll angles, catalogs', default=None)
        return parser


def get_entry(dict, entry_number):
    """Return a numbered entry from a dictionary that corresponds to the observataion_list.yaml.

    Parameters
    ----------
    dict
    entry_number

    Returns
    -------

    """
    entry_key = 'EntryNumber{}'.format(entry_number)
    for key, observation in dict.items():
        if entry_key in observation.keys():
            return observation[entry_key]


# if __name__ == '__main__':
#
#     usagestring = 'USAGE: apt_inputs.py NIRCam_obs.xml NIRCam_obs.pointing'
#
#     input = AptInput()
#     parser = input.add_options(usage=usagestring)
#     args = parser.parse_args(namespace=input)
#     input.create_input_table()
