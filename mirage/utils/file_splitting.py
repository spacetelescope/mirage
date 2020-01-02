#! /usr/bin/env python

"""Utility for determining how to split up seed images and dark files
that are too large. This is meant to roughly mimic the way that the
DMS pipeline will do the same job. Files that are spllit will have
``segNNN`` added to their names, where ``NNN`` is a number.
"""
import copy
import numpy as np
from mirage.utils.constants import FILE_SPLITTING_LIMIT


def find_file_splits(xdim, ydim, groups, integrations, frames_per_group=None,
                     pixel_limit=FILE_SPLITTING_LIMIT):
    """Determine the frame and/or integration numbers where a file
    should be split in order to keep the file size reasonable.

    Parameters
    ----------
    xdim : int
        Number of columns in the aperture

    ydim : int
        Number of rows in the aperture

    groups : int
        Number of groups (or frames) in an integration

    integrations : int
        Number of integrations in exposure

    frames_per_group : int
        Number of frames (inculding skipped frames) per
        group

    pixel_limit : int
        Proxy for file size limit. Number of pixel read outs to treat
        as the upper limit to be contained in a single file.

    Returns
    -------
    split : bool
        Whether or not the exposure needs to be split into multiple
        files

    group_list : numpy.ndarray
        1D array listing the beginning group number of each file split

    integration_list : numpy.ndarray
        1d array listing the beginning integration number of each file
        split
    """
    pix_per_group = ydim * xdim
    pix_per_int = groups * pix_per_group
    observation = pix_per_int * integrations

    # Default = no splitting
    split = False
    group_list = np.array([0, groups])
    integration_list = np.array([0, integrations])

    # Check for splitting between groups first
    # i.e. splitting within an integration
    if pix_per_int > pixel_limit:
        split = True
        delta_group = np.int(pixel_limit / pix_per_group)
        # If calculations are being done using frames rather than groups,
        # make sure that we don't split in the middle of a group. Force
        # the splits to be between groups.
        if frames_per_group is not None:
            if delta_group % frames_per_group != 0:
                closest_multiple = delta_group // frames_per_group
                delta_group = closest_multiple * frames_per_group

        group_list = np.arange(0, groups, delta_group).astype(np.int)
        group_list = np.append(group_list, groups)
        print('Splitting within each integration:')
        integration_list = np.arange(integrations + 1).astype(np.int)
        print('integration_list: ', integration_list)
        print('group_list: ', group_list)
    elif observation > pixel_limit:
        split = True
        print('Splitting by integration:')
        group_list = np.array([0, groups])
        delta_int = np.int(pixel_limit / pix_per_int)
        integration_list = np.arange(0, integrations, delta_int).astype(np.int)
        integration_list = np.append(integration_list, integrations)
        print('integration_list: ', integration_list)
        print('group_list: ', group_list)

    return split, group_list, integration_list


class SplitFileMetaData:
    def __init__(self, integration_splits, frame_splits, DMS_integration_splits, DMS_group_splits,
                 frames_per_integration, frames_per_group, frametime):
        """Calculate metadata values related to the position of the split files'
        data within the segment and exposure. These data will be needed when
        reconstructing `parts` into segments later

        Parameters
        ----------
        integration_splits : list
            List of the beginning integration numbers for file splitting
            when split using Mirage's parameters (i.e. including 'parts')
            This list is generated using ``find_file_splits``

        frame_splits : list
            List of the beginning frame numbers for file splitting when
            split using Mirage's parameters (i.e. including 'parts')
            This list is generated using ``find_file_splits``

        DMS_integration_splits : list
            List of the beginning integration numbers for file splitting
            when split using DMS parameters (i.e. no splits within an
            integration, splitting on groups rather than frames)
            This list is generated using ``find_file_splits``

        DMS_group_splits : list
            List of the beginning group numbers for file splitting
            when split using DMS parameters (i.e. no splits within an
            integration, splitting on groups rather than frames)
            This should always be [0, last_group]
            This list is generated using ``find_file_splits``

        frames_per_integration : int
            Number of frames per integration of the simulated data

        frames_per_group : int
            Number of frames per group in the simulated data

        frametime : float
            Exposure time of one frame (e.e. detector readout) of the
            simulated data
        """
        self.total_frames = []
        self.total_ints = []
        self.time_start = []
        self.frame_start = []
        self.segment_part_number = []
        self.segment_ints = []
        self.segment_frames = []
        self.segment_frame_start_number = []
        self.segment_int_start_number = []
        self.part_int_start_number = []
        self.part_frame_start_number = []
        self.segment_number = []
        segment_starting_int_number = 0
        segment_part_number = 0
        segment_frame_start_number = 0
        segment_int_start_number = 0

        i = 1
        previous_segment = 0
        for int_start in integration_splits[:-1]:
            int_end = integration_splits[i]

            j = 1
            for initial_frame in frame_splits[:-1]:
                frame_start = int_start * (frames_per_integration + 1) + initial_frame
                frame_end = frame_start + (frame_splits[j] - initial_frame)
                time_start = frame_start * frametime
                time_end = frame_end * frametime
                total_ints = int_end - int_start
                # total_frames should be the total number of frames
                # in the chunk to be worked on, including resets
                total_frames = (frame_end - frame_start) * total_ints
                if total_ints > 1:
                    total_frames += total_ints - 1

                #print('int_start, int_end, i:', int_start, int_end, i)
                #print('initial_frame, frame_start, frame_end, total_frames:', initial_frame, frame_start, frame_end, total_frames)

                self.total_frames.append(total_frames)
                self.total_ints.append(total_ints)
                self.time_start.append(time_start)
                self.frame_start.append(frame_start)

                segment_number = np.where(int_end <= DMS_integration_splits)[0][0]
                self.segment_ints.append(DMS_integration_splits[segment_number] - DMS_integration_splits[segment_number - 1])
                self.segment_frames.append((DMS_group_splits[1] - DMS_group_splits[0]) * frames_per_group)

                if segment_number == previous_segment:
                    segment_part_number += 1
                    part_int_start_number = int_start - segment_starting_int_number
                    part_frame_start_number = initial_frame
                else:
                    segment_part_number = 1
                    previous_segment = copy.deepcopy(segment_number)
                    segment_frame_start_number = initial_frame
                    segment_int_start_number = int_start
                    part_int_start_number = 0
                    part_frame_start_number = 0
                    segment_starting_int_number = copy.deepcopy(int_start)

                self.segment_part_number.append(segment_part_number)
                self.segment_frame_start_number.append(segment_frame_start_number)
                self.segment_int_start_number.append(segment_int_start_number)
                self.part_int_start_number.append(part_int_start_number)
                self.part_frame_start_number.append(part_frame_start_number)
                self.segment_number.append(segment_number)

                #print("\n\nSegment Numbers:")
                #print(int_start, initial_frame)
                #print(self.total_frames)
                #print(self.segment_number, self.segment_part_number)
                #print(self.segment_ints, self.segment_frames)
                #print(self.segment_int_start_number, self.segment_frame_start_number)
                #print(self.part_int_start_number, self.part_frame_start_number)
                #print('\n\n')

                j += 1

            i += 1
