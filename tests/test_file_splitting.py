'''Define unit tests that cover how large output exposures will be split
between files.

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
'''
import numpy as np

from mirage.utils import file_splitting

def test_find_file_splits():
    """Test the function that decides where to split exposures that are
    too large
    """
    ncols = 2048
    nrows = 2048
    ngroups = 5
    nints = 1000
    limit = 2048 * 2048 * 38

    # First try splitting using groups (rather than frames)
    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, ngroups, nints,
                                                                   frames_per_group=None,
                                                                   pixel_limit=limit)
    assert is_split
    assert np.all(np.array(grp_list) == np.array([0, 5]))
    assert np.all(np.array(int_list) == np.append(np.arange(0, 1000, 7), np.array([1000])))

    # Case where this is no splitting
    nints = 5
    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, ngroups, nints,
                                                                   frames_per_group=None,
                                                                   pixel_limit=limit)
    assert is_split is False
    assert np.all(np.array(grp_list) == np.array([0, 5]))
    assert np.all(np.array(int_list) == np.array([0, 5]))

    # Now split assuming frames. There should be no splits within a group
    # In this case, ngroup is treated as the number of frames.
    frm_per_group = 5  # BRIGHT readout pattern
    nints = 1000
    nframes = ngroups * frm_per_group

    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, nframes, nints,
                                                                   frames_per_group=frm_per_group,
                                                                   pixel_limit=limit)
    assert is_split
    assert np.all(np.array(grp_list) == np.array([0, 25]))
    assert np.all(np.array(int_list) == np.arange(0, 1001))

    # Case where the splitting will be between groups
    frm_per_group = 10  # Medium readout pattern
    nframes = ngroups * frm_per_group
    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, nframes, nints,
                                                                   frames_per_group=frm_per_group,
                                                                   pixel_limit=limit)
    assert is_split
    assert np.all(np.array(grp_list) == np.array([0, 30, 50]))
    assert np.all(np.array(int_list) == np.arange(0, 1001))

def test_SplitFileMetaData():
    """Test that the correct metadata are produced for given integration and
    group splitting lists
    """
    ncols = 2048
    nrows = 2048
    ngroups = 5
    nints = 1000
    limit = 2048 * 2048 * 38

    # Most common case. RAPID TSO observation, so no need to worry about
    # frames vs groups
    frm_per_int = ngroups
    frame_time = 10.

    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, ngroups, nints,
                                                                   frames_per_group=None,
                                                                   pixel_limit=limit)

    param = file_splitting.SplitFileMetaData(int_list, grp_list, int_list, grp_list, frm_per_int, 1, frame_time)

    compare_frames = np.zeros_like(param.total_frames) + 41
    compare_frames[-1] = 35
    compare_ints = np.zeros_like(param.total_ints) + 7
    compare_ints[-1] = 6

    assert np.all(param.total_frames == compare_frames)
    assert np.all(param.total_ints == compare_ints)
    assert np.all(param.time_start == np.arange(0, 420. * len(param.time_start), 420.))
    assert np.all(param.frame_start == np.arange(0, 42 * len(param.frame_start), 42))
    assert np.all(param.segment_ints == param.total_ints)
    assert np.all(param.segment_frames == np.zeros(len(param.segment_frames)).astype(int) + 5)
    assert np.all(param.segment_part_number == np.zeros(len(param.segment_part_number)).astype(int) + 1)
    assert np.all(param.segment_frame_start_number == np.zeros(len(param.segment_frame_start_number)).astype(int))
    assert np.all(param.segment_int_start_number == np.arange(0, 995, 7).astype(int))
    assert np.all(param.part_int_start_number == np.zeros(len(param.part_int_start_number)).astype(int))
    assert np.all(param.part_frame_start_number == np.zeros(len(param.part_frame_start_number)).astype(int))
    assert np.all(param.segment_number == np.arange(1, 144).astype(int))

    # More complex case, using a non-RAPID readout pattern. Let's assume
    # one of the BRIGHT patterns, which have 2 frames per group
    # FUTURE WORK: look into the irregular splitting here. It's not clear
    # why some segments are split into 2 parts while others are split into
    # 3 parts.
    frm_per_group = 2
    frm_per_int = ngroups * frm_per_group

    # Mirage splits
    is_split, grp_list, int_list = file_splitting.find_file_splits(ncols, nrows, ngroups * frm_per_group, nints,
                                                                   frames_per_group=frm_per_group,
                                                                   pixel_limit=limit)

    # DMS splits
    is_split, dms_grp_list, dms_int_list = file_splitting.find_file_splits(ncols, nrows, ngroups, nints,
                                                                   frames_per_group=None,
                                                                   pixel_limit=limit)

    param = file_splitting.SplitFileMetaData(int_list, grp_list, dms_int_list, dms_grp_list, frm_per_int,
                                             frm_per_group, frame_time)

    compare_frames = np.zeros_like(param.total_frames) + 32
    compare_frames[-1] = 10

    compare_ints = np.zeros_like(param.total_ints) + 3
    compare_ints[-1] = 1

    compare_seg_ints = np.zeros_like(param.segment_ints) + 7
    compare_seg_ints[-3:] = 6

    compare_part_num = np.array([1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1,
                                 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2,
                                 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1,
                                 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2,
                                 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1,
                                 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2,
                                 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3,
                                 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1,
                                 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2,
                                 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1,
                                 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2,
                                 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1,
                                 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2,
                                 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3,
                                 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 1,
                                 2, 1, 2, 3])

    compare_seg_int_start = np.array([0, 0, 6, 6, 12, 12, 12, 21, 21, 27, 27, 33, 33, 33, 42, 42, 48,
                                      48, 54, 54, 54, 63, 63, 69, 69, 75, 75, 75, 84, 84, 90, 90, 96,
                                      96, 96, 105, 105, 111, 111, 117, 117, 117, 126, 126, 132, 132,
                                      138, 138, 138, 147, 147, 153, 153, 159, 159, 159, 168, 168, 174,
                                      174, 180, 180, 180, 189, 189, 195, 195, 201, 201, 201, 210, 210,
                                      216, 216, 222, 222, 222, 231, 231, 237, 237, 243, 243, 243, 252,
                                      252, 258, 258, 264, 264, 264, 273, 273, 279, 279, 285, 285, 285,
                                      294, 294, 300, 300, 306, 306, 306, 315, 315, 321, 321, 327, 327,
                                      327, 336, 336, 342, 342, 348, 348, 348, 357, 357, 363, 363, 369,
                                      369, 369, 378, 378, 384, 384, 390, 390, 390, 399, 399, 405, 405,
                                      411, 411, 411, 420, 420, 426, 426, 432, 432, 432, 441, 441, 447,
                                      447, 453, 453, 453, 462, 462, 468, 468, 474, 474, 474, 483, 483,
                                      489, 489, 495, 495, 495, 504, 504, 510, 510, 516, 516, 516, 525,
                                      525, 531, 531, 537, 537, 537, 546, 546, 552, 552, 558, 558, 558,
                                      567, 567, 573, 573, 579, 579, 579, 588, 588, 594, 594, 600, 600,
                                      600, 609, 609, 615, 615, 621, 621, 621, 630, 630, 636, 636, 642,
                                      642, 642, 651, 651, 657, 657, 663, 663, 663, 672, 672, 678, 678,
                                      684, 684, 684, 693, 693, 699, 699, 705, 705, 705, 714, 714, 720,
                                      720, 726, 726, 726, 735, 735, 741, 741, 747, 747, 747, 756, 756,
                                      762, 762, 768, 768, 768, 777, 777, 783, 783, 789, 789, 789, 798,
                                      798, 804, 804, 810, 810, 810, 819, 819, 825, 825, 831, 831, 831,
                                      840, 840, 846, 846, 852, 852, 852, 861, 861, 867, 867, 873, 873,
                                      873, 882, 882, 888, 888, 894, 894, 894, 903, 903, 909, 909, 915,
                                      915, 915, 924, 924, 930, 930, 936, 936, 936, 945, 945, 951, 951,
                                      957, 957, 957, 966, 966, 972, 972, 978, 978, 978, 987, 987, 993,
                                      993, 993])

    compare_part_int_start = np.array([0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0,
                                       3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3,
                                       0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0,
                                       3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3,
                                       0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0,
                                       3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3,
                                       6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6,
                                       0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0,
                                       3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3,
                                       0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0,
                                       3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3,
                                       0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0,
                                       3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3,
                                       6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6,
                                       0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0, 3, 0, 3, 0, 3, 6, 0,
                                       3, 0, 3, 6])

    compare_seg_num = np.array([1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 9, 10, 10, 11,
                                11, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19,
                                19, 20, 20, 21, 21, 21, 22, 22, 23, 23, 24, 24, 24, 25, 25, 26, 26, 27, 27,
                                27, 28, 28, 29, 29, 30, 30, 30, 31, 31, 32, 32, 33, 33, 33, 34, 34, 35, 35,
                                36, 36, 36, 37, 37, 38, 38, 39, 39, 39, 40, 40, 41, 41, 42, 42, 42, 43, 43,
                                44, 44, 45, 45, 45, 46, 46, 47, 47, 48, 48, 48, 49, 49, 50, 50, 51, 51, 51,
                                52, 52, 53, 53, 54, 54, 54, 55, 55, 56, 56, 57, 57, 57, 58, 58, 59, 59, 60,
                                60, 60, 61, 61, 62, 62, 63, 63, 63, 64, 64, 65, 65, 66, 66, 66, 67, 67, 68,
                                68, 69, 69, 69, 70, 70, 71, 71, 72, 72, 72, 73, 73, 74, 74, 75, 75, 75, 76,
                                76, 77, 77, 78, 78, 78, 79, 79, 80, 80, 81, 81, 81, 82, 82, 83, 83, 84, 84,
                                84, 85, 85, 86, 86, 87, 87, 87, 88, 88, 89, 89, 90, 90, 90, 91, 91, 92, 92,
                                93, 93, 93, 94, 94, 95, 95, 96, 96, 96, 97, 97, 98, 98, 99, 99, 99, 100, 100,
                                101, 101, 102, 102, 102, 103, 103, 104, 104, 105, 105, 105, 106, 106, 107, 107,
                                108, 108, 108, 109, 109, 110, 110, 111, 111, 111, 112, 112, 113, 113, 114, 114,
                                114, 115, 115, 116, 116, 117, 117, 117, 118, 118, 119, 119, 120, 120, 120, 121,
                                121, 122, 122, 123, 123, 123, 124, 124, 125, 125, 126, 126, 126, 127, 127, 128,
                                128, 129, 129, 129, 130, 130, 131, 131, 132, 132, 132, 133, 133, 134, 134, 135,
                                135, 135, 136, 136, 137, 137, 138, 138, 138, 139, 139, 140, 140, 141, 141, 141,
                                142, 142, 143, 143, 143])

    assert np.all(param.total_frames == compare_frames)
    assert np.all(param.total_ints == compare_ints)
    assert np.all(param.time_start == np.arange(0., 330. * len(param.time_start), 330.))
    assert np.all(param.frame_start == np.arange(0, 33 * len(param.frame_start), 33))
    assert np.all(param.segment_ints == compare_seg_ints)
    assert np.all(param.segment_frames == np.zeros(len(param.segment_frames)).astype(int) + 10)
    assert np.all(param.segment_part_number == compare_part_num)
    assert np.all(param.segment_frame_start_number == np.zeros(len(param.segment_frame_start_number)))
    assert np.all(param.segment_int_start_number == compare_seg_int_start)
    assert np.all(param.part_int_start_number == compare_part_int_start)
    assert np.all(param.part_frame_start_number == np.zeros(len(param.part_frame_start_number)))
    assert np.all(param.segment_number == compare_seg_num)
