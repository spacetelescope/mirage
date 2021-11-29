"""Unit tests for optical ghost-related code

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
"""

from astropy.table import Table
import numpy as np
import os
import pytest

from mirage.ghosts import niriss_ghosts
from mirage.reference_files.downloader import download_file
from mirage.utils.constants import DEFAULT_NIRISS_PTSRC_GHOST_FILE, NIRISS_GHOST_GAP_FILE, \
                                   NIRISS_GHOST_GAP_URL


def test_determine_ghost_stamp_filename():
    # Table with niriss_ghost_stamp column - cases that should work
    tab = Table()
    source_types = ['point_source', 'point_source', 'galaxies', 'extended']
    stamp_files = ['my_stamp.fits', None, 'my_galstamp.fits', 'my_extstamp.fits']
    tab['niriss_ghost_stamp'] = stamp_files
    expected = ['my_stamp.fits', DEFAULT_NIRISS_PTSRC_GHOST_FILE, 'my_galstamp.fits', 'my_extstamp.fits']

    for row, source, expectation in zip(tab, source_types, expected):
        outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)
        assert outfile == expectation

    # Table with niriss_ghost_stamp column - case that should raise exception
    source_types = ['galaxy']
    tab = Table()
    exception_files = [None]
    tab['niriss_ghost_stamp'] = exception_files

    for row, source in zip(tab, source_types):
        with pytest.raises(ValueError):
            outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)

    # Table with niriss_ghost_stamp column - case where outfile should be None
    source_types = ['extended']
    tab = Table()
    exception_files = [None]
    tab['niriss_ghost_stamp'] = exception_files

    for row, source in zip(tab, source_types):
        outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)
        assert outfile is None

    # Table without niriss_ghost_stamp column - cases that should work
    tab = Table()
    tab['a'] = [1]
    source_types = ['point_source']
    expected = [DEFAULT_NIRISS_PTSRC_GHOST_FILE]

    for row, source, expectation in zip(tab, source_types, expected):
        outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)
        assert outfile == expectation

    # Table without niriss_ghost_stamp column - cases that should raise exceptions
    tab = Table()
    tab['a'] = [1]
    source_types = ['galaxy']

    for row, source in zip(tab, source_types):
        with pytest.raises(ValueError):
            outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)

    # Table without niriss_ghost_stamp column - cases where outfile should be None
    tab = Table()
    tab['a'] = [1]
    source_types = ['extended']

    for row, source in zip(tab, source_types):
        outfile = niriss_ghosts.determine_ghost_stamp_filename(row, source)
        assert outfile is None


def test_get_gap():
    filts = ['CLEAR']
    pupils = ['F090W']
    expected_xs = [1168.900]
    expected_ys = [937.200]
    expected_fracs = [1.100]

    if not os.path.isfile(NIRISS_GHOST_GAP_FILE):
        config_dir, ghost_file = os.path.split(NIRISS_GHOST_GAP_FILE)
        download_file(NIRISS_GHOST_GAP_URL, ghost_file, output_directory=config_dir, force=True)

    for filt, pup, ex_x, ex_y, ex_f in zip(filts, pupils, expected_xs, expected_ys, expected_fracs):
        x, y, f = niriss_ghosts.get_gap(filt, pup, NIRISS_GHOST_GAP_FILE)
        assert np.isclose(x, ex_x)
        assert np.isclose(y, ex_y)
        assert np.isclose(f, ex_f)

    # Invalid filter should return NaNs
    x, y, f = niriss_ghosts.get_gap('CLEAR', 'F999W', NIRISS_GHOST_GAP_FILE)
    assert not np.isfinite(x)
    assert not np.isfinite(y)
    assert not np.isfinite(f)


def test_get_ghost():
    x = 100
    y = 100
    flux = 500

    if not os.path.isfile(NIRISS_GHOST_GAP_FILE):
        config_dir, ghost_file = os.path.split(NIRISS_GHOST_GAP_FILE)
        download_file(NIRISS_GHOST_GAP_URL, ghost_file, output_directory=config_dir, force=True)

    filt = 'CLEAR'
    pupil = 'F090W'
    xghost, yghost, fluxghost = niriss_ghosts.get_ghost(x, y, flux, filt, pupil, NIRISS_GHOST_GAP_FILE, shift=0)

    assert np.isclose(xghost, 2237.8)
    assert np.isclose(yghost, 1774.4)
    assert np.isclose(fluxghost, 5.5)

    pupil = 'F999W'
    xghost, yghost, fluxghost = niriss_ghosts.get_ghost(x, y, flux, filt, pupil, NIRISS_GHOST_GAP_FILE, shift=0)
    assert not np.isfinite(xghost)
    assert not np.isfinite(yghost)
    assert not np.isfinite(fluxghost)
