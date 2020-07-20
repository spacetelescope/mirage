.. _imaging_simulator:

Simulating Imaging data with Mirage
===================================

The simulation of imaging mode observations (with NIRCam, NIRISS, or FGS) with Mirage is done using the *imaging_simulator.py* module. As discussed on the :ref:`3 stages <stages>` page, Mirage is broken down into three main stages. To streamline the creation of simulated data *imaging_simulator.py* wraps around these three stages. A single call to *imaging_simulator.py* with an :ref:`yaml parameter file <example_yaml>` will create a single simulated exposure from one detector.

The imaging simulator has only two possible inputs: the yaml parameter file, and an option to input a dark current ramp of the proper format. Here is an example call:

::

    from mirage.imaging_simulator import ImgSim

    sim = ImgSim(paramfile='my_yaml_file.yaml')
    sim.create()

.. _img_provide_segmented_darks:

If you have dark current products for your simulation from a previous run of Mirage, it is possible to provide these files to *imaging_simulator.py* as inputs using the **override_dark** parameter. In this case, the call to *dark_prep.py* will be skipped, which will save some computing time. Note that the dark current products must be specific to your exoposure, in that they must contain arrays of the proper shape. So in practice, this detail is useful if you are repeating a previous call to Mirage. In that case, the darks can be provided as shown below. In the case where a single dark file is needed, it can be provided as a string or a 1-element list. In cases where the exposure is broken into segments and there are multiple dark files, these files must be provided as a list.

::

    from mirage.imaging_simulator import ImgSim

    # Single dark file
    dark = 'jw09996001001_01101_00001_nrcb5_uncal_linear_dark_prep_object.fits'
    m = ImgSim(override_dark=dark)
    m.paramfile = 'jw09996001001_01101_00001_nrcb5.yaml'
    m.create()

    # Exposure broken into multiple segments
    dark = ['jw09996001001_01101_00001_nrcb5_uncal_seg001_linear_dark_prep_object.fits',
            'jw09996001001_01101_00001_nrcb5_uncal_seg002_linear_dark_prep_object.fits']
    m = ImgSim(override_dark=dark)
    m.paramfile = 'jw09996001001_01101_00001_nrcb5.yaml'
    m.create()

.. tip::

    For more examples of calling the imaging simulator, see the `Imaging Simulator example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb>`_, as well as the `Simulations from input mosaic image notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Simulated_data_from_mosaic_image.ipynb>`_.
