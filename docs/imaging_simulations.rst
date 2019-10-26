.. _imaging_simulator:

Simulating Imaging data with Mirage
===================================

The simulation of imaging mode observations (with NIRCam, NIRISS, or FGS) with Mirage is done using the *imaging_simulator.py* module. As discussed on the :ref:`3 stages <stages>` page, Mirage is broken down into three main stages. To streamline the creation of simulated data *imaging_simulator.py* wraps around these three stages. A single call to *imaging_simulator.py* with an :ref:`yaml parameter file <example_yaml>` will create a single simulated exposure from one detector.

The imaging simulator has only two possible inputs: the yaml parameter file, and an option to input a dark current ramp of the proper format. Here is an example call:

::

    from mirage.imaging_simulator import ImgSim

    sim = ImgSim(paramfile='my_yaml_file.yaml')
    sim.create()


If you have a fits file containing a dark current exposure that is the proper format (linearized, with the correct readout pattern and array size) for the simulated data you are creating, you can provide this via the ``override_dark`` keyword parameter. This will cause the dark current preparation step to be skipped, which will save some time. In practice, the only way to have a fits file with the properly formatted dark current exposure will be from previous runs of the imaging simulator (or :ref:`dark prep <dark_prep>` step).

.. tip::

    For more examples of calling the imaging simulator, see the `Imaging Simulator example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb>`_, as well as the `Simulations from input mosaic image notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Simulated_data_from_mosaic_image.ipynb>`_.
