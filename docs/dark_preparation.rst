.. _dark_prep:

Preparation of the Dark Current Exposure
========================================

After creating the seed image `mirage` moves on to preparing the input dark current exposure. This involves manipulating the input dark current exposure in several ways.

Adjust Number of Frames/Groups and Integrations
-----------------------------------------------

Frist the number of groups and integrations are adjusted to match the requested groups and integrations specified in the input yaml file. If the input dark has more groups or integrations than the requested output, the extras are removed. If additional groups or integrations are required above what is present in the dark file, then copies of existing groups/integrations are made and appended to the observation. In the case where extra groups are needed, copies of the existing groups are added to the signal in the final group of the existing dark, such that the dark current signals continue increasing in as they would in a longer integration.

Put into Requested Readout Pattern
----------------------------------

Next, the `readout pattern <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Detector+Readout+Patterns>`_ of the dark current exposure is adjusted to match the requested output. Effectively this only works in two cases. If the input dark and output exposure have the same readout pattern, then no action is taken. If the input dark is in the RAPID or NISRAPID (for NIRCam and NIRISS respectively) readout mode, where each group is composed of only one frame, then frames can be averaged/dropped in such a way to transform the dark into any other readout pattern. Input darks using readout patterns other than RAPID/NISRAPID cannot be translated into any other readout patterns, as all other readout patterns have groups composed of averaged frames, such that the original information from all frames is not accessible.

Calibration and Linearization
-----------------------------

If the input dark current exposure is "raw", meaning that no calibration steps have been performed on the data, then the reorganized dark is run through the superbias and reference pixel subtraction, as well as the non-linearity correction steps of the JWST calibration pipeline. The signal subtracted in the superbias and reference pixel subtraction steps is saved. Once the simulated data are created, this signal can be added back in order to create a "raw" exposure if requested.

There is also an option to provide a linearized, rather than raw, dark current exposure. In this case, the file containing the linearized dark current also must contain the subtracted superbias and reference pixel signal, and the pipeline steps are subsequently skipped. ]

Crop to Requested Subarray
--------------------------

Once the linearized dark current exposure is present, it is then optionally cut down to the requested `subarray <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Detector+Subarrays>`_ size. At this point, the dark current data are ready to have signals from simulated sources added.

.. tip::

    Both raw and linearized dark current exposures are provided in the collection of :ref:`reference files <reference_files>` that accompany `mirage`.