.. _obs_generator:

Generate Observation
====================

The final stage of the simulator is the observation generator. Here, the seed image is combined with the dark current exposure to create the final simulated data. This process includes several steps: addition of instrumental effects to the seed image, translation of the seed image from a count rate image to a signal ramp, the addition of noise sources, and finally summing the simulated signals with the dark current signals.


Translate Seed Image into a RAPID Exposure
------------------------------------------

The first task performed in the observation generator is to translate the seed image, which is a count rate image, into a "seed exposure" composed of multiple frames with signals in counts. Specifically, the seed image is translated into one or more integrations using the RAPID, NISRAPID, or FGSRAPID `readout pattern <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Detector+Readout+Patterns>`_. The exposure is the collection of these integrations.

Add Cosmic Rays
---------------

Cosmic rays are then added to the seed exposure. Cosmic ray details are controlled through the entries in the `cosmicRay` section of the input yaml file. The user must provide a pointer to the location of the cosmic ray library. The cosmic ray libraries that are included in *Mirage's* :ref:`reference files <reference_files>` are based on the cosmic ray library created by Massimo Robberto, which is detailed in his `2009 technical report <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-001928.pdf>`_.

Three separate cosmic ray libraries are provided in the reference file collection. These libraries are accessed using the :ref:`cosmicRay:library <library>` yaml file entry. The three options are **SUNMIN**, **SUNMAX**, and **FLARE**. Each library has an associated cosmic ray rate, shown below. Users can modify these rates through the multiplicative factor specified in the `cosmicRay:scale` entry.

Default cosmic ray rates (cosmic rays / pixel / second):

+--------+---------------------------+
| SUNMIN | .. math:: 6.153 * 10^{-5} |
+--------+---------------------------+
| SUNMAX | .. math:: 1.6955 * 10^{-4}|
+--------+---------------------------+
| FLARE  | .. math:: 0.10546         |
+--------+---------------------------+

Add Poisson Noise
-----------------

`Mirage` then adds Poisson noise to the seed exposure. This is done one frame at a time within the integration, with the signal for each frame dependent upon the signal and poisson noise from the previous frame. Through the `simSignals:poissonseed` entry in the input yaml file, users can provide a seed to the random number generator used to create the Poisson noise, in order to have reproducible results.

Flat Field Effects
------------------

If the user provides an illumination flat or pixel flat in the `Reffiles:illumflat` or `Reffiles:pixelflat` yaml file entries, these files are read in and multiplied into the seed exposure.

IPC
---

Interpixel capacitance (IPC) effects are added to the seed exposure at this point. IPC is an electronic coupling between adjacent pixels that causes some of the signal generated in a pixel to be measured by a neighboring pixel (`Donlon 2016 <https://ui.adsabs.harvard.edu/#abs/2016SPIE.9915E..2ID/abstract>`_). To add IPC, we convolve the image with a user-specified IPC kernel. Typically the kernel is a 3x3 array. The kernel lists the relative fraction of signal that will be measured in the 3x3 grid of pixels around a central pixel that sees some input flux. For example, in the case where 0.5% of the signal spreads to each of the 4 pixels immediately adjacent to the central pixel, we would use the kernel shown below.

+--------+--------+--------+
| 0.0002 | 0.0050 | 0.0002 |
+--------+--------+--------+
| 0.0050 | 0.9792 | 0.0050 |
+--------+--------+--------+
| 0.0002 | 0.0050 | 0.0002 |
+--------+--------+--------+

This is the inverse of the kernel that can be used by the JWST calibration pipeline to remove IPC effects, like that shown below. Therefore, `Mirage` has an option to accept a pipeline-formatted kernel that will then be inverted before being used to add IPC effects. It is also possible to use a file with a separate 3x3 kernel for each pixel. In this case the total kernel will have a shape of 2048 x 2048 x 3 x 3.


+--------+--------+---------+
| -0.0002| -0.0052| -0.0002 |
+--------+--------+---------+
| -0.0052| 1.0214 | -0.0052 |
+--------+--------+---------+
| -0.0002| -0.0052| -0.0002 |
+--------+--------+---------+

The latest IPC kernels for all instruments and detectors are provided in the reference files assocaited with `Mirage`.


Add Seed Exposure to Dark Current Exposure
------------------------------------------

At this point, the seed exposure is added to the dark current exposure from the previous stage. At this point, both are in the RAPID, NISRAPID, or FGSRAPID readout pattern, depending on instrument. They are added frame by frame, resulting in a single observation exposure.

Crosstalk
---------

Crosstalk effects can also be added to the observation. Crosstalk in this case is defined as the unintentional mirroring of a fraction of the measured signal across amplifiers. Bright sources located in one quadrant of a detector can create dimmed copies of themselves in the other three quadrants. Crosstalk on the NIRCam detectors has been quantified by Armin Rest. Results are detailed in his `2015 technical report <http://www.stsci.edu/files/live/sites/www/files/home/jwst/documentation/technical-documents/_documents/JWST-STScI-004361.pdf>`_. Crosstalk is added to the observation using an input file containing crosstalk coefficients. These coefficients control the magnitude and sign of the crosstalk signal between each combination of amplifiers. Currently, only the NIRCam detectors have non-zero crosstalk coefficients in the collection of *Mirage* :ref:`reference files <reference_files>`.

(Optionally) Return from Linearized to Raw State
------------------------------------------------

At this point, the observation contains linearized, partially calibrated (superbias and reference pixel subtracted) signal. This can be directly saved to an output fits file in the format used by the JWST calibration pipeline. There is also an option to return the observation exposure to a completely uncalibrated ("raw") state, such that the full calibration pipeline can be run on the data. The advantage of this is that pipeline parameters and settings can be customized/optimized by the user.

In order to return to a raw state, the linearization step of the pipeline is run in reverse on the data. The previously subtracted superbias and reference pixel signal is then added back in to the exposure.

Save Observation
----------------

With the construction of the observation exposure complete, the exposure is then saved in the requested format.

.. hint::
    The output format of the data is controlled by the **Output:datatype** line in the input yaml file. Possible values for this parameter include "raw", "linear", or "raw,linear". In the latter case, both the linearized and raw versions of the observation will be saved.