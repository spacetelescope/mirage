.. _logging:

Logging with Mirage
===================

Creation of log files
---------------------

Currently, logging in Mirage is a 2 step process, born out of the desire to have a single log file for each module called, regardless of whether that is a call to **catalog_seed_image.py** or a call to **imaging_simulator.py** (which in turn calls catalog_seed_image.py, dark_prep.py, and obs_generator.py).

While Mirage is running, it will write log messages to a file called **mirage_latest.log** in the working directory. Upon successful completion, this log file will be copied into a new location. The new location depends on exactly which Mirage modules are being called. In cases where Mirage is using a :ref:`yaml input file <example_yaml>` to create a simulated exposure (via imaging_simulator.py, catalog_seed_image.py, dark_prep.py, obs_generator.py, wfss_simulator.py, or grism_tso_simulator.py), the log file will be copied to a subdirectory called **mirage_logs** which will be located within the output directory specified in the yaml file. The name of the copied log file will also be updated to be the name of the yaml file with the suffix ".log". In practice, this means that in order to find the log files, simply go to the directory containing your newly created simulated data, and look in the **mirage_logs** subdirectory.

If you are calling the **yaml_generator.py** in order to create input yaml files, then the log file will be copied to a subdirectory called **mirage_logs** located in the directory containing the :ref:`APT xml and pointing files <from_apt>`. In this case, the log file name will match the name of the APT xml file, with the suffix "log".

Note that the initial log file, **mirage_latest.log** will remain in the working directory. In the event of a Mirage crash, **mirage_latest.log** will contain the traceback.

Customizing logging
-------------------

The default beehavior is for Mirage to print all messages with a level of INFO or above to the screen, and to print all messages with a level of DEBUG or above into the log file. Users can customize the logs created by Mirage using the log configuration file contained in the Mirage repository. This file is located in mirage/logging/logging_config.yaml and specifies which types of logs Mirage creates and the message levels that are printed to the logs. By changing these values, users can customize the output log files.

