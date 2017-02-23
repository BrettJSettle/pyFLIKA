# NOTE
This program is no longer developed and has been superceeded by the more general image processing program Flika (https://github.com/flika-org/flika) and its detect puff plugin (https://github.com/kyleellefsen/detect_puffs).

# pyFLIKA
An algorithm for automated detection, localization and measurement of local calcium signals from camera-based imaging


INTRODUCTION
------------
The pyFLIKA program is used to quickly locate and analyze local Calcium ion events, or puffs.  The algorithm is explained in detail in the JOVE paper here: (http://www.scitru.com/papers/2014_an_algorithm_for_automated_detection_localization_and_measurement.pdf).  Video recordings of cell activity are analyzed by light intensity to locate areas of illumination, normalized by spatial and temporal filtering, followed by analysis of the 'puffing' areas.

REQUIREMENTS
------------
Reliant on Python 3.4 or 2.7 (untested), preferably installed through Anaconda Scientific Library (http://continuum.io/downloads#py34) as well as several other python plugins:

Included in Anaconda:
*	PyQt 		v4	(will not install with 'pip')
*	numpy 	v1.9.2
*	scipy 	v0.16
*	scikit-learn	v0.16.1
*	xlrd		v0.9.3

Separate installs:
*	pyqtgraph	v0.9.10
*	biodocks	v0.7

To install plugins, open a command prompt and type 'pip install plugin-name', using the names as they appear above.

INSTALLATION
------------
To install the PyFLIKA Program, install all of the plugins listed above into your python environment, download the pyFLIKA from github (https://github.com/BrettJSettle/pyFLIKA) by clicking on the 'Download Zip' button on the right sidebar. Extract it to your computer and run by clicking the 'run.bat' file located inside the pyFLIKA Folder


USING THE PROGRAM
-----------------
IN PROGRESS

ABOUT
-----
*	Program: pyFLIKA
*	Author: Brett Settle
*	Date: August 13, 2015
*	Lab: UCI Parker Lab, McGaugh Hall
