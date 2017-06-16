# About PD3D

PeakDetect3D is an open-source python package developed to detects peaks in liquid or gas 
chromatography mass-spectrometry data. It uses profile data and the 3D nature of peaks
to reduce false positives and requires very few parameters specified by the user. Many
parameters are determined by the algorithm from the data itself making the program 
versatile as well as easy to use.

# Sub-packages: 
* peakdetectioncwt
* datamodel
* easyIOmassspec
* generalcurvetools
* writing
* plotting

## Descriptions 

### peakdetectioncwt
This is a peak detection program for finding features/peaks in LCMS or GCMS EICs though it could
also be used for feature/peak detection in many other contexts as well. It uses a continuous wavelet
transform method to determine peaks and their boundaries. False positives can be reduced with several
filtering options. For more information see README in peakdetectioncwt folder.

### datamodel
A couple of classes to help with the organization of data including the hard coded parameters.
#### datamodel contains:
* data_point.py
* result.py
* peak.py
* hard_parameter.py

### easyIOmassspec
A way to easily read and write different mass spec data file types.  Currently this is under
development and only allows for the reading of mzXML and CDF files and the writing of CDF files.
Mostly a wrapper for the mspy library (http://www.mmass.org/) as far as mzXML files are concerned.

### generalcurvetools
* curve_tools.py: Contains functions for a very particular type of similarity measure, an asymmetric
    Gaussian fit and a function that estimates the full width half max of the mass spectrum peaks.
* curvs.py: A place to keep useful functions for possibly fitting to peaks.

### writing
Separated the details of writing the results to file into these modules

### plotting
Everything for plotting final results of the exact piece of data that are considered peaks.

### misc
Miscellaneous scripts and information for post processing data or looking at some statistics of the data set.

### jython_turn_results_into_mzmine2_xml.py
Script for taking results and turning them into peak list xml file that can be read by
MZmine 2.

## Example

`python main.py -f <path to file>
        --absoluteintensitythresh 500
        --peakintensitythresh 5000
        --numinitpeaks 20
        -o <path of location to put results>
        -n <name of folder created to store results> `

`--absoluteintensitythresh` discard all data points with intensities less than this value.
`--peakintensitythresh` is a parameter specified by the user. Only peaks with intensities higher
than this number will be considered real.
`--numinitpeaks` sets the number of most intense peaks used in the determination of parameters.



