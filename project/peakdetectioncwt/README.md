

# About peakdetectioncwt
This is a peak detection program for finding 
featuresi/peaks in LCMS or GCMS EICs though it 
could also be used for feature/peak detection 
in many other contexts as well.

# Sub-packages:
* cwt 
* tests


## Descriptions

### peakdetector.py 
Contains the class for the user to interface with. It also contains 
some additional simple peak filtering steps like minimum height 
chech and a CWT coefficient over peak area check. 

### cwt
The continuous wavelet transform and ridgeline detection
can be found in here. The signal to noise estimation method
is also in here. Though not directly related to CWT, signal to 
noise estimation is an important part of peak detection.


### tests
Contains a couple of small and non-comprehensive tests.



## Example

What is shown and explained below can all
be found in the example_use_peakdetector.py script.

`f = open("TestEICs/testing_eic2.txt","r")` Open a test file containing a single EIC.

`d = pl.genfromtxt(f)` 

`intensities = d[0,:]` Parse out the intensities of the EIC,

`rt = d[1,:]` and the retention times.

`peakdetector = pd.PeakDetector()` make the interface object

`peakdetector.setX(rt)`

`peakdetector.setSignal(intensities)` Pass in the intensities and retention times. They must both be numpy arrays.

`peakdetector.setMinIntensityThreshold(100)` Peaks with intensities below this value will not be detected.

`peakdetector.setCoefOverAreaThreshold(10)` Peaks with C/A below this value will not be detected.

`peakdetector.setWaveletPeakHighestSimilarity("off")` Some thresholds can be turned off and completely disregarded.

`peakdetector.setSignalToNoiseThreshold(5)` Peaks with S/N below this value will not be detected.

`peakdetector.setVisualize(False)` Turns some very simple visualization plots off. If this is on they will pop up
through the process showing useful information.

`allLeftBounds, allRightBounds , allPeakPositions = peakdetector.findPeaks()` Do 
the peak detection and get the results! `allLeftBounds` (`allRightBounds`) is a simple list
of the indices in the `rt` and `signal` arrays of all the detected left (right) boundaries. 
`allPeakPositions` are the indices of the highest intensity point of each detected peak.

`allSNValues = peakdetector.getallSNVals()` More specific information can be gathered as well.
This returns the S/N values for each peak.

`peakdetector.showSNAnnotatedEIC()` An annotated plot of the signal with 
each of the S/N values of each peak will pop up.

