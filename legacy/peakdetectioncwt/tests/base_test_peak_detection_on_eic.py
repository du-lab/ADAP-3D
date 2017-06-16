
import unittest
import pylab as pl
from peakdetectioncwt import peakdetector as pd

class BaseTestPeakDetectionOnEIC(unittest.TestCase):
    def setUp(self):
        f = open('testing_eic.txt', 'r')
        d = pl.genfromtxt(f)
        intensities = d[0, :]
        rt = d[1, :]

        self.peakdetector = pd.PeakDetector()

        # intensities and rt (inputs) must be numpy arrays
        self.peakdetector.setSignal(intensities)
        self.peakdetector.setX(rt)
        self.peakdetector.setWaveletPeakHighestSimilarity("off")
        self.peakdetector.setSignalToNoiseThreshold(7)
        self.allLeftBounds, self.allRightBounds, self.allPeakPositions = self.peakdetector.findPeaks()

    def tearDown(self):
        self.peakdetector = None