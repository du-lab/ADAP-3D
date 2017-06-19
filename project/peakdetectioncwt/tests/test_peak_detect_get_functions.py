
import base_test_peak_detection_on_eic as btpdoe
import unittest
import numpy

class TestPeakDetectGetFunctions(btpdoe.BaseTestPeakDetectionOnEIC):

    def test_get_signal_to_noise_values(self):
        peak_signal_to_noise_vals = self.peakdetector.getallSNVals()

        self.assertTrue(type(peak_signal_to_noise_vals) is list,
                        'S/N values are not in a list')
        if len(peak_signal_to_noise_vals)>0:
            self.assertTrue(type(peak_signal_to_noise_vals[0]) is float,
                            'S/N values are not floats')
        self.assertEqual(len(self.allPeakPositions),len(peak_signal_to_noise_vals),
                         'different number of S/N values than number of peak positions')

    def test_get_get_all_coef_over_area_vals(self):
        coef_over_area_vals = self.peakdetector.getallCoefOverAreaVals()

        self.assertTrue(type(coef_over_area_vals) is list)
        if len(coef_over_area_vals)>0:
            self.assertTrue((type(coef_over_area_vals[0]) is float)
                            or (type(coef_over_area_vals[0]) is numpy.float64))
        self.assertEqual(len(self.allPeakPositions),len(coef_over_area_vals))

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPeakDetectGetFunctions)
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
    #unittest.main()


