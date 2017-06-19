# Algorithm of ADAP-3D python code

The conventional peak detection algorithms consist of creating a centroid data,
chromatogram building, and peak detection. In our algorithm, we skip
centroiding and chromatogram building. Instead, we work with profile raw data
and detect 3D peaks that contain both retention time and m/z information.

In order to make the algorithm parameter-free,
we __estimate__ _peak detection parameters_ as follows:

1. Set up initial parameters and detect several (e.g. 20) highest-intensity
peaks, which are usually easy to detect.

2. Estimate parameters from the detected peaks and perform peak detection on
the rest of the data.

In both steps above, we perform the same __peak detection algorithm__:

1. Find the highest-intensity point in the whole data. Say, this point has
coordinates _(mz, rt, int)_, which are the m/z value, retention time, and
intensity respectively. If *int* is lower than a threshold, stop the algorithm.

2. Get a slice of the data over the interval *(rt - rt_span, rt + rt_span)* at
fixed m/z value *mz*. Use continuous wavelet transform (CWT) to detect peaks in
the data slice.

    * If no peaks are found in the data, we set intensities of all data
    points in the slice to 0, so that they would not participate in the
    following peak detection. Then, continue to Step 1.

    * If peaks are detected, but neither of their apexes coincides with
    the highest-intensity point *(mz, rt, int)*, we set intensities of all data
    points except the detected peaks to 0 and continue to Step 1.

    * If a peak *P* is found in the slice
    and its apex coincides with the point *(mz, rt, int)*, then estimate how good
    the peak *P* is by looking into the adjacent slices at m/z values
    *mz + delta, mz - delta, mz + 2 delta, mz - 2 delta* and so on, where
    *delta* is the distance between two adjacent data points in m/z-direction.

    The peak *P* is considered to be good if we can detect similar peaks in
    *N* adjacent slices. The peak similarity is determined by calculating the
    normalized area between the elution profiles of two peaks. The constant *N*
    is estimated at the beginning of the algorithm by finding the average FWHM
    of the highest m/z-profile peaks in several randomly chosen scans.

    We also fit the asymmetric gaussian to estimate the shape of peak *P*. If
    the asymmetric gaussian doesn't fit well, we consider peak *P* to be bad.

    If the peak *P* is good (i.e. we detected similar peaks in the adjacent
    slices and we have a good fit of the asymmetric gaussian), then we save
    the information about that peak.

    Regardless if peak *P* is good or bad, we set intensities of the data points
    over the interval *(rt_start, rt_end)* in the processed slices to 0. Here, *rt_start, rt_end* are the boundaries of peak *P*. Next, we continue to Step 1.
