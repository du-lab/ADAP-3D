def define_variables():
    PLOT_ALL_PEAKS = True

    VERBOSE = False

    USE_HARD_CODED_DETECTION_PARAMETERS = False
    USE_ISOTOPE_PARAMETERS_FOR_ALL = False



    USE_SMALL_TEST_WINDOW = False







    # isotopes will be shown for the tmp_mz_of_peak value set
    ONLY_VISUALIZE_ISOTOPES_FOR_NORMAL_DETECTED_PEAK = False

    if USE_SMALL_TEST_WINDOW:
        tmp_mz_of_peak = 212.006378174  # 20.25 21.75 rt
        RT_MIN = 17.25
        RT_MAX = 21.75
        MZ_MIN = tmp_mz_of_peak - 10
        MZ_MAX = tmp_mz_of_peak + 10

    ##########################################################################
    ########### Important numbers that need to be set ########################
    ##########################################################################
    HP = hard_parameter.RequiredParameters()
    if USE_ISOTOPE_PARAMETERS_FOR_ALL:
        HP.use_isotope_parameters_for_all()

    ##########################################################################
    ############ Done setting important numbers ##############################
    ##########################################################################

    if ONLY_VISUALIZE_ISOTOPES_FOR_NORMAL_DETECTED_PEAK:
        investigate_specific_mz = [-1]
    else:
        investigate_specific_mz = [-1]
        # investigate_specific_mz = [tmp_mz_of_peak]

    RESULT_LIST = []
    ISOTOPE_RESULT_PEAK_LIST = []
