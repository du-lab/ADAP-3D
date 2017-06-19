Coppyright 2016, Owen Myers, All rights reserved
Author Owen Myers
------------------------------------------------

Name of project: easyIOmassspec
Sub-moduals: mspy, Tests


**********************************
********* Descriptions ***********
**********************************
easyIOmassspec:
A way to easily read and write different mass spec data file types. 
Currently this is under development and only alows for the reading of
mzXML files and the writing of CDF files.

mspy:
To read mzXML files we re-wrape a set of functions in the mMass (http://www.mmass.org/)
library, specificly the mspy folder is taken directly from the mMass project. We have propagated
their license and awknolege their hard work.

Tests:
simple scripts to test the easyIOmassspec modual.


**********************************
********** Examples **************
**********************************

For an example of writing to CDF see the test script in Test.


Example of reading mzXML:

    from easyIOmassspec import easyio as eio
    
    for i in range(dfread.get_num_scans()):
        cur_mz_vals,cur_int_vals = dfread.get_next_scan_mzvals_intensities()
    
        # do some stuff


**********************************
*********** LICENSE **************
**********************************
GNU GENERAL PUBLIC LICENSE

----->  See/read LICENSE file
