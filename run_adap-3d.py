#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:13:55 2017

@author: xdu4
"""




import os
import subprocess



code_dir = "/Users/xdu4/Documents/Duxiuxia/Bitbucket/adap-3d"
os.chdir(code_dir)


raw_file_dir = "/Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw"
raw_file_name = "DC_010814_StandardsMixTest1_34StandardMix_01.mzXML"
raw_file_full_name = os.path.join(raw_file_dir, raw_file_name)

absolute_intensity_thresh = 500
peak_intensity_thresh = 5000

num_initial_peaks = 20

path_to_results_dir = "/Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/"
results_dir = "results"

subprocess.call(["python", "main.py", \
                 "-f", raw_file_full_name,\
                 "--absoluteintensitythresh", str(absolute_intensity_thresh),\
                 "--peakintensitythresh", str(peak_intensity_thresh),\
                 "--numinitpeaks", str(num_initial_peaks),\
                 "-o", path_to_results_dir,\
                 "-n", results_dir])

subprocess.call(["python", "/Users/xdu4/Documents/Duxiuxia/Bitbucket/adap-3d/main.py", "-f", raw_file_full_name, "--absoluteintensitythresh", str(absolute_intensity_thresh), "--peakintensitythresh", str(peak_intensity_thresh), "--numinitpeaks", str(num_initial_peaks), "-o", path_to_results_dir, "-n", results_dir])



#python main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw/DC_010814_StandardsMixTest1_34StandardMix_01.mzXML --absoluteintensitythresh 500 --peakintensitythresh 5000 --numinitpeaks 20 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d -n results

# jython

%run -d main.py -f /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d/raw/DC_010814_StandardsMixTest1_34StandardMix_01.mzXML --absoluteintensitythresh 500 --peakintensitythresh 5000 --numinitpeaks 20 -o /Users/xdu4/Documents/Duxiuxia/Analysis/my_projects/adap-3d -n results