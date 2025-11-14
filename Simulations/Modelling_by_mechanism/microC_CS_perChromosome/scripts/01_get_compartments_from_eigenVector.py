#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import seaborn as sns

# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.1'):
    raise AssertionError("tutorials rely on cooltools version 0.5.1 or higher,"+
                         "please check your cooltools version and update to the latest")

import sys
import argparse
import bioframe

# Input files from Flora

for minE1 in [0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]:
    for minLength in [0, 5]:
        for condition in ["DMSO","TSA","TSA24hREC"]:
            
            print("Condition ", condition, "minE1 ", minE1, "minLength ", minLength)
            
            eigenvector_track_cis_eigs_100kb = bioframe.read_table("/zssd/scratch/mdistefano/2023_10_22_Project_hyperacetylation/02_AB_compartments/eigenvector_track_cis_eigs_100kb_clr_s%s.bedgraph" % condition, header=0)

            try:
                # Extract and merge A compartments
                A_compartments_cis_eigs_100kb = bioframe.merge(eigenvector_track_cis_eigs_100kb[eigenvector_track_cis_eigs_100kb["E1"] > +minE1], min_dist=0)    
                # Select compartments at least 500kb in size
                A_compartments_cis_eigs_100kb = A_compartments_cis_eigs_100kb[A_compartments_cis_eigs_100kb["n_intervals"] >= minLength]
                A_compartments_cis_eigs_100kb[["chrom", "start", "end"]].to_csv("A_compartments_cis_s%s_minLength_%dkb_minE1_%.2f.bed" % (condition,minLength*100,minE1), sep="\t", header=None, index=False)
                
                # Extract and merge B compartments
                B_compartments_cis_eigs_100kb = bioframe.merge(eigenvector_track_cis_eigs_100kb[eigenvector_track_cis_eigs_100kb["E1"] < -minE1], min_dist=0)
                # Select compartments at least 500kb in size
                B_compartments_cis_eigs_100kb = B_compartments_cis_eigs_100kb[B_compartments_cis_eigs_100kb["n_intervals"] >= minLength]
                B_compartments_cis_eigs_100kb[["chrom", "start", "end"]].to_csv("B_compartments_cis_s%s_minLength_%dkb_minE1_%.2f.bed" % (condition,minLength*100,minE1), sep="\t", header=None, index=False)
            except:
                pass
                
            print('finished script and saved compartment tracks')
