#! /usr/bin/env python
"""
Created on Weds Sep 13 12:03:22 2023
@author: Isabel Ann Wingert, DiCIPHR Lab 
CONTACT: isabel.wingert@pennmedicine.upenn.edu
"""

import os
import numpy
import pandas as pd

# =============================================================================
# DESCRIPTION = '''
#    Compiles streamlength statistics.
#'''
# PROTOCOL_NAME='Conchshell'
#
# Will make changes to argsparser and use this to build onto for later, just drafting functions 
#   (down below) to return metrics of interest (9/13/2023)

# def builArgsParser():
#     p = argparse.ArgumentParser(description=DESCRIPTION)
#     p.add_argument('-s', '--subjects', action='store', 
#                    metavar='<csv>', dest='subjectfile',
#                    type=str, required=True,
#                    help='A text file containing subject IDs'
#                    )
#     p.add_argument('-f', '--filename', action='store', 
#                    metavar='<str>', dest='filelame_template',
#                    type=str, required=True,
#                    help="A template for the output statistics / Z-Scores, with {s}" 
#                    "to be replaced by the subject ID"
#                    )
#     p.add_argument('-o', '--outdir', action='store',
#                    metavar="<str>", dest='outdir',
#                    type=str, required=False,
#                    help="The output directory. Defaults to {}'.format"(os.path.join('$PWD', 'Protocols', PROTOCOL_NAME))
#                    )
# =============================================================================

def nseeds(subjectConnectomePath):
    nseeds_file = "nseeds.txt"
    nseeds_path = "\\".join([subjectConnectomePath, nseeds_file])
    with open(nseeds_path, "r") as file:
        n_seeds = int(file.read().replace('\n', ''))
    return(n_seeds)

def tckstats_statistics(subjectConnectomePath):
    tckstats_file = "tckstats.txt"
    tckstats_path = "\\".join([subjectConnectomePath, tckstats_file])
    with open(tckstats_path, "r") as file:
        tckstats_df = pd.read_csv("tckstats.txt", 
                                  skiprows=2, nrows=1, 
                                  header=None, 
                                  names=["mean", "median", "stdev", "min", "max", "count"], 
                                  skipinitialspace=True, delim_whitespace=True, index_col=False)
        
    # stores all metrics into dataframe, if user wants specific metric,
    # functions are below and named accordingly
    return(tckstats_df)

def mean_streamlength(tckstats_df):
    mean_streamlength = float(tckstats_df["mean"])
    # mean length of connections
    return(mean_streamlength)
    
def median_streamlength(tckstats_df):
    median_streamlength = float(tckstats_df["median"])
    # median length of connections
    return(median_streamlength)

def stdev_streamlength(tckstats_df):
    stdev_streamlength = float(tckstats_df["stdev"])
    # standarad deviation of connection lengths
    return(stdev_streamlength)

def min_streamlength(tckstats_df):
    min_streamlength = float(tckstats_df["min"])
    # minimum length of connections measured
    return(min_streamlength)

def max_streamlength(tckstats_df):
    max_streamlength = float(tckstats_df["max"])
    # maximum length of connections measured
    return(max_streamlength)

def count(tckstats_df):
    count = float(tckstats_df["count"])
    # represents the number of streamlines in the tractogram
    return(count)
    
def ratioCN(n_seeds, count):
    ratioCN = count / n_seeds
    return(ratioCN)

# Next steps: What else do we need to do to organize after we gather all of the metrics above?

