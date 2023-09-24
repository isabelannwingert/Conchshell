#! /usr/bin/env python

import sys
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import streamline_statistics as slstats
from plotter import plot_heatmap, plot_hist

DESCRIPTION = '''
    Conchshell pipeline for DTI quality control, including the examination of connectome heatmaps, histograms, and derived metrics.
'''

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    
    p.add_argument('-i', '--input', action='store', metavar='<txt>', dest='inputtxt',
                    type=str, required=True, 
                    help='A two-column comma-delimited txt file containing subjectIDs and paths to each connectome matrix that will undergo examination for quality control'
                    )
    p.add_argument('-o', '--output', action='store', metavar='<dir>', dest='outputdir',
                    type=str, required=True, 
                    help='Path to the folder for the outputs'
                    )
    p.add_argument('-t', '--tck', action='store', metavar='<txt>', dest='tcktxt',
                   type=str, required=False,
                   help='The txt file containing a list of paths to folders that store track statistics data (tckstats.txt, nseeds.txt). If not specified, the pipeline will by default search the folder where each connectome matrix is located',
                   )
    p.add_argument('-c', '--covars', action='store', metavar='<csv>', dest='covarscsv',
                    type=str, required=False, default=None,
                    help='A csv file containing covariate features (e.g. Age, Group, Site) with subjectID in the leftmost column'
                    )
    p.add_argument('-C', '--columns', action='store', metavar='<txt>', dest='columnstxt',
                    type=str, required=False, default=None,
                    help='A txt file containing the column names of covariates to adjust for during the quality control process'
                    )
    p.add_argument('-f', '--formula', action='store', metavar='<str>', dest='formula',
                    type=str, required=False, default=None,
                    help='A patsy-like formula that defines the covariates to adjust for. Default is no covariates. Example: Age+Sex+Group+Site'
                    ) 
    return p

def main(argv):
    parser = buildArgsParser()
    args = parser.parse_args(argv)

    # Import connectome data and compute QC measures
    subjlist, connpaths = np.genfromtxt(args.inputtxt, dtype='str', delimiter=',', unpack=True)
    if args.tcktxt is not None:
        tckpaths = [r.strip() for r in open(args.tcktxt,'r').readlines()]
    else:
        tckpaths = [os.path.dirname(i) for i in connpaths]
    if not os.path.isdir(args.outputdir):
        os.makedirs(args.outputdir)
    pp_heatmap = PdfPages(os.path.join(args.outputdir, 'heatmaps.pdf'))
    pp_histogram = PdfPages(os.path.join(args.outputdir, 'histograms.pdf'))

    measures = {'Subject':[], 'file': [],
                'density':[], 'num_connected_components':[], 'avg_network_strength':[],
                'avg_interhemispheric_strength':[], 'avg_intrahemispheric_strength':[],
                'count':[], 'nseeds':[], 'ratioCN':[],
                'mean_streamlength':[], 'median_streamlength':[], 'stdev_streamlength':[],
                'min_streamlength':[], 'max_streamlength':[]
                }

    FIG_ROW = 4
    FIG_COL = 5
    fig1 = plt.figure(figsize=(25, 20))
    fig2 = plt.figure(figsize=(25, 20))
    for i, subject in enumerate(subjlist):
        connmat = slstats.read_connmat(connpaths[i])
        measures['Subject'].append(subject)
        measures['file'].append(connpaths[i])
        measures['density'].append(slstats.density(connmat))
        measures['num_connected_components'].append(slstats.num_conn_comp(connmat))
        measures['avg_network_strength'].append(slstats.average_network_strength(connmat))
        measures['avg_interhemispheric_strength'].append(slstats.average_interhemispheric_strength(connmat))
        measures['avg_intrahemispheric_strength'].append(slstats.average_intrahemispheric_strength(connmat))

        tckstats_df = slstats.tckstats_statistics(tckpaths[i])
        measures['count'].append(slstats.count(tckstats_df))
        measures['nseeds'].append(slstats.nseeds(tckpaths[i]))
        measures['ratioCN'].append(slstats.ratioCN(measures['count'][-1] / measures['nseeds'][-1]))
        measures['mean_streamlength'].append(slstats.mean_streamlength(tckstats_df))
        measures['median_streamlength'].append(slstats.median_streamlength(tckstats_df))
        measures['stdev_streamlength'].append(slstats.stdev_streamlength(tckstats_df))
        measures['min_streamlength'].append(slstats.min_streamlength(tckstats_df))
        measures['max_streamlength'].append(slstats.max_streamlength(tckstats_df))

        # Plot heatmaps and histograms to pdf file
        findex = i%(FIG_COL*FIG_ROW)
        if  (i > 0) and (findex == 0):
            pp_heatmap.savefig(fig1)
            pp_histogram.savefig(fig2)
            fig1.clf()
            fig2.clf()
        ax = fig1.add_subplot(FIG_ROW, FIG_COL, findex+1)
        plot_heatmap(connmat, ax=ax, title=subject)
        ax = fig2.add_subplot(FIG_ROW, FIG_COL, findex+1)
        plot_hist(connmat, ax=ax, title=subject)

    pp_heatmap.savefig(fig1)
    pp_histogram.savefig(fig2)
    pp_heatmap.close()
    pp_histogram.close()

    # Save QC measures to csv file
    measures_df = pd.DataFrame(measures)
    measures_df.to_csv(os.path.join(args.outputdir, 'QC_measures.csv'))
    
    
    #TO DO: adjust for covariates and compute z-scores
    covars = None
    if args.covarscsv is not None:
        covars = pd.read_csv(args.covarscsv, index_col=0, na_values=[' ','na','nan','NaN','NAN','NA','#N/A','.','NULL'])
    columns = None
    if args.columnstxt is not None:
        columns = [r.strip() for r in open(args.columnstxt,'r').readlines()]
    
    #TO DO: plot subject-wise variance and detect outliers

if __name__ == '__main__':
    main(sys.argv[1:])
