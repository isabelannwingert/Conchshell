#! /usr/bin/env python

import sys
import argparse, textwrap
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import streamline_statistics as slstats
import zscores
from plotter import plot_heatmap, plot_hist, plot_measure

DESCRIPTION = '''
    Conchshell pipeline for DTI quality control, including the examination of connectome heatmaps, histograms, and derived metrics.
'''

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter)
    
    p.add_argument('-s', '--subject', action='store', metavar='<txt>', dest='subjecttxt',
                    type=str, required=True, 
                    help='A txt file containing subjectIDs that will undergo examination for quality control.'
                    )
    p.add_argument('-p', '--path', action='store', metavar='<txt>', dest='pathtxt',
                    type=str, required=True, 
                    help=textwrap.dedent('''\
                                        A comma-delimited txt file containing:
                                        column 1 (required): atlas name,
                                        column 2 (required): path template to the connectome file with {s} token for subjectIDs,
                                        column 3 (optional): path template to the nifti atlas in subject space with {s} token for subjectIDs,
                                        column 4 (optional): path to the atlas node order txt file.
                                        Example txt file:
                                        desikan,/path/to/{s}_desikan_connmat.txt,/path/to/{s}_desikan.nii.gz,desikan_lut.txt
                                        Schaefer100,/path/to/{s}_Schaefer100_connmat.txt,/path/to/{s}_Schaefer100.nii.gz,schaefer_node_order.txt
                                        ''')
                    )
    p.add_argument('-o', '--output', action='store', metavar='<dir>', dest='outputdir',
                    type=str, required=True, 
                    help='Path to the folder for the outputs.'
                    )
    p.add_argument('-t', '--tck', action='store', metavar='<str>', dest='tckstr',
                    type=str, required=False,
                    help=textwrap.dedent('''\
                                        The path template to the folders that store track statistics data (tckstats.txt, nseeds.txt).
                                        If not specified, the pipeline will by default search the folder where each connectome matrix is located. If not found, skip those measures.
                                        ''')
                    )
    p.add_argument('-c', '--covars', action='store', metavar='<csv>', dest='covarscsv',
                    type=str, required=False, default=None,
                    help='A csv file containing covariate features (e.g. Age, Group, Site) with subjectID in the leftmost column'
                    )
    # p.add_argument('-C', '--columns', action='store', metavar='<txt>', dest='columnstxt',
                    # type=str, required=False, default=None,
                    # help='A txt file containing the column names of covariates to adjust for during the quality control process'
                    # )
    p.add_argument('-f', '--formula', action='store', metavar='<str>', dest='formula',
                    type=str, required=False, default=None,
                    help='A patsy-like formula that defines the covariates to adjust for. Default is no covariates. Example: Age+Sex+Group+Site'
                    ) 
    return p

def main(argv):
    parser = buildArgsParser()
    args = parser.parse_args(argv)
    # Parse argements
    subjlist = [ a.strip() for a in open(args.subjecttxt, 'r').readlines()]
    atlaslist, connpaths, *params = np.genfromtxt(args.pathtxt, dtype='str', delimiter=',', unpack=True)
    if isinstance(atlaslist, str):  
        # if args.pathtxt has only 1 line genfromtxt returns a string 
        atlaslist = [atlaslist]
        connpaths = [connpaths]
        params = [[a] for a in params]
    if not os.path.isdir(args.outputdir):
        os.makedirs(args.outputdir)
    covars = None
    if args.covarscsv is not None:
        covars = pd.read_csv(args.covarscsv, index_col=0, na_values=[' ','na','nan','NaN','NAN','NA','#N/A','.','NULL'])
    # columns = None
    # if args.columnstxt is not None:
        # columns = [r.strip() for r in open(args.columnstxt,'r').readlines()]

    FIG_ROW = 4
    FIG_COL = 5
    outlier_threshold = 1
    # compute QC measures for each atlas
    for t, atlas in enumerate(atlaslist):
        pp_heatmap = PdfPages(os.path.join(args.outputdir, atlas+'_heatmaps.pdf'))
        pp_histogram = PdfPages(os.path.join(args.outputdir, atlas+'_histograms.pdf'))

        measures = {'Subject':[], 'file': [],
                    'density':[], 'num_connected_components':[], 'avg_network_strength':[],
                    'avg_interhemispheric_strength':[], 'avg_intrahemispheric_strength':[]
                    }

        fig1 = plt.figure(figsize=(25, 20))
        fig2 = plt.figure(figsize=(25, 20))
        max_edges = [slstats.read_connmat(connpaths[t].format(s=subject)).max() for subject in subjlist] 
        max_edges = np.asarray(max_edges)
        vmax = np.log(np.percentile(max_edges, 99)+1)
        
        for i, subject in enumerate(subjlist):
            print('Checking atlas:', atlas, 'connectome:', i+1, 'out of', len(subjlist))
            connmat = slstats.read_connmat(connpaths[t].format(s=subject))
            measures['Subject'].append(subject)
            measures['file'].append(connpaths[t].format(s=subject))
            measures['density'].append(slstats.density(connmat))
            measures['num_connected_components'].append(slstats.num_conn_comp(connmat))
            measures['avg_network_strength'].append(slstats.average_network_strength(connmat))
            measures['avg_interhemispheric_strength'].append(slstats.average_interhemispheric_strength(connmat))
            measures['avg_intrahemispheric_strength'].append(slstats.average_intrahemispheric_strength(connmat))

            # Plot heatmaps and histograms to pdf file
            findex = i % (FIG_COL*FIG_ROW)
            if  (i > 0) and (findex == 0):
                pp_heatmap.savefig(fig1)
                pp_histogram.savefig(fig2)
                fig1.clf()
                fig2.clf()
            ax = fig1.add_subplot(FIG_ROW, FIG_COL, findex+1)
            plot_heatmap(np.log(connmat+1), ax=ax, title=subject, vmax=vmax)
            ax = fig2.add_subplot(FIG_ROW, FIG_COL, findex+1)
            plot_hist(connmat, ax=ax, title=subject, bins=50, log=True)

        fig1.tight_layout()
        pp_heatmap.savefig(fig1)
        fig2.tight_layout()
        pp_histogram.savefig(fig2)
        pp_heatmap.close()
        pp_histogram.close()

        # Save QC measures to csv file
        measures_df = pd.DataFrame(measures)
        measures_df.to_csv(os.path.join(args.outputdir, atlas+'_measures_QC.csv'))
        #Adjust for covariates and compute z-scores
        zscores_df = zscores.corrected_zscores(measures_df.set_index('Subject'), covars=covars, formula=args.formula)  # columns=columns, 
        zscores_df.to_csv(os.path.join(args.outputdir, atlas+'_measures_zscore.csv'))

        # Plot subject-wise variance and detect outliers
        pp_zscore = PdfPages(os.path.join(args.outputdir, atlas+'_measures_zscore.pdf'))
        fig3 = plt.figure(figsize=(40, 5))
        for i, col in enumerate(zscores_df.columns):
            findex = i % FIG_COL
            if (i > 0) and (findex == 0):
                pp_zscore.savefig(fig3)
                fig3.clf()
            ax = fig3.add_subplot(1, FIG_COL, findex+1)
            plot_measure(zscores_df[col], ax=ax, title=col, violin=True, outlier_threshold=outlier_threshold, subject=zscores_df.index)
        fig3.tight_layout()
        pp_zscore.savefig(fig3)
        pp_zscore.close()

    #compute QC measures for tckstats.txt and nseeds.txt
    if args.tckstr is not None:
        tckpaths = args.tckstr
    else:
        tckpaths = os.path.dirname(connpaths[0])
    measures = {'Subject':[], 'count':[], 'nseeds':[], 'ratioCN':[],
                'mean_streamlength':[], 'median_streamlength':[], 'stdev_streamlength':[],
                'min_streamlength':[], 'max_streamlength':[]
                }
    tck_found = True
    for subject in subjlist:
        try:
            tckstats_df = slstats.tckstats_statistics(tckpaths.format(s=subject))
        except:
            tck_found = False
            break
        measures['Subject'].append(subject)
        measures['count'].append(slstats.count(tckstats_df))
        measures['nseeds'].append(slstats.nseeds(tckpaths.format(s=subject)))
        measures['ratioCN'].append(measures['count'][-1] / measures['nseeds'][-1])
        measures['mean_streamlength'].append(slstats.mean_streamlength(tckstats_df))
        measures['median_streamlength'].append(slstats.median_streamlength(tckstats_df))
        measures['stdev_streamlength'].append(slstats.stdev_streamlength(tckstats_df))
        measures['min_streamlength'].append(slstats.min_streamlength(tckstats_df))
        measures['max_streamlength'].append(slstats.max_streamlength(tckstats_df))
    
    if tck_found:
        measures_df = pd.DataFrame(measures)
        measures_df.to_csv(os.path.join(args.outputdir, 'tckstats_QC.csv'))
        #Adjust for covariates and compute z-scores
        zscores_df = zscores.corrected_zscores(measures_df.set_index('Subject'), covars=covars, formula=args.formula) # columns=columns, 
        zscores_df.to_csv(os.path.join(args.outputdir, 'tckstats_zscore.csv'))

        # Plot subject-wise variance and detect outliers
        pp_zscore = PdfPages(os.path.join(args.outputdir, atlas+'tckstats_zscore.pdf'))
        fig3 = plt.figure(figsize=(40, 5))
        for i, col in enumerate(zscores_df.columns):
            findex = i % FIG_COL
            if (i > 0) and (findex == 0):
                pp_zscore.savefig(fig3)
                fig3.clf()
            ax = fig3.add_subplot(1, FIG_COL, findex+1)
            plot_measure(zscores_df[col], ax=ax, title=col, violin=True, outlier_threshold=outlier_threshold, subject=zscores_df.index)
        fig3.tight_layout()
        pp_zscore.savefig(fig3)
        pp_zscore.close()
    print('Done!')

if __name__ == '__main__':
    main(sys.argv[1:])
