#! /usr/bin/env python

import sys, argparse
import numpy as np
import pandas as pd 
import statsmodels.formula.api as smf

DESCRIPTION = '''
    Calculates z-scores from data features including correction for covariates of interest.
'''

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    
    p.add_argument('-i', '--input', action='store', metavar='<csv>', dest='inputcsv',
                    type=str, required=True, 
                    help='The csv file containing sample features (e.g. roi stats) with subjectID in the leftmost column and, optionally, covariate columns'
                    )
    p.add_argument('-o', '--output', action='store', metavar='<csv>', dest='outputcsv',
                    type=str, required=True, 
                    help='The csv file to write the results'
                    )
    p.add_argument('-c', '--covars', action='store', metavar='<csv>', dest='covarscsv',
                    type=str, required=False, default=None,
                    help='A csv file containing covariate features (e.g. Age, Group, Site) with subjectID in the leftmost column, if input csv does not contain these columns'
                    )
    p.add_argument('-C', '--columns', action='store', metavar='<txt>', dest='columnstxt',
                    type=str, required=False, default=None,
                    help='A txt file of column names to adjust. If not provided, will correct all columns that are detected as containing numerical data'
                    )
    p.add_argument('-f', '--formula', action='store', metavar='<str>', dest='formula',
                    type=str, required=False, default=None,
                    help='A patsy-like formula that defines the covariates to correct for. Default is no covariates. Example: Age+Sex+Group+Site'
                    ) 
    return p

def corrected_zscores(dataframe, columns=None, covars=None, formula=None):
    if columns is None:
        columns = list(dataframe.select_dtypes(include=[np.number]).columns.values)
    if covars is not None:
        dataframe = dataframe.join(covars, how='inner')
    if formula is not None:
        # add a non-numeric character to column names because statsmodels cant interpret names like 001_SPG_L
        rcolumns = ['r'+c for c in columns]
        rename_dict = dict(zip(columns, rcolumns))
        rdataframe = dataframe.rename(rename_dict, axis='columns')
        resids = pd.DataFrame(0.0, columns=rcolumns, index=rdataframe.index)
        for col in rcolumns:
            model = smf.ols('{0} ~ {1}'.format(col, formula), data=rdataframe)
            res = model.fit()
            pred = res.predict(rdataframe)
            resids[col] = rdataframe[col] - pred 
        mean_ = np.array(np.mean(resids, axis=0))
        std_ = np.array(np.std(resids, axis=0))
        zscores = pd.DataFrame((resids.values - mean_[None,:])/std_[None,:], columns=rcolumns, index=dataframe.index)
    else:
        # no correction just z score 
        mean_ = np.array(np.mean(dataframe[columns], axis=0))
        std_ = np.array(np.std(dataframe[columns], axis=0))
        zscores = pd.DataFrame((dataframe[columns].values - mean_[None,:])/std_[None,:], columns=columns, index=dataframe.index)
        rcolumns = columns
    ret = dataframe[columns].copy()
    for c, rc in zip(columns, rcolumns):
        ret.loc[:,c] = zscores.loc[:,rc]
    return ret 

def main(argv):
    parser = buildArgsParser()
    args = parser.parse_args(argv)
    dataframe = pd.read_csv(args.inputcsv, index_col=0, na_values=[' ','na','nan','NaN','NAN','NA','#N/A','.','NULL'])
    covars = None
    if args.covarscsv is not None:
        covars = pd.read_csv(args.covarscsv, index_col=0, na_values=[' ','na','nan','NaN','NAN','NA','#N/A','.','NULL'])
    columns = None
    if args.columnstxt is not None:
        columns = [r.strip() for r in open(args.columnstxt,'r').readlines()]
    zscores_df = corrected_zscores(dataframe, columns=columns, covars=covars, formula=args.formula)
    zscores_df.to_csv(args.outputcsv)
    
if __name__ == '__main__':
    main(sys.argv[1:])
