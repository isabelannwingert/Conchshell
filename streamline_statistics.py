#! /usr/bin/env python

import os
import numpy as np
import pandas as pd

def read_connmat(txtfile, delimiter=' '):
    return np.genfromtxt(txtfile, delimiter=delimiter)

def symmetrize_connmat(connmat):
    c = connmat.copy()
    np.fill_diagonal(c, 0)
    c = (c.T + c)/2.0 
    return c 

def density(connmat):
    numnodes = connmat.shape[0]
    denom = (numnodes-1)*numnodes / 2 
    numer = (connmat[np.triu_indices(numnodes,1)] > 0).sum()
    return numer/denom 
    
def num_conn_comp(connmat):
    # binarize the connectome 
    connmat = (symmetrize_connmat(connmat)>0)*1 
    # make a dictionary of direct connections from each edge 
    edges = np.arange(connmat.shape[0])
    graph = {}
    for e in edges:
        row = connmat[e]
        graph[e] = np.where(row)[0]
    visited = set()
    NCC = 0
    def dfs(visited, graph, node): 
        # function for depth-first search 
        if node not in visited:
            visited.add(node)
            for neighbor in graph[node]:
                dfs(visited, graph, neighbor)
    while len(set(edges).difference(visited)) > 0:
        # grab the first node from the list of unvisited 
        node = list(set(edges).difference(visited))[0]
        dfs(visited, graph, node)
        NCC+=1 
    return NCC 
    
def average_network_strength(connmat):
    return connmat[np.triu_indices(connmat.shape[0],1)].mean()
    
def average_self_edges_strength(connmat):
    return np.diagonal(connmat).mean()
    
def average_interhemispheric_strength(connmat):
    # assumes atlas has even number of nodes 
    # first half of which are left hemisphere
    numnodes = connmat.shape[0]
    if numnodes %2 != 0:
        raise ValueError('Only even number of nodes supported')
    halfn = int(numnodes / 2)
    connmat = symmetrize_connmat(connmat)
    submat = connmat[halfn:numnodes, 0:halfn]
    return np.mean(submat)
    
def average_intrahemispheric_strength(connmat):
    # assumes atlas has even number of nodes 
    # first half of which are left hemisphere
    numnodes = connmat.shape[0]
    if numnodes %2 != 0:
        raise ValueError('Only even number of nodes supported')
    halfn = int(numnodes / 2)
    connmat = symmetrize_connmat(connmat)
    submatL = connmat[0:halfn, 0:halfn]
    submatR = connmat[halfn:numnodes, halfn:numnodes]
    submatLut = submatL[np.triu_indices(halfn, 1)]
    submatRut = submatR[np.triu_indices(halfn, 1)]
    #print(submatLut.shape)
    return np.mean(np.concatenate((submatLut, submatRut)))

def nseeds(subjectConnectomePath, nseeds_file="nseeds.txt"):
    nseeds_path = os.path.join(subjectConnectomePath, nseeds_file)
    with open(nseeds_path, "r") as file:
        n_seeds = int(file.read().replace('\n', ''))
    return(n_seeds)

def tckstats_statistics(subjectConnectomePath, tckstats_file="tckstats.txt"):
    tckstats_path = os.path.join(subjectConnectomePath, tckstats_file)
    tckstats_df = pd.read_csv(tckstats_path,
                              skiprows=2, nrows=1, header=None,
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