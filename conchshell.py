#! /usr/bin/env python

import numpy as np

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
    print(submatLut.shape)
    return np.mean(np.concatenate((submatLut, submatRut)))
    