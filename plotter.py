#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def plot_heatmap(connmat, ax=None, title=None, row_labels=None, col_labels=None, colorbar=False, annot=False, **kwargs):
    if ax is None:
        ax = plt.gca()
    im = ax.imshow(connmat, **kwargs)
    #ax.set_xticks(np.arange(connmat.shape[0]))
    #ax.set_yticks(np.arange(connmat.shape[1]))
    if title is not None:
        ax.set_title(title)
    if col_labels is not None:
        ax.set_xticklabels(col_labels)
        plt.setp(ax.get_xticklabels(), rotation=30, ha="center", rotation_mode="anchor")
    if row_labels is not None:
        ax.set_yticklabels(row_labels)
    if colorbar:
        ax.figure.colorbar(im, ax=ax)
    if annot:
        for i in range(connmat.shape[0]):
            for j in range(connmat.shape[1]):
                ax.text(j, i, int(connmat[i, j]), ha="center", va="center", color="w")
    return im


def plot_hist(connmat, ax=None, title=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    im = ax.hist(connmat[np.triu_indices(connmat.shape[0],1)], **kwargs)
    ax.set_xticklabels(ax.get_xticks(), rotation=30)
    if title is not None:
        ax.set_title(title)
    return im
