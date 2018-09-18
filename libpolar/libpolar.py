#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 18:00:53 2018

@author: antony
"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append('/ifs/scratch/cancer/Lab_RDF/abh2138/scripts/python/')
import libplot
import lib10x
import libcluster
import math
from scipy.interpolate import UnivariateSpline

BINS = 50
K = 5
TWO_PI = 2 * math.pi

def time2theta(time):
    x = time

    x = (x - x.min()) / (x.max() - x.min())
    return (x * 2 * math.pi)


def bin_time_data(data, time):
    time_idx = np.argsort(time)
    
    ts = time[time_idx]
    ds = data[time_idx]

    d_bins = np.empty((0, 1))
    t_bins = np.empty((0, 1))
    d_std_bins = np.empty((0, 1))
    
    idx = 0
    
    tcurrent = 0
    tinc = 2 * math.pi / BINS
    
    while idx < len(time_idx):
        #print(idx)
        tnext = tcurrent + tinc

        e = np.empty((0, 1))
        
        while ts[idx] < tnext:
            v = ds[idx]
            
            if v > 0:
                e = np.append(e, v)
            
            idx += 1
            
            if idx == len(time_idx):
                break
        
        if e.size > 0:
            d_bins = np.append(d_bins, np.mean(e))
            d_std_bins = np.append(d_std_bins, np.std(e))
        else:
            d_bins = np.append(d_bins, 0)
            d_std_bins = np.append(d_std_bins, 0)
        
        t_bins = np.append(t_bins, tcurrent)
        
        tcurrent = tnext
        
    return d_bins, d_std_bins, t_bins

def new_fig():
    fig = libplot.polar_fig(w=8, h=8)
    ax = libplot.polar_clock_ax(fig)
    
    return fig, ax

def set_rho_lim(rho, ax):
    ax.set_rlim(0, 10 * int(np.ceil(np.max(rho) / 10.0)))

def cluster_polar_plot(theta, rho, clusters, marker='o', s=libplot.MARKER_SIZE, edgecolors='none', label='', fig=None, ax=None):
    """
    Create a tsne plot without the formatting
    """
    
    if fig is None:
        fig, ax = new_fig()
    
    colors = libcluster.colors()
      
    ids = list(sorted(set(clusters['Cluster'])))
      
    for i in range(0, len(ids)):
        l = ids[i]
        
        indices = np.where(clusters['Cluster'] == l)[0]
        
        n = len(indices)
        
        label = 'C{} ({:,})'.format(l, n)
        
        t = theta[indices]
        r = rho[indices]
          
        ax.scatter(t, r, c=colors[i], edgecolors=edgecolors, s=s, marker=marker, alpha=libplot.ALPHA, label=label)

    set_rho_lim(rho, ax)

    libcluster.format_legend(ax, cols=6, markerscale=2)    
       
    return fig, ax



def exp_polar_plot(theta, rho, exp, marker='o', edgecolors='none', s=libplot.MARKER_SIZE, cmap=plt.cm.plasma, norm=None, fig=None, ax=None):
    """
    Create a polar tsne plot for gene exp
    """
    
    new_mode = False
    
    if ax is None:
        fig, ax = new_fig()
        new_mode = True
    
    ax.scatter(theta, rho, c=exp, cmap=cmap, edgecolors=edgecolors, s=s, marker=marker, alpha=libplot.ALPHA, norm=norm)

    if new_mode:
        set_rho_lim(rho, ax)

    return fig, ax


def plot_gene_exp_times(data, theta, rho, genes, name, cmap=plt.cm.plasma):
    if type(genes) is pd.core.frame.DataFrame:
        genes = genes['Genes'].values
        
    w, h, rows, cols = lib10x.expr_grid_size(genes)

    print(w, h, rows, cols)

    fig = libplot.polar_fig(w=w, h=h) # libplot.newfig(w=16, h=16, subplot=441)
    
    si = 1
    
    for gene in genes:
        if type(gene) == list:
            exp = np.zeros(data.shape[1])
            
            for g in gene:
                print('g ' + g)
                
                inv = g.endswith('-')
                
                if g.endswith('-') or g.endswith('+'):
                    g = g[0:(len(g) - 1)]
                    
                v = data.loc[data.index.str.endswith(g), :].T.iloc[:, 0].values
                
                if inv:
                    exp += (4 - v)
                else:
                    exp += v
            
            #exp /= len(gene)
        else:
            exp = data.loc[data.index.str.endswith(gene), :].T.iloc[:, 0].values
    
    
        # sort by exp so high exp always on top
        idx = np.argsort(exp)
        
        t = theta[idx]
        r = rho[idx]
        e = exp[idx]
        
        ax = libplot.polar_clock_ax(fig, subplot=(rows, cols, si)) #fig.add_subplot(si)
        exp_polar_plot(t, r, e, cmap=cmap, fig=fig, ax=ax)
        set_rho_lim(r, ax)
        
        if type(gene) == list:
            ax.set_title(', '.join(gene))    
        else:
            ax.set_title(gene)
        
        si += 1
    
    libplot.add_colorbar(fig, cmap)
    
    return fig


def plot_genes_longitudinal(data, time, genes):
    if type(genes) is pd.core.frame.DataFrame:
        genes = genes['Genes'].values
        
    time_idx = np.argsort(time)
    
    ts = time[time_idx]
        
    rows = int(np.ceil(np.sqrt(len(genes))))
    
    w = 4 * rows
    h = 3 * rows
    
    fig = libplot.new_base_fig(w=w, h=h)
    
    for i in range(0, len(genes)):
        gene = genes[i]
        
        if type(gene) == list:
            d = np.zeros(data.shape[1])
            
            for g in gene:
                print(g)
                d += data.loc[data.index.str.endswith(g), :].T.iloc[:, 0].values
            
            #exp /= len(gene)
        else:
            d = data.loc[data.index.str.endswith(gene), :].T.iloc[:,0].values
        
        ds = d[time_idx]
        
        ax = libplot.new_ax(fig, subplot=(rows, rows, i + 1))
        
        for j in range(1, 7):
            ax.axvline(x=j, linewidth=0.5, color=libplot.TRANS_GRAY, linestyle='--')
        
        libplot.scatter(ts, ds, ax=ax, c=libplot.TRANS_GRAY, alpha=0.25)
    
        ge_bins, ge_std_bins, t_bins = bin_time_data(ds, ts)
        
        xs = np.linspace(0, 7, BINS)
        
        # interpolate the means and the standard deviation
        
        spl = UnivariateSpline(t_bins, ge_bins, ext='extrapolate', k=K)
        ys = spl(xs)
        
        spl = UnivariateSpline(t_bins, ge_bins - ge_std_bins, ext='extrapolate', k=K)
        ys1 = spl(xs)
        
        spl = UnivariateSpline(t_bins, ge_bins + ge_std_bins, ext='extrapolate', k=K)
        ys2 = spl(xs)
        
        libplot.plot(xs, ys, label='Gene Expression', ax=ax)
        #libplot.plot(t_bins, ge_bins, label='Gene Expression', ax=ax)
        ax.fill_between(xs, ys1, ys2, alpha=0.2)
        
        if type(gene) == list:
            ax.set_title(', '.join(gene))    
        else:
            ax.set_title(gene)
        
        ax.set_xlim([0, TWO_PI])
        ax.set_ylim([0, 4])
        ax.set_xlabel(r'Time ($\theta$)')
        ax.set_ylabel('Gene Expression')
        
    return fig

def to_polar(data, r=0, m=1):
    """
    Convert tsne coordinates to polar
    """
    
    n = data.shape[0]
    
    theta = np.zeros(n)
    rho = np.zeros(n)
    
    for i in range(0, n):
        x = libplot.cart2pol(data.iloc[i, 0], data.iloc[i, 1])
        
        theta[i] = x[0]
        
        if theta[i] < 0:
            theta[i] = TWO_PI + theta[i]
        
        rho[i] = x[1] * m + r
        
    return theta, rho