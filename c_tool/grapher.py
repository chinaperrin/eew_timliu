# grapher
# generic plotter with functions for plotting different sets of data
# can plot errors (1-D), correlations (2-D) and correlations with 
# another variable denoted by color (3-D)
#
# Revision History
# 6/27/17    Tim Liu    started file
# 6/27/17    Tim Liu    wrote one_d
# 6/27/17    Tim Liu    updated HOME to reflect new file system
# 6/28/17    Tim Liu    wrote two_d
#
# Table of Contents
# one_d    plots one dimensional data
# two_d    plots two dimensional data as a scatter plot
# f2l      converts a file into a list of data
# scale    creates bin sizes and bin range for data



import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import os

HOME = '/Users/Timothy/Desktop/SURF2017/eew_timliu'


def one_d(err, title, x_label, subplot=0, color = 'deepskyblue', path = 'errors'):
    '''plots one dimensional data - most likely error data.
    inputs: err - file with a list of errors
            path - path to data to graph
            fig_name - name of figure
            title - tile of figure
            x_label - x axis label
            subplot - index of subplot if multiple plots in figure
            color - color of graph; default set'''
    
    #find data to graph
    os.chdir(os.path.join(HOME, path))
    
    data = f2l(err)                 #convert file to a list
    abs_er = [abs(x) for x in data] #list of absolute errors
    
    mean_err = np.mean(abs_er)
    std_err = np.std(abs_er)
    
    #create figure
    plt.figure(err[:-4])
    if subplot != 0:
        plt.subplot(subplot)   #optionally make subplots
    
    #call function to determine scale of bins
    dat_min, dat_max, dat_stp = scale(data)
    plt.hist(data, bins = range(dat_min, dat_max, dat_stp), color=color)
    plt.title(title, fontsize = 18)    
    
    plt.ylabel('Counts', fontsize = 18)
    plt.xlabel(x_label, fontsize = 18)    
    plt.grid(True)

    os.chdir(os.path.join(HOME, 'g_out'))  #directory to dump images
    plt.savefig()
    os.chdir(os.path.join(HOME, 'c_tool')) #return to directory
    
    #function returns - does not actually show plot!
    return

def two_d(varx, vary, title, xlab, ylab, subplot=0, color = 'deepskyblue', pathy = 'errors'):
    '''plots two different variables
    inputs: varx - .txt file with x variable to plot
            vary - .txt file with y variable to plot
            xlab - label for x axis
            ylab - label for y axis
            subplot - index of subplot if multiple plots in figure
            color - color of graph; default set
            pathx - path to x variable'''
    
    #find data to x variable
    os.chdir(os.path.join(HOME, 'ex_data', varx[:-7]))
    datax = f2l(varx)    
    
    #find data to y variable
    os.chdir(os.path.join(HOME, pathy))   
    datay = f2l(vary)                 #convert file to a list  
    
    #create figure
    plt.figure(vary[:-4] + varx[:-4])
    if subplot != 0:
        plt.subplot(subplot)   #optionally make subplots    
        
        
    plt.scatter(datax, datay, color = color, s = 0.1)
    plt.title(title, fontsize = 18)    
    
    plt.ylabel(ylab, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)    
    plt.grid(True)

    os.chdir(os.path.join(HOME, 'g_out'))  #directory to dump images
    plt.savefig(vary[:-4] + varx[:-4])
    os.chdir(os.path.join(HOME, 'c_tool')) #return to directory
            
    return
    
def f2l(f):
    '''helper function converts a file to a list of data.
    inputs: f - name of file
    outputs l - list of floats from the file passed through.'''
    
    f_h = open(f, 'r')                 #file handle
    f_sl = f_h.readlines()             #file as a list of strings
    
    l = [float(x[:-1]) for x in f_sl]  #strip newline char and convert to float
    f_h.close()
    
    
    return l

def scale(data):
    '''calculates appropriate min, max, and bin size
    inputs: data - list of data
    outputs: bin_par - (min, max, step)'''
    
    if max(data) < 3:   #get range of data and find how to scale it
        dec = 1
        dat_stp = .1
    elif max(data) < 30 and max(data) > 3:
        dec = 0
        dat_stp = 1
    elif max(data) > 30:
        dec = -1
        dat_stp = 10
    
    dat_min = int(round(1.05 * min(data), dec))  #find lower bound of chart
    dat_max = int(round(1.05 * max(data), dec))   #find upper bound of chart
    
    return (dat_min, dat_max, dat_stp)