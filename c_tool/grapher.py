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
# 6/29/17    Tim Liu    started three_d - still not complete
# 7/7/17     Tim Liu    added regression line to two_d
# 7/10/17    Tim Liu    finished three_d and added 10% and 90% lines
#
# Table of Contents
# one_d    plots one dimensional data
# two_d    plots two dimensional data as a scatter plot
# f2l      converts a file into a list of data
# scale    creates bin sizes and bin range for data
# three_d  plots 3-d data (not finished)


import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os
from scipy import stats
import pandas as pd


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
    plt.savefig(err[:-4])
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
        
        
    plt.scatter(datax, datay, color = color, s = 0.5)
    plt.title(title, fontsize = 18)    
    
    plt.ylabel(ylab, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)    
    plt.grid(True)


    
    # Calculate best fit line
    (slope, intercept, r_value, p_value, std_err) = stats.linregress(datax, datay)
    coefficients = ([slope, intercept])
    p = np.poly1d(coefficients)
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    t1 = np.arange(xmin, xmax, 1.0)
    
    '''
    # Display equation and r-value of graph
    plt.plot(t1, p(t1), 'b-', linewidth = 3.0)
    plt.text(0.7*xmax, 0.9*ymax, r"y = " + str(round(slope, 3)) + "x + " + str(round(intercept, 2)), fontsize=14)
    plt.text(0.7*xmax, 0.8*ymax, r"r =" + str(round(r_value, 3)), fontsize = 14)
    r = np.poly1d([1, 0])
    '''
    
    '''plot the percentiles'''
    a = np.arange(0, 100, 5)
    b = [x + 2.5 for x in a]
    up = []
    down = []
    mid = []
    for j in range(len(a)-1):
        chunks = [datay[x] for x in range(len(datax)) if (datax[x] > a[j] and datax[x] <a[j+1])]
        up.append(np.percentile(chunks, 90))
        down.append(np.percentile(chunks, 10))
        mid.append(np.percentile(chunks, 50))
    plt.plot(b[:-1], up, '-', color='g')
    plt.plot(b[:-1], down, '', color='r')
    plt.plot(b[:-1], mid, '-', color='b')
    plt.plot([0, 100], [0, 0], '-', color='black')
    os.chdir(os.path.join(HOME, 'g_out'))  #directory to dump images
    plt.savefig('GBA_errorv. GBA est' + title[19:22] + '.png', format='png')
    os.chdir(os.path.join(HOME, 'c_tool')) #return to directory
            
    return

def three_d(varx, vary, varz, xlab, ylab, zlab, title, filename, subplot=0, pathy='errors'):
    #find data to x variable
    os.chdir(os.path.join(HOME, 'ex_data', varx[:-7]))
    datax = f2l(varx)    
    
    #find data to y variable
    os.chdir(os.path.join(HOME, pathy))   
    datay = f2l(vary)                 #convert file to a list
    
    #find data to z variable
    os.chdir(os.path.join(HOME, 'ex_data', varz[:-7]))
    dataz = f2l(varz)     
    
    #create figure
    plt.figure(vary[:-4] + varx[:-4])
    #plt.figure(1)
    if subplot != 0:
        plt.subplot(subplot)   #optionally make subplots
        print(subplot)
    
    df = pd.DataFrame({"x":datax, "y":datay, "c":dataz})   
    cmap = plt.cm.rainbow
    norm=matplotlib.colors.Normalize(vmin=4, vmax=6.0)    
      
    plt.scatter(df.x, df.y, color=cmap(norm(df.c)), s = 1)
    
    plt.title(title, fontsize = 18)        
    plt.ylabel(ylab, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)    
    plt.grid(True)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cmap = plt.colorbar(sm)
    cmap.set_label(zlab)
    
    a = np.arange(0, 100, 5)
    b = [x + 2.5 for x in a]
    up = []
    down = []
    mid = []
    for j in range(len(a)-1):
        chunks = [datay[x] for x in range(len(datax)) if (datax[x] > a[j] and datax[x] <a[j+1])]
        up.append(np.percentile(chunks, 90))
        down.append(np.percentile(chunks, 10))
        mid.append(np.percentile(chunks, 50))
    plt.plot(b[:-1], up, '-', color='g')
    plt.plot(b[:-1], down, '', color='r')
    plt.plot(b[:-1], mid, '-', color='b')
    plt.plot([0, 100], [0, 0], '-', color='black')
    plt.xlim(0, 120)
    plt.ylim(-100, 100)

    os.chdir(os.path.join(HOME, 'g_out'))  #directory to dump images
    plt.savefig(filename + '.png', format='png')
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