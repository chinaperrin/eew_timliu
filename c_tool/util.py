# util.py
# separator that filters out certain properties from data_lists
# 
# Revision History
# 6/26/17    Tim Liu    started file
# 6/26/17    Tim Liu    wrote separ
# 6/26/17    Tim Liu    added compare function
# 6/27/17    Tim Liu    updated compare and separ function to navigate
#                       directories
# 6/27/17    Tim Liu    added option for corrections and other operations
#                       to the compare function
# 6/27/17    Tim Liu    wrote bazi_ad to handle bazi wraparound
# 6/27/17    Tim Liu    changed name to "util.py" for utilities
# 6/27/17    Tim Liu    updated HOME to reflect new file system
# 6/27/17    Tim Liu    compare also writes file with absolute value errors
# 6/28/17    Tim Liu    compare can optionally write to summary file
# 6/28/17    Tim Liu    fixed bug where compare overwrote summary file
# 6/28/17    Tim Liu    rewrote compare to enhance speed
# 6/29/17    Tim Liu    added comp argument to compare
# 6/30/17    Tim Liu    updated comments in compare
# 7/06/17    Tim Liu    added convert function
# 7/07/17    Tim Liu    added by_event function for splitting by event
# 7/08/17    Tim Liu    changed compare to throw out estimates of 1.000
# 7/08/17    Tim Liu    fixed bug in separ that would throw invalid condition
#                       if cond != 2; changed if to elif


# home directory
HOME = '/Users/Timothy/Desktop/SURF2017/eew_timliu'

import os
import numpy as np
import math

def separ(cond, value, single, prop, comp = 0):
    '''separator function. Creates a subset of a single data file with
    a certain property.
    inputs: cond - condition to test for 0 <, 1 >, 2 ==
            value- value to compare to
            single - name of single to sort through
            prop - name of file with property we are interested in.
            comp - optional argument; if !=0 perform some operation on 
                   data before copying over'''
    
    con_map = {0: 'g', 1:'l', 2:'is'}  #mapping between cond and strings
    
    #go to the single we're interested in
    #os.chdir(os.path.join(HOME, 'ex_data', single[:-7]))
    os.chdir(os.path.join(HOME, 'errors'))
    
    singlef = open(single, 'r')  #open single file data
    s_list = singlef.readlines() #read in lines from single file
    
    singlef.close()              #close the single file
    
    #go to folder with property we're interested in - strip _01.txt
    os.chdir(os.path.join(HOME, 'ex_data', prop[:-12]))    
    
    propf = open(prop, 'r')      #file with property we care about
    prop_list = propf.readlines() #read the property file
    propf.close()    
    
    #see which data points we're interested in
    if cond == 0:
        keep = [i for i in range(len(prop_list)) if float(prop_list[i][:-1]) < value]
    elif cond == 1:
        keep = [i for i in range(len(prop_list)) if float(prop_list[i][:-1]) > value]
    elif cond == 2:
        keep = [i for i in range(len(prop_list)) if prop_list[i][:-1] in value]
    else:
        print("Invalid Condition!!")
        os.chdir(os.path.join(HOME, 'c_tool'))        
        return
        
    #name of folder with subset
    fpath = os.path.join(HOME, 'ex_data', single[:-7] + '_' + prop[:-7] + '_'\
                         + con_map[cond] + '_' + str(value))
    
    #make new directory if it doesn't already exist
    if not os.path.exists(fpath):
        os.makedirs(fpath)    
    os.chdir(fpath)
    #create the new list of data that meets set condition
    
    n_list = []  #new list
    
    for i in range(len(s_list)):  #go through s_list and copy the value if
        if len(keep) == 0:        #index is the first value in keep
            break        
        if i == keep[0]:
            n_list.append(s_list[i])
            del keep[0]   #remove first element of keep so first element
                          #is always the next to be added to s_list
    
    #perform computation when separating - comp = 1 means take log10 of data        
    if comp == 1:
        n_list = [str(math.log10(float(x[:-1]))) + '\n' for x in n_list]
        
    #new file that's subset of original single
    newf = open(single[:-7] +'_' + prop[:-7] + '_' + con_map[cond] +\
                '_' + str(value) + '_'+ single[-6:-4] + '.txt', 'w')
    for line in n_list:
        newf.write(line)
    newf.close()
    #change directory back to where we started
    os.chdir(os.path.join(HOME, 'c_tool'))
    
    return

def compare(est, act, comp = 0, sum_file = ''):
    '''compares the estimated and actual values of some parameter
    inputs: est - file of estimated values
            act - file of actual values
            comp - integer indicating if more computing is needed on the 
                  error for example modding
                  1 - mod 180 for back azimuth
    outputs: file with list of estimation errors'''
    
    #directory with estimates
    os.chdir(os.path.join(HOME, 'ex_data', est[:-7]))
    #os.chdir(os.path.join(HOME, 'ex_data', est[:-12]))
    ef = open(est, 'r')  #estimated file
    el = ef.readlines()  #read all the estimates
    eflt = [float(x[:-1]) for x in el]  #convert to floats
    ef.close()
    #directory with actual value
    os.chdir(os.path.join(HOME, 'ex_data', act[:-7]))    
    af = open(act, 'r')  #actual file
    al = af.readlines()  #read all the actual measurements
    aflt = [float(x[:-1]) for x in al]   #convert to floats
    af.close()
    #directory to store error files
    os.chdir(os.path.join(HOME, 'errors'))
    erf = open(est[:-4] + '_err.txt', 'w')
    erfa = open(est[:-4] + '_abs_err.txt', 'w')
    e_list = [(eflt[i] - aflt[i]) for i in range(len(eflt)) if eflt[i] != 1.0]
    
    #do additional computations if spcified
    if comp == 1:             #1 indicates handle bazi wraparound
        e_list = bazi_ad(e_list)
        
    abs_list = [abs(x) for x in e_list]        #find absolute values        
    ea_list = ['%.4f\n' %x for x in abs_list]  #write absolute values
    e_list = ['%.4f\n' %x for x in e_list]     #write as a string
    
    for line in e_list:
        erf.write(line)
    erf.close()
    
    for line in ea_list:
        erfa.write(line)
    erfa.close()
    
    #write to summary file - optional
    if sum_file != '':
        s_file = open(sum_file, 'a')
        s_file.write('%s\n' %est)
        s_file.write('Mean Error: %.4f\n'%np.mean(abs_list))
        s_file.close()
        
    
    os.chdir(os.path.join(HOME, 'c_tool'))   #return to where we started
    
    return

def bazi_ad(baz_list):
    '''adjusts back azimuth so error spans -180 to 180 instead of
    -360 to 360
    inputs:  baz_list - list of back azimuths unadjusted
    outputs: baz_cor - list of corrected back azimuths'''
    
    baz_cor = []
    
    for err in baz_list:
        if err > 180.0:                      #handle wraparound
            err -= 360.0
        if err < -180.0:
            err += 360.0
        baz_cor.append(err)
        
        
    return baz_cor

def convert(filename, comp):
    '''opens a textfile with data and modifies data. Creates new file
    with the modified data without modifiying existing file.
    inputs: filename - name of file to modify
            comp - string encoding what operation to perform'''
    #go to where extracted data lives    
    os.chdir(os.path.join(HOME, 'ex_data', filename[:-7]))
    
    fo = open(filename, 'r')                    #original file
    dat_o = fo.readlines()                      #read in the file
    dat_float = [float(x[:-1]) for x in dat_o]  #convert to floats
    
    #perform computation and open and name file
    if comp == 'ex10':
        #exponentiate to power 10
        dat_n = ['%.3f\n' %(10**x) for x in dat_float]      
        fn = open(filename[:-4] + '_ex10' + '.txt', 'w')
     
    #write all the data to the new file   
    for line in dat_n:
        fn.write(line)
        
    fo.close()    #close original file
    fn.close()    #close new file
    #change directory back to where we started
    os.chdir(os.path.join(HOME, 'c_tool'))
    
def by_event(filename, title, event_list = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/event/event.txt'):
    '''breaks data by event. All reporting stations from an event are grouped
    together
    inputs: filename - name of file to break up
            event_list - .txt file with list of event names
    outputs: writes one file collecting all from same event'''
    
    #go to where extracted data lives    
    os.chdir(os.path.join(HOME, 'ex_data', filename[:-7]))
    #open the list of events
    fe = open(event_list, 'r')
    #open the file with the data
    f = open(filename, 'r')
    
    events = fe.readlines()   #read off all the events
    data_list = f.readlines()      #read all the dat to be split
    
    #current event we're reporting on
    current = events[0]
    se = open(current + title + '.txt', 'w')
    
    #go through the data list
    for i in range(len(data_list)):
        #event isn't the same - switch events
        if events[i] != current:
            current = event[i]
            se.close()
            se = open(current + title + '.txt', 'w')
        #write to the file
        se.write(data_list[i])
        
    fe.close()
    f.close()
    se.close()
        
    #change directory back to where we started
    os.chdir(os.path.join(HOME, 'c_tool'))    
    
    return