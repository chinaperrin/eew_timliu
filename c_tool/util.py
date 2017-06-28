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


# home directory
HOME = '/Users/Timothy/Desktop/SURF2017/eew_timliu'

import os

def separ(cond, value, single, prop):
    '''separator function. Creates a subset of a single data file with
    a certain property.
    inputs: cond - condition to test for 0 <, 1 >, 2 ==
            value- value to compare to
            single - name of single to sort through
            prop - name of file with property we are interested in.'''
    
    con_map = {0: 'g', 1:'l', 2:'is'}  #mapping between cond and strings
    
    #go to the single we're interested in
    os.chdir(os.path.join(HOME, 'ex_data', single[:-7]))
    
    singlef = open(single, 'r')  #open single file data
    s_list = singlef.readlines() #read in lines from single file
    
    singlef.close()              #close the single file
    
    #go to folder with property we're interested in - strip _01.txt
    os.chdir(os.path.join(HOME, 'ex_data', prop[:-7]))    
    
    propf = open(prop, 'r')      #file with property we care about
    prop_list = propf.readlines() #read the property file
    propf.close()    
    
    #see which data points we're interested in
    if cond == 0:
        keep = [i for i in range(len(prop_list)) if float(prop_list[i][:-1]) < value]
    if cond == 1:
        keep = [i for i in range(len(prop_list)) if float(prop_list[i][:-1]) > value]
    if cond == 2:
        keep = [i for i in range(len(prop_list)) if prop_list[i][:-1] == value]
    else:
        print("Invalid Condition!!")
        return
        
    #name offolder with subset
    fpath = os.path.join(HOME, 'ex_data', single[:-7] + '_' + prop[:-7] + '_'\
                         + con_map[cond] + '_' + value)
    
    #make new directory if it doesn't already exist
    if not os.path.exists(fpath):
        os.makedirs(fpath)
        
    os.chdir(fpath)
    #create the new list of data that meets set condition
    n_list = [s_list[i] for i in range(len(s_list)) if i in keep]
    #new file that's subset of original single
    newf = open(single[:-7] +'_' + prop[:-7] + '_' + con_map[cond] +\
                '_' + value + '_'+ single[-6:-4] + '.txt', 'w')
    for line in n_list:
        newf.write(line)
    newf.close()
    #change directory back to where we started
    os.chdir(os.path.join(HOME, 'c_tool'))
    return

def compare(est, act, comp = 0):
    '''compares the estimated and actual values of some parameter
    inputs: est - file of estimated values
            act - file of actual values
            comp - integer indicating if more computing is needed on the 
                  error for example modding
                  1 - mod 180 for back azimuth
    outputs: file with list of estimation errors'''
    
    #directory with estimates
    os.chdir(os.path.join(HOME, 'ex_data', est[:-7]))
    ef = open(est, 'r')  #estimated file
    el = ef.readlines()  #read all the estimates
    eflt = [float(x[:-2]) for x in el]  #convert to floats
    ef.close()
    #directory with actual value
    os.chdir(os.path.join(HOME, 'ex_data', act[:-7]))    
    af = open(act, 'r')  #actual file
    al = af.readlines()  #read all the actual measurements
    aflt = [float(x[:-2]) for x in al]   #convert to floats
    af.close()
    #directory to store error files
    os.chdir(os.path.join(HOME, 'errors'))
    erf = open(est[:-4] + '_err.txt', 'w')    
    e_list = [(eflt[i] - aflt[i]) for i in range(len(eflt))]
    
    #do additional computations if spcified
    if comp == 1:             #1 indicates handle bazi wraparound
        e_list = bazi_ad(e_list)
        
    e_list = ['%.4f\n' %x for x in e_list]  #write as a string
    for line in e_list:
        erf.write(line)
    erf.close()
    
    os.chdir(os.path.join(HOME, 'laz'))   #return to where we started
    
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