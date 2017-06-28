# extra.py
# extractor pulls certain subsets of data from the matlab tracelists

# Revision History
# 6/26/17    Tim Liu    started file
# 6/26/17    Tim Liu    wrote main and gen_script
# 6/27/17    Tim Liu    updated comments and renamed "main" to "extract"
# 6/27/17    Tim Liu    updated HOME to reflect new directory structure
#
# Table of Contents
# extract     extract data with certain object field under a specified name
# gen_script  generates a script for extracting data

import os
import sys
import matlab.engine

HOME = '/Users/Timothy/Desktop/SURF2017/eew_timliu'



def extract(var, name, i):
    '''extract function. Calls other functions to extract data from matlab
    tracelist.
    input: var - variable to be extracted     
           name - name of file
           i - number of subfields in var
    output: none'''
    
    gen_script(var, name, i)    #call function to generate extraction script
    
    #go to directory with extraction script
    os.chdir(os.path.join(HOME, 'code'))
    
    eng = matlab.engine.start_matlab()
    eng.extractor(nargout=0)           #run extractor script
        
    os.remove("extractor.m")           #delete the script
    
    os.chdir(os.path.join(HOME, 'c_tool')) #switch back to directory
        
    return

def gen_script(var, name, i):
    '''generates a script for extracting specified piece of data
    input: var - variable the script is interested in extracting
           name - name of variable and title of files
           i - number of different subfields
    output: matlab script for extracting data'''
    
    os.chdir('/Users/Timothy/Desktop/SURF2017/tracelists/code/')
    
    temp = open('extract_temp.m', 'r')    #open extractor template
    
    extract = open('extractor.m', 'w')    #open extractor script
    t_string = temp.read()
    q_ind = t_string.find('Q')   #find name field
    j_ind = t_string.find('J')   #find length field
    k_ind = t_string.find('K')   #find data field
    
    script = t_string[0:q_ind]         #beginning of script up to Q
    script += "'%s'"%name              #add name of variable
    script += t_string[q_ind+1:j_ind]  #next chunk of script
    script += str(i)                   #add number of different fields
    script += t_string[j_ind+1:k_ind]  #next chunk of script
    script += var                      #write class subfield
    script += t_string[k_ind+1:]       #last chunk of script
    
    
    extract.write(script)              #write the script
    
    temp.close()                       #close out
    extract.close()
    
    #return to original directory
    os.chdir('/Users/Timothy/Desktop/SURF2017/c_tool')
    
    return