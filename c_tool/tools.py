# tools.py
# used to repeatedly call functions in separ.py, grapher.py, and extra.py
#
# Revision History
# 6/27/17    Tim Liu    started file

from util import *
from grapher import *
from extra import *

def call():
    '''used to repeatedly call different comparator utilities'''
    
    y_list = ['deichmann_sentype_is_H_01_abs_err.txt',\
              'la05_sentype_is_H_01_abs_err.txt',\
              'onsiteBazi_sentype_is_H_01_abs_err.txt',\
              'sac_sentype_is_H_01_abs_err.txt']
    x_list = ['epi_sentype_is_H_01.txt', 'mag_sentype_is_H_01.txt', \
              'snr_sentype_is_H_01.txt']
    i = 0
    for y in y_list:
        for x in x_list:
            two_d(x, y, y[:-28] + 'vs.' + x[:4], x[:4], 'Back Azimuth Error')
            i += 1
            print('Graph: ', i)
        
        
    return


'''
        a1 = 'onsiteBazi_%d.txt' %i
        a2 = 'deichmann_%d.txt' %i
        a3 = 'sac_%d.txt' %i
        a4 = 'la05_%d.txt' %i        
        separ(2, 'LN', a1, 'sentype_01.txt')
        separ(2, 'LN', a2, 'sentype_01.txt')
        separ(2, 'LN', a3, 'sentype_01.txt')
        separ(2, 'LN', a4, 'sentype_01.txt')
        print('Completed: ', i)
        '''