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
    
    g_list = []
    for i in range(1,10):
        a1 = 'onsiteBazi_sentype_is_LN_0%d.txt' %i
        a2 = 'deichmann_sentype_is_LN_0%d.txt' %i
        a3 = 'sac_sentype_is_LN_0%d.txt' %i
        a4 = 'la05_sentype_is_LN_0%d.txt' %i 
        
        
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