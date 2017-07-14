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
    for i in range(10):
        f = 'GBA_onset_0%d_ex10_err.txt' %i
        g = 'GBA_dis_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
        h = 'GBA_mag_act_GBA_onset_0%d_e_l_1.0_01.txt' %(i)
        three_d(g, f, h, 'Actual Distance (km)', 'Estimate - Actual(km)', 'Actual Mag.', 'GBA Distance Error %1.1f sec Since Onset' %(i*0.5), 'GBA_diserr_actmag%1.1f' %(i*0.5))
        print(i)
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

'''    y_list = ['deichmann_sentype_is_H_01_abs_err.txt',\
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
            print('Graph: ', i)'''

'''    charts = ['onsiteBazi_sentype_is_H_01_err.txt', \
              'deichmann_sentype_is_H_01_err.txt', \
              'sac_sentype_is_H_01_err.txt', \
              'la05_sentype_is_H_01_err.txt']
    
    
    for dat in charts:
        one_d(dat, 'Back Azimuth Errors', 'Degrees')'''





'''two_d('GBA_dis_act_01.txt', 'GBA_hat2_10_ex10.txt', 'Estimated v. Actual Distance' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_hat2')      
two_d('GBA_dis_act_01.txt', 'GBA_hat2_20_ex10.txt', 'Estimated v. Actual Distance' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_hat2')   
two_d('GBA_dis_act_01.txt', 'GBA_hat2_30_ex10.txt', 'Estimated v. Actual Distance' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_hat2')
two_d('GBA_dis_act_01.txt', 'GBA_hat2_40_ex10.txt', 'Estimated v. Actual Distance' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_hat2')'''

'''    for i in range(10):
        f = 'GBA_onset_0%d.txt' %i
        convert(f, 'ex10')
    f = 'GBA_onset_10.txt'
    convert(f, 'ex10')'''

'''    two_d('GBA_dis_act_01.txt', 'GBA_onset_00_ex10.txt', 'Estimated v. Actual Distance 0 sec' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_onset')
    two_d('GBA_dis_act_01.txt', 'GBA_onset_04_ex10.txt', 'Estimated v. Actual Distance 2 sec' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_onset') 
    two_d('GBA_dis_act_01.txt', 'GBA_onset_08_ex10.txt', 'Estimated v. Actual Distance 4 sec' ,'Actual (km)', 'Estimated (km)', pathy = '/Users/Timothy/Desktop/SURF2017/eew_timliu/ex_data/GBA_onset') '''

'''    for i in range(10):
        f = 'GBA_onset_0%d_ex10_err.txt' %i
        g = 'GBA_dis_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
        two_d(g, f, 'GBA Distance Error %1.1f sec Since Onset' %(i*0.5), 'Actual\
        Distance (km)', 'Estimate Error')
            
    return'''

'''    for i in range(10):
        f = 'GBA_onset_0%d_ex10_err.txt' %i
        g = 'GBA_dis_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
        two_d(g, f, 'GBA Distance Error %1.1f sec Since Onset' %(i*0.5), 'Actual Distance (km)', 'Estimate Error (km)')
        print(i)
            
    return'''

'''    for i in range(10):
        f = 'GBA_onset_0%d_ex10.txt' %i    
        separ(1, 1.00, 'GBA_mag_act_01.txt', f)'''

'''    for i in range(10):
        f = 'GBA_onset_0%d_ex10_err.txt' %i
        g = 'GBA_dis_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
        h = 'GBA_mag_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
        three_d(g, f, h, 'Actual Distance (km)', 'Estimate Error (km)', 'Magnitude', 'GBA Distance Error %1.1f sec Since Onset' %(i*0.5), )
        print(i)'''

'''    for i in range(10):
        f = 'GBA_mag_onset_0%d.txt' %i
        compare(f, 'GBA_mag_act_01.txt', sum_file = 'sum_mag_onset.txt')'''



'''for a in range(4) :
    b = [0, 2, 6, 8]
    subs = [221, 222, 223, 224]
    i = b[a]
    sub = subs[a]
    f = 'GBA_onset_0%d_ex10_err.txt' %i
    g = 'GBA_dis_act_GBA_onset_0%d_e_l_1.0_01.txt' %i
    h = 'GBA_mag_onset_0%d__GBA_onset_0%d_e_l_1.0_rr.txt' %(i, i)
    three_d(g, f, h, 'Actual Distance (km)', 'Estimate - Actual(km)', 'Est-Actual Mag.', 'GBA Distance Error %1.1f sec Since Onset' %(i*0.5), 'GBA_dismagerr_act%1.1f' %(i*0.5), subplot = sub)
    print(i)
return'''