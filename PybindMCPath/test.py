import numpy as np

from datetime import datetime as dtm
from datetime import timedelta as td

import time
import os
import sys
from multiprocessing import Pool

sys.path.append("bin")
import MCPath 

today = (24,2,2020)
todaydt = dtm(2020,2,24)
num = 250000
steps = 366
tenor = steps/365

ir_type = 1
ir_term = np.array([0,31,94,182,276,367])
ir_data = np.array([2.61,2.61,2.76,2.79,2.81,2.83])/100
ir_dc = 0

d_type = 0
d_term = np.array([0])
d_data = np.array([0.01])
d_dc = 0

v_type = 1
v_term = np.array([1,30,63,92,124,154,183,215,245,274,306,336,366])
v_data = np.array([36.97,36.97,32.88,32.60,31.62,29.16,28.44,27.57,29.29,28.19,27.70,27.45,27.55])/100
v_dc = 0

proc_type = 0

upout_obday = np.array([29,60,91,121,151,182,213,245,274,304,336,366])
upout_obidx = np.isin(np.arange(0,steps+1),upout_obday)
upout_barrier = np.zeros(steps+1)
upout_barrier[upout_obidx] = np.array([1.03,1.03,1.03,1.02,1.02,1.02,1.01,1.01,1.01,1.00,1.00,1.00])

downout_obday = np.array([x for x in range(1,steps+1) if (todaydt+td(days=x)).isoweekday()<6])
downout_obidx = np.isin(np.arange(0,steps+1),downout_obday)
downout_barrier = np.array([0.7])


def MPMC(task):
    array = np.zeros((task[0],steps+1))
    return MCPath.GeneratePath(today,task[0],steps,tenor,
                               ir_type,ir_term,ir_data,ir_dc,
                               d_type,d_term,d_data,d_dc,
                               v_type,v_term,v_data,v_dc,
                               task[1],upout_obidx,upout_barrier,
                               task[2],downout_obidx,downout_barrier,
                               proc_type, array,
                               True,0,42)

if __name__ == '__main__':
    proc_type = 0
    input_array = np.zeros((num,steps+1))

    print("Test generating MC paths...")
    upout_type,downout_type = 0,0
    t0 = time.time()
    res=MCPath.GeneratePath(today,num,steps,tenor,
                            ir_type,ir_term,ir_data,ir_dc,
                            d_type,d_term,d_data,d_dc,
                            v_type,v_term,v_data,v_dc,
                            upout_type,upout_obidx,upout_barrier,
                            downout_type,downout_obidx,downout_barrier,
                            proc_type,input_array,
                            True,0,np.random.randint(0,99999))
    t1 = time.time()
    print(" [Result]: ",t1-t0)

    print("Test generating MC paths with early stop barrier ...")
    upout_type,downout_type = 2,1
    t2 = time.time()    
    res=MCPath.GeneratePath(today,num,steps,tenor,
                            ir_type,ir_term,ir_data,ir_dc,
                            d_type,d_term,d_data,d_dc,
                            v_type,v_term,v_data,v_dc,
                            upout_type,upout_obidx,upout_barrier,
                            downout_type,downout_obidx,downout_barrier,
                            proc_type,input_array,
                            True,0,np.random.randint(0,99999))
    t3 = time.time()
    print(" [Result]: ",t3-t2)

    #=========================
    #  Multi Processing Test
    #=========================

    n_proc = 4

    upout_type,downout_type = 0,0
    print(f"Test generating MC paths with {n_proc} processes...")
    t4 = time.time()
    with Pool(n_proc) as p:
        res_list = p.map(MPMC, [(int(num/n_proc),upout_type,downout_type) for x in range(n_proc)])
    res = np.concatenate(res_list,axis=0)
    print(" [Result]: ",time.time()-t4) 

    
    upout_type,downout_type = 2,1
    print(f"Test generating MC paths with {n_proc} processes and early stop barrier...")
    t5 = time.time()
    with Pool(n_proc) as p:
        res_list = p.map(MPMC, [(int(num/n_proc),upout_type,downout_type) for x in range(n_proc)])
    res = np.concatenate(res_list,axis=0)
    print(" [Result]: ",time.time()-t5)
    print(res[:,-1].mean())
    os.system("pause")