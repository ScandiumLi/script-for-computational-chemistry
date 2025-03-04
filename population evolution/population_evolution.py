# -*- coding: utf-8 -*-
"""
@Time : 2024/06/23/ 21:40
@Auth : Scandium, kanglihn@163.com
@File :  population_evolution.py
@IDE : Pycharm

"""

import numpy as np
from scipy.integrate import odeint
import json


total_time = 5.0e-5 # 总的模拟时间, 单位为s
population_0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # 光致发光设S1初始值为1
file_name = "rate_constant.json" # 提供速率常数的文件
out_file = "Au28.txt" # 保存结果文件


def lorenz2(electric_state, time):
    state_S1 = electric_state[0]
    state_T1 = electric_state[1]
    state_T2 = electric_state[2]
    state_S0 = electric_state[3]
	
    # 读取速率常数的文件
    with open(file_name, 'r') as file:
        rate_constant = json.load(file)

    # 获取速率常数, 单位为s^-1
    # 辐射速率常数
    k_f = rate_constant["k_f"]
    k_pT1 = rate_constant["k_pT1"]
    
     # 非辐射速率常数
    k_S1T1 = rate_constant["k_S1T1"]
    k_T1S1 = rate_constant["k_T1S1"]
    k_S1T2 = rate_constant["k_S1T2"]
    k_T2S1 = rate_constant["k_T2S1"]
    k_T2T1 = rate_constant["k_T2T1"]
    k_T1T2 = rate_constant["k_T1T2"]  
    k_S1S0 = rate_constant["k_S1S0"]  # S1->S0内转换(IC)速率常数
    k_T1S0 = rate_constant["k_T1S0"]  # T1->S0系间窜越(ISC)速率常数
    
    # 质量作用定律方程
    d_S1_dt = k_T1S1 * state_T1 + k_T2S1 * state_T2 - (k_S1T1 + k_S1T2 + k_S1S0 + k_f) * state_S1
    d_T1_dt = k_T2T1 * state_T2 + k_S1T1 * state_S1 - (k_T1S1 + k_pT1 + k_T1T2 + k_T1S0) * state_T1
    d_T2_dt = k_S1T2 * state_S1 + k_T1T2 * state_T1 - (k_T2S1 + k_T2T1) * state_T2
    d_S0_dt = (k_f + k_S1S0) * state_S1 + (k_pT1 + k_T1S0) * state_T1
    d_S0S1_dt = (k_f + k_S1S0) * state_S1
    d_S0T1_dt = (k_pT1 + k_T1S0) * state_T1

    return [d_S1_dt, d_T1_dt, d_T2_dt, d_S0_dt, d_S0S1_dt, d_S0T1_dt]

time_0 = 0.0
time_end = total_time
time = np.linspace(time_0, time_end, 100000)

solution = odeint(lorenz2, population_0, time)

np.savetxt(out_file, np.column_stack((time, solution)), delimiter='\t',
           header='Time\t\tS1\t\tT1\t\tT2\t\t total S0\t\tS0 from S1\t\tS0 from T1', comments='')
