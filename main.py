import torch
import numpy as np


# 导弹参数
M0  =   320
P   =   2000
S_lef   =   0.45
L_ref   =   2.5

# 初始运动状态
V0  =   250
x0  =   0
H0  =   7000
theta0  =   0
phi0    =   0
alpha0  =   0
    

# 大气参数
def air (High):
    rho0 =1.2495
    T0  = 288.15
    Temp = T0 - 0.0065*High
    rho = rho0 * np.exp(4.25588*np.log(Temp / T0))

    return (Temp , rho)

# 升力系数

