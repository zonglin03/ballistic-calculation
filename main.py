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


# 导弹状态定义
class statu():
    time = 0
    mass = 0
    # 位置
    x = 0
    y = 0
    z = 0
    # 速度
    vx = 0
    vy = 0
    vz = 0
    # 欧拉角
    theta = 0
    phi = 0
    alpha = 0
    # 角加速度
    theta_a =0
    # 舵偏角
    


# 大气参数
def air (High):
    rho0 =1.2495
    T0  = 288.15
    Temp = T0 - 0.0065*High
    rho = rho0 * np.exp(4.25588*np.log(Temp / T0))

    return (Temp , rho)

# 升力系数

