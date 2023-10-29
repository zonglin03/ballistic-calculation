"""""
弹道计算程序
"""

import numpy as np
from matplotlib import pyplot as plt

# 展示高清图
from matplotlib_inline import backend_inline
backend_inline.set_matplotlib_formats('svg')

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False 

# 导弹参数
S_ref   =   0.45
L_ref   =   2.5

# 放大系数
K_phi = -0.6
K_phi_dot= 0.5 * K_phi
K_q = 5


# 仿真时间步
timestep = 0.001

# 导弹状态定义
class statu():
    __slot__=['Time','X','H','V','theta','mass','alpha','deltaz']
    
    # 位置
    # 速度
    # 欧拉角
    # 角加速度
    # 舵偏角

    # 初始化
    def __init__(self, Time, X=0, H=0, V=0, theta=0, mass=0):
        self.Time = Time
        self.X = X
        self.H = H
        self.V = V
        self.theta = theta
        self.mass = mass
        self.alpha = 0
        self.deltaz = 0
        self.q = 0

    # 显式Euler法，给定飞行高度
    def Euler(self, before, dmass):
        self.Time = before.Time + timestep
        
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep

        self.deltaz = K_phi * (self.H - High_goal(self.X)) + K_phi_dot* (before.V * np.sin(before.theta) - High_goal_dot(self.X))

        if self.deltaz > 30:
            self.deltaz = 30
        elif self.deltaz < -30:
            self.deltaz = -30

        self.alpha =  0.24 *self.deltaz

        Y = (0.25 * self.alpha + 0.05* self.deltaz) * 0.5 * air(self.H) * before.V * before.V * S_ref

        X = (0.005 * self.alpha * self.alpha + 0.2) * 0.5 * air(self.H) * before.V * before.V * S_ref
        
        self.mass = before.mass - dmass * timestep
        if dmass == 0:
            P = 0
        else:
            P = 2000
        
        self.V = before.V + (P*np.cos(self.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
        self.theta = before.theta + (P*np.sin(self.alpha*3.14159625/180) + Y - self.mass*9.8*np.cos(before.theta)) /self.mass/self.V*timestep

    # 比例导引法，给定目标位置
    def Euler2(self, before, Xm, Ym):
        self.Time = before.Time + timestep
        
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass 
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        
        self.q = np.arctan(( Ym - self.H )/(Xm - self.X))
        self.dq = - before.V * np.sin(before.theta - self.q)/ self.r

        self.theta = before.theta + K_q * self.dq * timestep
        
        P = 0
        

        self.alpha = (self.mass* before.V * K_q * self.dq + self.mass * 9.8 *np.cos(self.theta))/(P +  (0.25 + 0.05/0.24) * 0.5 * air(self.H) * before.V * before.V * L_ref) /3.14159*180

        self.deltaz = self.alpha / 0.24
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha =  0.24 *self.deltaz
        
        X = (0.005 * before.alpha * before.alpha + 0.2) * 0.5 * air(self.H) * before.V * before.V * S_ref
        self.V = before.V + (P*np.cos(self.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
    
# 大气参数
def air (High):
    rho0 =1.2495
    T0  = 288.15
    Temp = T0 - 0.0065*High
    rho = rho0 * np.exp(4.25588*np.log(Temp / T0))
    return rho

# 飞行方案
def High_goal(X):
    if X <= 9100:
        return 2000 * np.cos(0.000314 * 1.1 * X) + 5000
    elif X <= 24000:
        return 3050
    else:
        return 0
    
def High_goal_dot(X):
    if X <= 9100:
        return -2000 * 0.000314 * np.sin(0.000314 * 1.1 * X)
    elif X <= 24000:
        return 0
    else:
        return 0


# 飞行初始状态
statu_n = [statu(0, 0, 7000, 250, 0, 320)]
statu_n[0].alpha = 0
statu_n[0].deltaz = 0

Time_goal = np.arange(0,200,timestep)
X_goal = np.arange(0,24000,10)
H_goal = [High_goal(i) for i in X_goal]
plt.plot(X_goal,H_goal, 'b--', alpha=0.5, linewidth=1, label='飞行方案高度')

while statu_n[-1].X < 9100:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler(statu_n[-2],0)
    #print(statu_n[-1].alpha)

while statu_n[-1].X <= 24000:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler(statu_n[-2],0.46)
    #print(statu_n[-1].theta)


statu_n[-1].Euler2(statu_n[-2],30000,-30000)
while statu_n[-1].X <= 30000 and statu_n[-1].H > 0:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler2(statu_n[-2],30000,0)
    #print(statu_n[-1].V)


X_data = [n.X for n in statu_n]
H_data = [n.H for n in statu_n]
plt.plot(X_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行高度')

# 飞行高度绘图
plt.legend()  #显示上面的label
plt.xlabel('X') #x_label
plt.ylabel('H')#y_label
plt.ylim(0,8000)
plt.xlim(0,30000) #仅设置y轴坐标范围


plt.figure(2)
T_data = [n.Time for n in statu_n]
deltaz_data = [n.deltaz for n in statu_n]
plt.plot(T_data,deltaz_data, 'r-.', alpha=0.5, linewidth=1, label='舵偏角$\delta z$')
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
"""
plt.figure(3)
M_data = [n.mass for n in statu_n]
plt.plot(T_data,M_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('速度V')#y_label
plt.ylim(250,350)
plt.xlim(0,200)
"""
plt.figure(4)
V_data = [n.V for n in statu_n]
plt.plot(T_data,V_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('速度V')#y_label
plt.ylim(100,250)
plt.xlim(0,200)

plt.show()
