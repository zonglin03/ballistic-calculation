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
S_ref   =   0.13
L_ref   =   5.31

# 放大系数
K_phi = - 0.5
K_phi_dot= 0.6* K_phi
K_q = 9

# 系数
Cy_a = 0.25
Cy_b = 0.05

Cx = 0.25
Cx_a = 0.05

Mz_a = 0.1
Mz_b = 0.024
# 仿真时间步
timestep = 0.01

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
        self.dq = 0

    # 比例导引法，给定目标位置
    def Euler1(self, before, Xm, Ym, P=0):
        dmass = P/2400
        self.Time = before.Time + timestep
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass - dmass * timestep
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        self.dq = - before.V * np.sin(before.theta - np.arctan(( self.H - Ym)/(self.X - Xm)))/ self.r
        self.theta = before.theta + K_q * self.dq * timestep
        self.alpha = (self.mass* before.V * K_q * self.dq + self.mass * 9.8 * np.cos(self.theta))/(P +  (Cy_a - Cy_b/ Mz_b * Mz_a) * 0.5 * air(self.H) * before.V * before.V * S_ref) /3.14159*180

        self.deltaz = - self.alpha * Mz_a / Mz_b 
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha = -self.deltaz / Mz_a * Mz_b
        
        X = (Cx_a * self.alpha * self.alpha + Cx) * 0.5 * air(self.H) * before.V * before.V * S_ref
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
    
    # 比例导引法，给定目标位置
    def Euler2(self, before, Xm, Ym):
        P = 0
        dmass=P/2400
        self.Time = before.Time + timestep
        
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass - dmass * timestep
        
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        
        self.dq = - before.V * np.sin(before.theta - np.arctan(( self.H - Ym)/(self.X - Xm)))/ self.r

        self.theta = before.theta + K_q * self.dq * timestep
        
        
        

        self.alpha = (self.mass* before.V * K_q * self.dq + self.mass * 9.8 * np.cos(self.theta))/(P +  (0.25 + 0.05/0.24) * 0.5 * air(self.H) * before.V * before.V * S_ref) /3.14159*180

        self.deltaz = self.alpha / 0.24
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha =  0.24 *self.deltaz
        
        X = (0.005 * self.alpha * self.alpha + 0.2) * 0.5 * air(self.H) * before.V * before.V * S_ref
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep

    def Euler3(self, before, Xm, Ym):
        self.Time = before.Time + timestep
        
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass 
        
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        
        self.dtheta = - K_q * before.V * np.sin(before.theta - np.arctan(( self.H - Ym)/(self.X - Xm)))/ self.r

        self.theta = before.theta + self.dtheta * timestep
        
        P = 0

        self.alpha = (self.mass* before.V * self.dtheta + self.mass * 9.8 * np.cos(self.theta))/(P +  (0.25 + 0.05/0.24) * 0.5 * air(self.H) * before.V * before.V * S_ref) /3.14159*180

        self.deltaz = self.alpha / 0.24
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha =  0.24 * self.deltaz
        
        X = (0.005 * self.alpha * self.alpha + 0.2) * 0.5 * air(self.H) * before.V * before.V * S_ref
        self.V = before.V + ((P*np.cos(self.alpha*3.14159/180) - X)/self.mass-9.8*np.sin(self.theta)) *timestep
    
# 大气参数
def air (High):
    rho0 =1.2495
    T0  = 288.15
    Temp = T0 - 0.0065*High
    rho = rho0 * np.exp(4.25588*np.log(Temp / T0))
    return rho


# 飞行初始状态
statu_n = [statu(0, 0, 0, 0, 38, 914)]
statu_n[0].alpha = 0
statu_n[0].deltaz = 0

# 第一阶段
while statu_n[-1].mass > 914-600 and statu_n[-1].H < 25000 and 0:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],30000,30000,30000)
    print(statu_n[-1].X,statu_n[-1].H)

print(statu_n[-1].Time)

# 第二阶段
while statu_n[-1].X < 95000 and statu_n[-1].Time < 150 and 0:
    K_q=8
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,16000,0)
    #print(statu_n[-1].V)
print(statu_n[-1].Time)

# 第三阶段
while statu_n[-1].X < 98712 and statu_n[-1].Time < 23 :
    K_q=15
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,16000,8e4)
    #print(statu_n[-1].V)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')

while statu_n[-1].X < 98712 and statu_n[-1].Time < 90 :
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,16000,0)
    x = [statu_n[-1].X, statu_n[-1].H]
    #print(statu_n[-1].V)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')

# 绘图
X_data = [n.X for n in statu_n]
H_data = [n.H for n in statu_n]
plt.plot(X_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行高度')
plt.title("弹道铅垂平面轨迹")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,None)
plt.xlim(0,100000) #仅设置y轴坐标范围
plt.savefig('img2/飞行轨迹.png', dpi=300)
plt.clf()

T_data = [n.Time for n in statu_n]
deltaz_data = [n.deltaz for n in statu_n]
plt.plot(T_data,deltaz_data, 'r-.', alpha=0.5, linewidth=1, label='舵偏角$\delta z$')
plt.title("飞行方案舵偏角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,None)
plt.savefig('img2/飞行舵偏角.png', dpi=300)
plt.clf()

plt.plot(T_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.title("飞行高度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,16000)
plt.xlim(0,None)
plt.savefig('img2/飞行高度.png', dpi=300)
plt.clf()

V_data = [n.V for n in statu_n]
plt.plot(T_data,V_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.title("飞行速度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('速度V')#y_label
plt.ylim(0,None)
plt.xlim(0,None)
plt.savefig('img2/飞行速度.png', dpi=300)
plt.clf()

theta_data = [n.theta*180/3.14159 for n in statu_n]
plt.plot(T_data,theta_data, 'r-.', alpha=0.5, linewidth=1, label=r'舵偏角$\theta$')
plt.title("飞行方案弹道倾角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'$\theta$')#y_label
plt.ylim(-50,50)
plt.xlim(0,None)
plt.savefig('img2/飞行弹道倾角.png', dpi=300)
plt.clf()

M_data = [n.mass for n in statu_n]
plt.plot(T_data,M_data, 'r-.', alpha=0.5, linewidth=1, label=r'舵偏角$\theta$')
plt.title("质量变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'$\theta$')#y_label
plt.ylim(0,917)
plt.xlim(0,None)
plt.savefig('img2/质量.png', dpi=300)
plt.clf()