"""
弹道计算程序
分析放大系数的影响
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
        self.q = 0

    # 显式Euler法，给定飞行高度
    def Euler(self, before, dmass,K_phi, K_phi_dot):
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
        
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
        self.theta = before.theta + (P*np.sin(self.alpha*3.14159625/180) + Y - self.mass*9.8*np.cos(before.theta)) /self.mass/self.V*timestep

    # 比例导引法，给定目标位置
    def Euler2(self, before, Xm, Ym, K_q):
        self.Time = before.Time + timestep
        
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass 
        
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        
        self.dq = - before.V * np.sin(before.theta - np.arctan(( self.H - Ym)/(self.X - Xm)))/ self.r

        self.theta = before.theta + K_q * self.dq * timestep
        
        P = 0
        

        self.alpha = (self.mass* before.V * K_q * self.dq + self.mass * 9.8 * np.cos(self.theta))/(P +  (0.25 + 0.05/0.24) * 0.5 * air(self.H) * before.V * before.V * S_ref) /3.14159*180

        self.deltaz = self.alpha / 0.24
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha =  0.24 *self.deltaz
        
        X = (0.005 * self.alpha * self.alpha + 0.2) * 0.5 * air(self.H) * before.V * before.V * S_ref
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
    
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

# 飞行方案的时间导数
def High_goal_dot(X):
    if X <= 9100:
        return -2000 * 0.000314 * np.sin(0.000314 * 1.1 * X)
    elif X <= 24000:
        return 0
    else:
        return 0

def calculate(K_phi ,K_phi_dot, K_q):

    # 飞行初始状态
    statu_n = [statu(0, 0, 7000, 250, 0, 320)]
    statu_n[0].alpha = 0
    statu_n[0].deltaz = 0

    # 第一阶段
    while statu_n[-1].X < 9100:
        statu_n.append(statu(statu_n[-1].Time + timestep))
        statu_n[-1].Euler(statu_n[-2],0,K_phi ,K_phi_dot )

    # 第二阶段
    while statu_n[-1].X <= 24000:
        statu_n.append(statu(statu_n[-1].Time + timestep))
        statu_n[-1].Euler(statu_n[-2],0.46,K_phi ,K_phi_dot)


    # 第三阶段
    while statu_n[-1].X <= 30000 and statu_n[-1].H > 0:
        statu_n.append(statu(statu_n[-1].Time + timestep))
        statu_n[-1].Euler2(statu_n[-2],30000,0,K_q)

    return statu_n

# 飞行方案
X_goal = np.arange(0,24000,10)
H_goal = [High_goal(i) for i in X_goal]

# 放大系数
statu_1=calculate(-0.2,-0.5,2)
statu_2=calculate(-0.8,-0.5,2)

# 绘图

## 第一个放大系数
X_data_1 = [n.X for n in statu_1]
H_data_1 = [n.H for n in statu_1]
X_data_2 = [n.X for n in statu_2]
H_data_2 = [n.H for n in statu_2]

plt.plot(X_goal,H_goal, 'c--', alpha=0.5, linewidth=1, label='飞行方案高度')

plt.plot(X_data_1,H_data_1, 'r-.', alpha=0.5, linewidth=1, label=r'$k_\varphi=-0.2$')
plt.plot(X_data_2,H_data_2, 'b-.', alpha=0.5, linewidth=1, label=r'$k_\varphi=-0.8$')

plt.title("弹道铅垂平面轨迹第一阶段对比")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(3000,7000)
plt.xlim(0,10000) #仅设置y轴坐标范围
plt.savefig('img/飞行轨迹2.png', dpi=300)
plt.clf()


T_data_1 = [n.Time for n in statu_1]
T_data_2 = [n.Time for n in statu_2]
deltaz_data_1 = [n.deltaz for n in statu_1]
deltaz_data_2 = [n.deltaz for n in statu_2]

plt.plot(T_data_1,deltaz_data_1, 'r-.', alpha=0.5, linewidth=1, label=r'$k_\varphi=-0.2$')
plt.plot(T_data_2,deltaz_data_2, 'b-.', alpha=0.5, linewidth=1, label=r'$k_\varphi=-0.8$')

plt.title("飞行方案舵偏角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
plt.savefig('img/飞行舵偏角2.png', dpi=300)
plt.clf()

## 第二个放大系数
statu_3 = calculate(-0.6,-0.3,3)
statu_4 = calculate(-0.6,-0.5,3)

X_data_3 = [n.X for n in statu_3]
H_data_3 = [n.H for n in statu_3]
X_data_4 = [n.X for n in statu_4]
H_data_4 = [n.H for n in statu_4]

plt.plot(X_goal,H_goal, 'c--', alpha=0.5, linewidth=1, label='飞行方案高度')

plt.plot(X_data_3,H_data_3, 'r-.', alpha=0.5, linewidth=1, label=r'$\dot{k}_\varphi=-0.3$')
plt.plot(X_data_4,H_data_4, 'b-.', alpha=0.5, linewidth=1, label=r'$\dot{k}_\varphi=-0.5$')

plt.title("弹道铅垂平面轨迹阶跃处对比")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(2800,3200)
plt.xlim(8500,11000) #仅设置y轴坐标范围
plt.savefig('img/飞行轨迹3.png', dpi=300)
plt.clf()

T_data_3 = [n.Time for n in statu_3]
T_data_4 = [n.Time for n in statu_4]
deltaz_data_3 = [n.deltaz for n in statu_3]
deltaz_data_4 = [n.deltaz for n in statu_4]

plt.plot(T_data_3,deltaz_data_3, 'r-.', alpha=0.5, linewidth=1, label=r'$\dot{k}_\varphi=-0.3$')
plt.plot(T_data_4,deltaz_data_4, 'b-.', alpha=0.5, linewidth=1, label=r'$\dot{k}_\varphi=-0.5$')

plt.title("飞行方案舵偏角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
plt.savefig('img/飞行舵偏角3.png', dpi=300)
plt.clf()


## 第三个放大系数
statu_5 = calculate(-0.6, -0.5, 3)
statu_6 = calculate(-0.6, -0.5, 6)

X_data_5 = [n.X for n in statu_5]
H_data_5 = [n.H for n in statu_5]
X_data_6 = [n.X for n in statu_6]
H_data_6 = [n.H for n in statu_6]

plt.plot(X_goal,H_goal, 'c--', alpha=0.5, linewidth=1, label='飞行方案高度')

plt.plot(X_data_5,H_data_5, 'r-.', alpha=0.5, linewidth=1, label=r'$k_3=3$')
plt.plot(X_data_6,H_data_6, 'b-.', alpha=0.5, linewidth=1, label=r'$k_3=5$')


plt.title("弹道铅垂平面第三阶段轨迹对比")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,3500)
plt.xlim(24000,30000) #仅设置y轴坐标范围
plt.savefig('img/飞行轨迹4.png', dpi=300)
plt.clf()

T_data_5 = [n.Time for n in statu_5]
T_data_6 = [n.Time for n in statu_6]
deltaz_data_5 = [n.deltaz for n in statu_5]
deltaz_data_6 = [n.deltaz for n in statu_6]

plt.plot(T_data_5,deltaz_data_5, 'r-.', alpha=0.5, linewidth=1, label=r'$k_3=3$')
plt.plot(T_data_6,deltaz_data_6, 'b-.', alpha=0.5, linewidth=1, label=r'$k_3=5$')
plt.title("飞行方案舵偏角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
plt.savefig('img/飞行舵偏角4.png', dpi=300)
plt.clf()
