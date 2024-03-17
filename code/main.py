"""
    一级半运载火箭弹道计算
"""


###############  环境准备  ###############
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d, interp1d,CloughTocher2DInterpolator

# 展示高清图
from matplotlib_inline import backend_inline
backend_inline.set_matplotlib_formats('svg')

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False 

# 忽略提示
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


###############  定义数值  ###############
# 导弹参数
S_ref   =   19.6

## 质量
# 构型：芯一级＋2个大助推器＋2个小助推器
m_0 = 645000 # 总质量
m0 = 6600 # 有效载荷质量

# 芯一级
m1_0 = 173000   # 芯一级起飞质量
m1_1 = 19500    # 芯一级结构质量 
qm_1 = 334    # 芯一级秒流量

# 大助推器
m2_0 = 158000   
m2_1 = 14500    
qm_2 = 820

# 小助推器(单个助推器)：
m3_0 = 72700
m3_1 = 7100
qm_3 = 410

# 整流罩
m4 = 4000

# 经纬度，单位度
latitude_0 = 19.6  # 维度
longitude_0 = 111   # 经度

# 攻角参数
alpha_max = 0.65 #单位弧度
alpha_costant = 0.1 # 转弯段参数

# 时间
t1 = 15 # 垂直起飞段结束，进入转弯段
t2 = 120    # 120s整流罩分离。
t3 = 160    # 160秒两小助推关机分离，大助推外侧关机（推力减半，秒流量减半）
t4 = 190 # 190秒两大助推内侧关机分离
t5 = 459.58 # 芯一级关机

# 目标轨道
Target_orbital_altitude = 500e3
Target_orbital_inclination = 60 #单位 度
#轨道高度偏差1km，轨道倾角偏差0.5°

# 放大系数
K_phi = 1
K_alpha = 1

# 仿真时间步
timestep = 0.01
###############  类的定义  ###############
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
        
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
        self.theta = before.theta + (P*np.sin(self.alpha*3.14159625/180) + Y - self.mass*9.8*np.cos(before.theta)) /self.mass/self.V*timestep

    # 比例导引法，给定目标位置
    def Euler2(self, before, Xm, Ym):
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


###############  部分函数  ###############
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

# 读入数据
def read_file(file_path):
    # 读取文档中的数据
    with open(file_path, 'r', encoding='UTF-8') as f:
        lines = f.readlines()

    # 提取关键数据
    data = []
    for line in lines:
        items = line.strip().split()
        if items:
            data.append(items)
    return data

def alpha_cornering( time ):
    b = np.exp(alpha_constant*(t1-time))
    return 4* alpha_max * b *(b-1)


###############  插值函数  ###############
# cx插值数据预处理
data = read_file('data/cx.txt')
cx_high_numbers = np.array(data[1], dtype=float)
cx_mach_numbers = np.array(data[3], dtype=float)
cx_values = [row for row in np.array(data[5:], dtype=float)]
# kind参数决定了插值的类型，'linear'代表线性插值
cx_interp_func = interp2d(cx_high_numbers, cx_mach_numbers, cx_values, kind='linear')
# print(cx_interp_func(0.1,0.1),cn_alpha_interp_func(0.3))

# cx插值数据预处理
data = read_file('data/cn_alpha.txt')
cn_alpha_mach = np.array(data[1], dtype=float)
cn_alpha_values = [row for row in np.array(data[3:], dtype=float)]
# 允许外插
cn_alpha_interp_func = interp1d(cn_alpha_mach, cn_alpha_values, kind='linear', fill_value='extrapolate')

## 芯一级推力P_1 插值预处理
data = read_file('data/P_1.txt')
P_1_high = np.array(data[1], dtype=float)
P_1_values = [row for row in np.array(data[3:], dtype=float)]

P_1_interp_func = interp1d(P_1_high, P_1_values, kind='linear', fill_value='extrapolate')

## 两大助推器总推力P_2 插值预处理
data = read_file('data/P_2.txt')
P_2_high = np.array(data[1], dtype=float)
P_2_values = [row for row in np.array(data[3:], dtype=float)]

P_2_interp_func = interp1d(P_2_high, P_2_values, kind='linear', fill_value='extrapolate')
print(P_2_interp_func(20))

## 两小助推器总推力P_3 插值预处理
data = read_file('data/P_1.txt')
P_3_high = np.array(data[1], dtype=float)
P_3_values = [row for row in np.array(data[3:], dtype=float)]

P_3_interp_func = interp1d(P_3_high, P_3_values, kind='linear', fill_value='extrapolate')
print(P_3_interp_func(20))

###############  计算初始化  ###############
# 发射方位角
A_0 = 

# 飞行初始状态
statu_n = [statu(0, 0, 7000, 250, 0, 320)]
statu_n[0].alpha = 0
statu_n[0].deltaz = 0

X_goal = np.arange(0,24000,10)
H_goal = [High_goal(i) for i in X_goal]
plt.plot(X_goal,H_goal, 'b--', alpha=0.5, linewidth=1, label='飞行方案高度')

# 方位角
###############  计算开始  ###############
# 第一阶段
while statu_n[-1].Time < 15:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler(statu_n[-2],0)


###############  绘图可视化  ###############
# 绘图
X_data = [n.X for n in statu_n]
H_data = [n.H for n in statu_n]
plt.plot(X_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行高度')
plt.title("弹道铅垂平面轨迹")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,8000)
plt.xlim(0,30000) #仅设置y轴坐标范围
plt.savefig('img/飞行轨迹.png', dpi=300)
plt.clf()

T_data = [n.Time for n in statu_n]
deltaz_data = [n.deltaz for n in statu_n]
plt.plot(T_data,deltaz_data, 'r-.', alpha=0.5, linewidth=1, label='舵偏角$\delta z$')
plt.title("飞行方案舵偏角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('$\delta z$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
plt.savefig('img/飞行舵偏角.png', dpi=300)
plt.clf()

plt.plot(T_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.title("飞行高度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,7000)
plt.xlim(0,200)
plt.savefig('img/飞行高度.png', dpi=300)
plt.clf()

V_data = [n.V for n in statu_n]
plt.plot(T_data,V_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.title("飞行速度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('速度V')#y_label
plt.ylim(100,250)
plt.xlim(0,200)
plt.savefig('img/飞行速度.png', dpi=300)
plt.clf()

theta_data = [n.theta*180/3.14159 for n in statu_n]
plt.plot(T_data,theta_data, 'r-.', alpha=0.5, linewidth=1, label=r'舵偏角$\theta$')
plt.title("飞行方案弹道倾角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'$\theta$')#y_label
plt.ylim(-50,50)
plt.xlim(0,200)
plt.savefig('img/飞行弹道倾角.png', dpi=300)
plt.clf()