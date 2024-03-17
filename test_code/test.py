"""""
弹道计算程序
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d

# 展示高清图
from matplotlib_inline import backend_inline
backend_inline.set_matplotlib_formats('svg')

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False 

#忽略interp2d的提示
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# 导弹参数
S_ref   =   0.13
L_ref   =   5.31

# 系数
Cy_a = 7
Cy_b = 0.07

Mz_a = -0.1
Mz_b = 0.024
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
        self.dq = 0
        self.ma = 0

    # 比例导引法，给定目标位置
    def Euler1(self, before, Xm, Ym, P=0):
        dmass = P/2347
        pho, temp =  air(self.H)

        v_c=20.05*np.sqrt(temp)

        self.Time = before.Time + timestep
        self.X = before.X + before.V * np.cos(before.theta) * timestep
        self.H = before.H + before.V * np.sin(before.theta) * timestep
        self.mass = before.mass - dmass * timestep
        self.r = np.sqrt((self.X - Xm)*(self.X - Xm) + (self.H - Ym)*(self.H - Ym))
        self.dq = - before.V * np.sin(before.theta - np.arctan(( self.H - Ym)/(self.X - Xm)))/ self.r
        self.theta = before.theta + K_q * self.dq * timestep 
        self.alpha = (self.mass* before.V * K_q * self.dq + self.mass * 9.8 * np.cos(self.theta))/(P +  (Cy_a - Cy_b/ Mz_b * Mz_a) * 0.5 * pho * before.V * before.V * S_ref) /3.14159*180

        self.deltaz = - self.alpha * Mz_a / Mz_b 
        
        if self.deltaz > 30:
            self.deltaz = 30
        if self.deltaz < -30:
            self.deltaz = -30

        self.alpha = -self.deltaz / Mz_a * Mz_b

        Cx = C_X(before.ma, self.alpha)
        
        X = Cx * 0.5 * pho * before.V * before.V * S_ref
        self.V = before.V + (P*np.cos(before.alpha*3.14159625/180) - X - self.mass*9.8*np.sin(before.theta)) /self.mass*timestep
        self.ma = self.V/v_c
    
    
    
# 大气参数
def air (High):
    rho0 =1.2495
    T0  = 288.15
    Temp = T0 - 0.0065*High
    rho = rho0 * np.exp(4.25588*np.log(Temp / T0))
    return rho,Temp

def C_X(Mach, alpha):    
    # 假设你有一个数据表，这里用numpy数组表示
    # data_table 是一个二维数组，其中第一列是x值，第二列是y值
    cx0=np.array([[0.181,0.157,0.147,0.142,0.215,0.233,0.204,0.168,0.146],
        [0.175,0.151,0.141,0.136,0.209,0.226,0.195,0.158,0.135],
        [0.171,0.147,0.137,0.132,0.205,0.221,0.189,0.152,0.128],
        [0.169,0.145,0.135,0.13,0.203,0.218,0.185,0.148,0.125],
        [0.168,0.144,0.134,0.129,0.202,0.217,0.184,0.147,0.124],
        [0.169,0.145,0.135,0.13,0.203,0.218,0.185,0.148,0.125],
        [0.171,0.147,0.137,0.132,0.205,0.221,0.189,0.152,0.128],
        [0.175,0.151,0.141,0.136,0.209,0.226,0.195,0.158,0.135],
        [0.181,0.157,0.147,0.142,0.215,0.233,0.204,0.168,0.146],
        [0.189,0.165,0.155,0.15,0.223,0.244,0.218,0.184,0.162],
        [0.200 ,0.176,0.166,0.16,0.235,0.259,0.237,0.206,0.183],
        [0.213,0.19,0.18,0.174,0.249,0.278,0.262,0.236,0.207],
        [0.229,0.206,0.197,0.19,0.267,0.303,0.295,0.271,0.237],
        [0.249,0.226,0.217,0.21,0.289,0.335,0.338,0.313,0.272],
        [0.271,0.249,0.24,0.233,0.315,0.374,0.391,0.358,0.314],])
    # 提取x和y值
    x_values =np.array( [0.1, 0.3, 0.5, 0.7,1,2,3,5,8])
    y_values = np.array( [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # 创建一个插值函数
    # kind参数决定了插值的类型，'linear'代表线性插值
    interp_func = interp2d(x_values, y_values, cx0, kind='linear')
    # 现在你可以通过调用interp_func来进行插值
    # 例如，我们要在x=2.5的位置进行插值
    cx_value = interp_func(Mach, alpha)

    return cx_value.item()/1.5


# 飞行初始状态
statu_n = [statu(0, 0, 0, 0, 38/180*3.14159, 914)]
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

K_q = 3
# 第三阶段
while statu_n[-1].X < 98712 and statu_n[-1].Time < 18 :

    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],58712,36000,8e4)
    #print(statu_n[-1].V)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')
print("飞行速度：",statu_n[-1].V,'m/s')

while statu_n[-1].X < 88712 and 1 :

    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,18000, 0)
    #print(statu_n[-1].V)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')
print("飞行速度：",statu_n[-1].V,'m/s')


while statu_n[-1].X < 98712 and statu_n[-1].mass > 914 - 700:
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,16000,8e4)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')
print("飞行速度：",statu_n[-1].V,'m/s')

while statu_n[-1].X < 98712 and statu_n[-1].H > 0 :
    statu_n.append(statu(statu_n[-1].Time + timestep))
    statu_n[-1].Euler1(statu_n[-2],98712,16000,0)
print("位置",statu_n[-1].X, statu_n[-1].H)
print("飞行时间：",statu_n[-1].Time,'seconds')
print("飞行速度：",statu_n[-1].V,'m/s')

# 绘图
X_data = [n.X for n in statu_n]
H_data = [n.H for n in statu_n]
plt.plot(X_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行高度')
plt.plot([98712],[16000],'*r',markersize=12,label='目标位置')
plt.title("弹道铅垂平面轨迹")
plt.legend()  #显示上面的label
plt.xlabel('X(m)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,max(H_data)+2000)
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

plt.plot(T_data,H_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行高度m')
plt.title("飞行高度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('H(m)')#y_label
plt.ylim(0,16000)
plt.xlim(0,None)
plt.savefig('img2/飞行高度.png', dpi=300)
plt.clf()

V_data = [n.ma for n in statu_n]
plt.plot(T_data,V_data, 'r-.', alpha=0.5, linewidth=1, label='实际飞行速度V')
plt.title("飞行速度的时间变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel('速度V(m/s)')#y_label
plt.ylim(0,None)
plt.xlim(0,None)
plt.savefig('img2/飞行速度.png', dpi=300)
plt.clf()

theta_data = [n.theta/3.14159*180 for n in statu_n]
plt.plot(T_data,theta_data, 'r-.', alpha=0.5, linewidth=1, label=r'弹道倾角角$\theta$')
plt.title("飞行方案弹道倾角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'$\theta$')#y_label
#plt.ylim(-90,90)
plt.xlim(0,None)
plt.savefig('img2/飞行弹道倾角.png', dpi=300)
plt.clf()

alpha_data = [n.alpha for n in statu_n]
plt.plot(T_data,alpha_data, 'r-.', alpha=0.5, linewidth=1, label=r'攻角$\alpha$')
plt.title("飞行方案攻角")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'$\alpha$')#y_label
#plt.ylim(-90,90)
plt.xlim(0,None)
plt.savefig('img2/飞行攻角.png', dpi=300)
plt.clf()

M_data = [n.mass for n in statu_n]
plt.plot(T_data,M_data, 'r-.', alpha=0.5, linewidth=1, label=r'质量$m$')
plt.title("质量变化曲线")
plt.legend()  #显示上面的label
plt.xlabel('Time(s)') #x_label
plt.ylabel(r'm(kg)')#y_label
plt.ylim(0,917)
plt.xlim(0,None)
plt.savefig('img2/质量.png', dpi=300)
plt.clf()