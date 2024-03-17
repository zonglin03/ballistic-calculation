import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 展示高清图
from matplotlib_inline import backend_inline
backend_inline.set_matplotlib_formats('svg')

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False 

# 发动机参数
initial_volume = 25e-6  # 初始自由容积，m^3
throat_diameter = 7.75e-3  # 喷管喉径，m
propellant_density = 1660  # 推进剂密度，kg/m^3
characteristic_velocity = 1525  # 特征速度，m/s
burning_rate_coefficient = 4.93  # 燃速系数，P^0.35 mm/s (P压强单位kgf/cm2)
k = 1.2  # 燃气参数
ignition_pressure = 2e6  # 点火压强，Pa

# 管状药信息
num_grains = 9
outer_radius = 5.3e-3  # 外径，m
inner_radius = 3.7e-3  # 内径，m
grain_length = 60e-3  # 长度，m

# 计算燃烧室条件
burning_surface_area = num_grains * np.pi * grain_length * (2 * outer_radius)  # 燃烧表面积，m^2
initial_pressure = ignition_pressure  # 初始压强等于点火压强
initial_density = propellant_density
initial_temperature = (initial_pressure / (initial_density * 287.05))  # 气体常数R=287.05 J/(kg·K)
initial_temperature = initial_temperature * initial_volume / burning_surface_area  # 初始温度，K

# 计算喷嘴出口条件
exit_pressure = initial_pressure * (1 + (k - 1) / 2) ** (k / (k - 1))  # 喷嘴出口压强，Pa
exit_temperature = initial_temperature * (1 + (k - 1) / 2)  # 喷嘴出口温度，K
exit_velocity = np.sqrt(2 * k * 287.05 * initial_temperature * (1 - (initial_pressure / exit_pressure) ** ((k - 1) / k)))  # 喷嘴出口速度，m/s
# 记录时间和燃燃室压强
times = []  # 存储时间
pressures = []  # 存储燃燃室压强

# 使用R-K方法计算内弹道
def rocket_equations(state, t):
    chamber_pressure, throat_area_ratio = state
    exit_area_ratio = (1 / throat_area_ratio)  # 喷嘴出口面积比

    # 使用R-K方程
    mdot = burning_surface_area * burning_rate_coefficient * chamber_pressure ** 0.35  # 质量流率
    throat_area = (mdot / (initial_density * exit_velocity))  # 喉面积
    throat_area_ratio = throat_area / throat_diameter ** 2  # 喉面积比

    # 记录时间和燃燃室压强
    times.append(t)
    pressures.append(chamber_pressure)

    return [
        mdot - propellant_density * burning_surface_area * throat_area * throat_area_ratio ** 0.5,
        exit_area_ratio - throat_area_ratio * (1 + (k - 1) / 2) ** (1 / (1 - k))
    ]

# 初始猜测值
initial_guess = [initial_pressure, 1.0]

# 时间步长和总时间
dt = 0.001  # 时间步长，秒
total_time = 10.0  # 总时间，秒

# 初始化时间和初始状态
t = 0
state = initial_guess

while t < total_time:
    # 使用fsolve函数求解非线性方程组，传递当前时间作为参数
    state = fsolve(rocket_equations, state, args=(t,), maxfev=1000)
    t += dt


# 转换时间和压强列表为NumPy数组
times = np.array(times)
pressures = np.array(pressures)

# 绘制压强与时间的关系图
plt.plot(times, pressures / 1e6)  # 压强以MPa为单位
plt.xlabel('时间 (秒)')
plt.ylabel('燃烧室压强 (MPa)')
plt.title('燃烧室压强与时间的关系')
# plt.grid(True)
plt.show()