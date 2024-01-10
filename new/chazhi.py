import numpy as np
from scipy.interpolate import interp1d
# 假设你有一个数据表，这里用numpy数组表示
# data_table 是一个二维数组，其中第一列是x值，第二列是y值
data_table = np.array(
    
)
# 提取x和y值
x_values = data_table[:, 0]
y_values = data_table[:, 1]
# 创建一个插值函数
# kind参数决定了插值的类型，'linear'代表线性插值
interp_func = interp1d(x_values, y_values, kind='linear')
# 现在你可以通过调用interp_func来进行插值
# 例如，我们要在x=2.5的位置进行插值
interp_value = interp_func(2.5)
print(interp_value)  # 这将输出在x=2.5时的y值估计