import numpy as np
from scipy.interpolate import interp2d,CloughTocher2DInterpolator

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

# 假设你有一个数据表，这里用numpy数组表示
# data_table 是一个二维数组，其中第一列是x值，第二列是y值
cx0=np.array([[0.253, 0.25, 0.247, 0.245, 0.319, 0.309, 0.257, 0.184, 0.146],
        [0.247, 0.244, 0.241, 0.239, 0.313, 0.302, 0.248, 0.174, 0.135],
        [0.243, 0.24, 0.237, 0.236, 0.309, 0.297, 0.242, 0.168, 0.128],
        [0.241, 0.238, 0.235, 0.234, 0.307, 0.295, 0.239, 0.164, 0.125],
        [0.24, 0.238, 0.235, 0.233, 0.307, 0.294, 0.238, 0.163, 0.124],
        [0.241, 0.238, 0.235, 0.234, 0.307, 0.295, 0.239, 0.164, 0.125],
        [0.243, 0.24, 0.237, 0.236, 0.309, 0.297, 0.242, 0.168, 0.128],
        [0.247, 0.244, 0.241, 0.239, 0.313, 0.302, 0.248, 0.174, 0.135],
        [0.253, 0.25, 0.247, 0.245, 0.319, 0.309, 0.257, 0.184, 0.146],
        [0.261, 0.258, 0.255, 0.253, 0.327, 0.32, 0.271, 0.2, 0.162],
        [0.271, 0.269, 0.266, 0.263, 0.338, 0.334, 0.29, 0.222, 0.183],
        [0.284, 0.282, 0.279, 0.276, 0.352, 0.353, 0.314, 0.252, 0.207],
        [0.3, 0.298, 0.295, 0.292, 0.369, 0.378, 0.347, 0.286, 0.237],
        [0.319, 0.317, 0.315, 0.311, 0.391, 0.409, 0.389, 0.328, 0.272],
        [0.341, 0.34, 0.337, 0.334, 0.416, 0.447, 0.442, 0.373, 0.314]])
# 提取x和y值
x_values =np.array( [0.1, 0.3, 0.5, 0.7,1,2,3,5,8])
y_values = np.array( [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
# 创建一个插值函数
# kind参数决定了插值的类型，'linear'代表线性插值
interp_func = interp2d(x_values, y_values, cx0, kind='linear')
# 现在你可以通过调用interp_func来进行插值
# 例如，我们要在x=2.5的位置进行插值
interp_value = interp_func(0.2,2)
print(interp_value)  # 这将输出在x=2.5时的y值估计