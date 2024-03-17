import numpy as np
# 读取文档中的数据
with open('data/cx.txt', 'r', encoding='UTF-8') as f:
    lines = f.readlines()

# 提取关键数据
data = []
for line in lines:
    items = line.strip().split()
    if items:
        data.append(items)

# 数据预处理，提取马赫数、轴向力系数等数据
high_numbers = np.array(data[1], dtype=float)
mach_numbers = np.array(data[3], dtype=float)
cx_values = [row for row in np.array(data[5:], dtype=float)]

# 输出结果
print("海拔(km)：", high_numbers)
print("马赫数：", mach_numbers)
print("轴向力系数 C_x：")
for i, row in enumerate(cx_values):
    print("海拔{}km: {}".format(mach_numbers[i], row))
2