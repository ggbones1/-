# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 20:39:30 2023

@author: 24397
"""

import pandas as pd
import numpy as np
# 读取Excel文件
df = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\u 子通道.xls', sheet_name='Sheet1')
# 获取指定位置的5*5的数据
data = df.iloc[0, 2:27]
data1 = df.iloc[1:26,1]
x = np.matrix(data) #x轴的相对坐标
y=np.matrix(data1)  #y轴的相对坐标
r_square=np.zeros((25,25)) #r的平方的数组
for i in range(25):
    r_square[i] = np.square(x)+np.square(y)  # 逐行计算结果数组
r=np.sqrt(r_square)
print(r_square,r)

#%%
import pandas as pd
import numpy as np
# 读取Excel文件
df = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\v 子通道.xls', sheet_name='Sheet1')
# 获取指定位置的5*5的数据
data2 = df.iloc[0, 2:27]
data3 = df.iloc[1:26,1]
x = np.matrix(data2) #x轴的相对坐标
y=np.matrix(data3)  #y轴的相对坐标
r_square=np.zeros((25,25)) #r的平方的数组
for i in range(25):
    r_square[i] = np.square(x)+np.square(y)  # 逐行计算结果数组
r=np.sqrt(r_square)
df = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\v 子通道.xls', sheet_name='Sheet1')
# 获取指定位置的5*5的数据
data4 = df.iloc[1:27, 2:27]
v=np.matrix(data4)
print(v)

#%%
import pandas as pd
import numpy as np
# 读取Excel文件
df = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\u 子通道.xls', sheet_name='Sheet1')
# 获取指定位置的5*5的数据
data5 = df.iloc[1:27, 2:27]
u=np.matrix(data5)
print(u)

#%%
import numpy as np
import pandas as pd


# 将三个数组都转化为625*1的数组
new_r = r.reshape((625, 1))
new_v = v.reshape((625, 1))
new_u = u.reshape((625, 1))
new_r_square=r_square.reshape((625, 1))
# 将数据转换成pandas数据框
output  = np.hstack((new_r, new_r_square,new_u,new_v))
df = pd.DataFrame(output)
df.columns = ['r', 'r_square', 'u', 'v']
df.to_csv('F:\\ADM02_横截面_稳态数据\\10Dh\\1mps_连续\\use in sum.csv', index=False)



