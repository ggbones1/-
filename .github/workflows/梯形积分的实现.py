# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 09:40:24 2023

@author: 24397
"""

import pandas as pd
import numpy as np
from scipy.integrate import trapz

# 读取Excel表格中的'A'列（第一列）数据到数组里
dx = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\计算公式.xlsx', usecols=[0])
dy = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\计算公式.xlsx', usecols=[4])
merged_df = pd.concat([dx, dy], axis=1, ignore_index=True)
merged_array = merged_df.to_numpy()
y = merged_array[:, 1]
x = merged_array[:, 0]

SM_up = trapz(y=y, x=x)
#%%
import pandas as pd
import numpy as np
from scipy.integrate import trapz

# 读取Excel表格中的'A'列（第一列）数据到数组里
dx = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\计算公式.xlsx', usecols=[0])
dy = pd.read_excel(r'F:\ADM02_横截面_稳态数据\10Dh\1mps_连续\计算公式.xlsx', usecols=[5])
merged_df = pd.concat([dx, dy], axis=1, ignore_index=True)
merged_array = merged_df.to_numpy()
y = merged_array[:, 1]
x = merged_array[:, 0]

integral_value = trapz(y=y, x=x)
SM_down=integral_value*8.32
SM=SM_up/SM_down
print("涡流因子为",SM)
#%%

