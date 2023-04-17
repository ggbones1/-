# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 17:39:44 2023

@author: 24397
"""

import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'SimHei' # 设置字体为SimHei

#取数据文件
data = pd.read_excel(r'C:\Users\24397\Documents\WPSDrive\621846420\WPS云盘\平均网格质量.xlsx', sheet_name='Sheet1')

#出第一列和第二列数据
x = data.iloc[:, 0]
y = data.iloc[:, 1]

#图
plt.figure(figsize=(8, 6)) # 设置图形大小
plt.plot(x, y, 'o-')
plt.grid() # 显示网格线

#置图表标题、x轴标题和y轴标题
plt.title('网格无关性分析')
plt.xlabel('网格数量')
plt.ylabel('平均网格质量')

#示图表
plt.show()