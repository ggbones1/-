# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:28:17 2023

@author: 24397
"""

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import numpy as np

# 输入变量 X 和输出变量 y：
X = np.array([[0.1, 0.2, 0.3, 0.4],
              [0.5, 0.6, 0.7, 0.8],
              [0.9, 1.0, 1.1, 1.2],
              [1.3, 1.4, 1.5, 1.6]])

y = np.array([[0.05,-0.15],
              [-0.25,-0.45],
              [-0.65,-0.85],
              [-1,-35]])

# 核函数：
kernel = RBF(length_scale=1e-2)+ WhiteKernel(noise_level=1)

# 高斯过程回归模型：
gp = GaussianProcessRegressor(kernel=kernel, alpha=0.1,
                              n_restarts_optimizer=10, normalize_y=True)

# 拟合数据：
gp.fit(X, y)

# 在新的输入值处估计输出值并计算置信区间：
X_new = np.array([[1.7, 1.8, 1.9, 2.0],
                  [2.1, 2.2, 2.3, 2.4]])

y_pred, sigma = gp.predict(X_new, return_std=True)

lower_bound = y_pred - 1.96 * sigma
upper_bound = y_pred + 1.96 * sigma

print("Predicted Output: ", y_pred)
print("Confidence Interval Lower Bound: ", lower_bound)
print("Confidence Interval Upper Bound: ", upper_bound)