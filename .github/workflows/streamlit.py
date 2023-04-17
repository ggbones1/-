# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 17:35:29 2023

"""
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
# 设置页面标题和描述
st.title('BP神经网络')
st.write('这是一个四输入两输出的神经网络示例，您可以在页面中调整参数并进行训练和预测。')

# 设置输入和输出数据
X = np.random.rand(100, 4)
y = np.random.rand(100, 2)

# 设置默认参数
hidden_layer_sizes = (5, 2)
activation = 'relu'
solver = 'adam'
learning_rate_init = 0.001
max_iter = 200
alpha = 0.0001

# 创建神经网络对象
mlp = MLPRegressor(hidden_layer_sizes=hidden_layer_sizes, activation=activation,
                   solver=solver, learning_rate_init=learning_rate_init,
                   max_iter=max_iter, alpha=alpha)

# 创建页面上的参数调整器
st.sidebar.title('参数调整')
st.sidebar.write('使用这些滑块来调整神经网络的参数。')
hl1 = st.sidebar.slider('第一层隐藏层节点数', 1, 10, 5)
hl2 = st.sidebar.slider('第二层隐藏层节点数', 1, 10, 2)
hidden_layer_sizes = (hl1, hl2)
activation = st.sidebar.selectbox('激活函数', ['identity', 'logistic', 'tanh', 'relu'])
solver = st.sidebar.selectbox('优化器', ['lbfgs', 'sgd', 'adam'])
learning_rate_init = st.sidebar.slider('学习率', 0.0001, 0.01, 0.001, 0.0001)
max_iter = st.sidebar.slider('最大迭代次数', 100, 1000, 200, 50)
alpha = st.sidebar.slider('正则化系数', 0.00001, 0.001, 0.0001, 0.00001)

# 训练神经网络
mlp.set_params(hidden_layer_sizes=hidden_layer_sizes, activation=activation,
               solver=solver, learning_rate_init=learning_rate_init,
               max_iter=max_iter, alpha=alpha)
mlp.fit(X, y)

# 进行预测
y_pred = mlp.predict(X)

# 绘制图像
fig, ax = plt.subplots()
ax.plot(y, label='真实值')
ax.plot(y_pred, label='预测值')
ax.legend()
st.pyplot(fig)