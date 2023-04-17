# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 12:44:07 2023

@author: 24397
"""

import numpy as np
import matplotlib.pyplot as plt

# 定义sigmoid函数
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

class FourLayerBP:
    def __init__(self, input_size, hidden_sizes, output_size):
        self.input_size = input_size
        self.hidden_sizes = hidden_sizes
        self.output_size = output_size

        # 初始化权重矩阵
        self.IH_weights = np.random.uniform(low=-0.5, high=0.5, size=(input_size, hidden_sizes[0]))
        self.HH_weights = np.random.uniform(low=-0.5, high=0.5, size=(hidden_sizes[0], hidden_sizes[1]))
        self.HO_weights = np.random.uniform(low=-0.5, high=0.5, size=(hidden_sizes[1], output_size))

    def forward(self, inputs):
        # 计算输入到第一层隐藏层的输出
        hidden_input_1 = inputs.dot(self.IH_weights)
        hidden_output_1 = sigmoid(hidden_input_1)

        # 计算第一层隐藏层到第二层隐藏层的输出
        hidden_input_2 = hidden_output_1.dot(self.HH_weights)
        hidden_output_2 = sigmoid(hidden_input_2)

        # 计算第二层隐藏层到输出层的输出
        output_input = hidden_output_2.dot(self.HO_weights)
        outputs = sigmoid(output_input)

        return outputs

    def train(self, inputs, targets, learning_rate=0.1):
        # 前向传播计算预测值
        hidden_input_1 = inputs.dot(self.IH_weights)
        hidden_output_1 = sigmoid(hidden_input_1)
        hidden_input_2 = hidden_output_1.dot(self.HH_weights)
        hidden_output_2 = sigmoid(hidden_input_2)
        output_input = hidden_output_2.dot(self.HO_weights)
        outputs = sigmoid(output_input)

        # 计算输出层的误差和梯度
        output_error = targets - outputs
        output_gradient = output_error * (outputs * (1 - outputs))

        # 计算第二层隐藏层的误差和梯度
        hidden_error_2 = output_gradient.dot(self.HO_weights.T)
        hidden_gradient_2 = hidden_error_2 * (hidden_output_2 * (1 - hidden_output_2))

        # 计算第一层隐藏层的误差和梯度
        hidden_error_1 = hidden_gradient_2.dot(self.HH_weights.T)
        hidden_gradient_1 = hidden_error_1 * (hidden_output_1 * (1 - hidden_output_1))

        # 更新权重矩阵
        self.HO_weights += learning_rate * np.outer(hidden_output_2, output_gradient)
        self.HH_weights += learning_rate * np.outer(hidden_output_1, hidden_gradient_2)
        self.IH_weights += learning_rate * np.outer(inputs, hidden_gradient_1)

if __name__ == '__main__':
    # 设置输入、输出、隐藏层神经元数量，以及训练数据和目标值
    input_size = 4
    output_size = 2
    hidden_sizes = [8, 6]
    X_train = np.array([[0.5, 0.3, 0.7, 0.1], [0.2, 0.5, 0.1, 0.9], [0.3, 0.8, 0.6, 0.4]])
    y_train = np.array([[0, 1], [1, 0], [1, 1]])

    # 创建BP神经网络
    net = FourLayerBP(input_size, hidden_sizes, output_size)

    # 训练神经网络
    for i in range(10000):
        for j in range(X_train.shape[0]):
            net.train(X_train[j:j+1], y_train[j:j+1])

    # 测试神经网络
    X_test = np.array([[0.6, 0.2, 0.8, 0.3], [0.2, 0.7, 0.4, 1]])
    y_test = np.array([[1, 1], [1, 0]])
    for k in range(X_test.shape[0]):
        y_pred = net.forward(X_test[k:k+1])
        print('预测值:', y_pred)
        print('目标值:', y_test[k:k+1])
        plt.plot(y_pred, label='True value')
        plt.plot(y_test, label='Predicted value')
        plt.legend()
        plt.show()
        #%%
        from sklearn.neural_network import MLPRegressor
        import numpy as np
        import matplotlib.pyplot as plt
                                                                         
        X = np.array([[0., 0., 0., 0.], [1., 1., 1., 1.]])
        y = np.array([0, 1])

        clf = MLPRegressor(hidden_layer_sizes=(2,), max_iter=1000)
        clf.fit(X, y)

        print(clf.predict(X))

        plt.plot(y, label='True value')
        plt.plot(clf.predict(X), label='Predicted value')
        plt.legend()
        plt.show()