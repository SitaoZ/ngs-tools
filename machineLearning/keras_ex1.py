#-*-coding:utf-8-*-

from keras.models import Sequential
from keras.layers import Dense, Activation
model = Sequential()

model.add(Dense(32, activation='relu', input_dim=100))
model.add(Dense(1, activation='sigmoid'))
model.compile(optimizer='rmsprop',
              loss='binary_crossentropy',
              metrics=['accuracy'])

# 生成虚假数据
import numpy as np
data = np.random.random((1000, 100))
labels = np.random.randint(2, size=(1000, 1))

# 训练模型，迭代数据（每个batch包含32个样本）
model.fit(data, labels, epochs=10, batch_size=32)