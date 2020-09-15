# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 19:39:13 2020

计算高斯干涉光在QPD的平均相位随偏转角的变化
@author: cate_liu
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate

## 高斯光参数定义区
w_0 = 0.3125e-3
lambda1 = 1064e-9
pi = math.pi
k = 2*pi/lambda1
z_r = pi*w_0**2/lambda1
## QPD尺寸参数
r = 0.5e-3
# 函数定义区

w = lambda z:w_0*(1+(z/z_r)**2)**0.5
R_z = lambda z:z/2/(z**2 + z_r**2)
phi_z = lambda z:np.arctan(z/z_r)
A = lambda x,y,z:1/w(z)*np.exp(-(x**2+y**2)/w(z)**2)
Phi = lambda x,y,z:k*(-(x**2+y**2)*R_z(z)+z)-phi_z(z)

# 角度变化量
theta = np.arange(-10e-3,10e-3,1e-4)

# 变换矩阵
Trans = lambda alpha:np.array([np.cos(alpha),0,np.sin(alpha)],[0,1,0],[-np.sin(alpha),0,np.cos(alpha)])

for i in range(0,np.size(theta)):
    T1 = Trans(theta[i])