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
#%%
# 字体设置
fontSongti = {'family': 'serif',
              'size':20,
    }
#%%
## 高斯光参数定义区
w_0 = 0.3125e-3
lambda1 = 1064e-9
pi = math.pi
k = 2*pi/lambda1
z_r = pi*w_0**2/lambda1
## QPD尺寸参数
r = 0.5e-3
# CCD像元尺寸
l1 = 5e-6
l2 = 50e-6
# 函数定义区

w =     lambda z:       w_0*(1+(z/z_r)**2)**0.5
R_z =   lambda z:       z/2/(z**2 + z_r**2)
phi_z = lambda z:       np.arctan(z/z_r)
A =     lambda x,y,z:   1/w(z)*np.exp(-(x**2+y**2)/w(z)**2)
Phi =   lambda x,y,z:   k*(-(x**2+y**2)*R_z(z)+z)-phi_z(z)

# 角度变化量
theta = np.arange(-10e-3,10e-3,1e-4)

# 变换矩阵
# Trans = lambda alpha:       np.array([np.cos(alpha),0,np.sin(alpha)],[0,1,0],[-np.sin(alpha),0,np.cos(alpha)])
E_1 =   lambda x,y:         A(x,y,0)**2
I1, err1 = scipy.integrate.dblquad(E_1, 0, r, 0, lambda x: (r**2-x**2)**0.5)
I11, err11 = scipy.integrate.dblquad(E_1, -l1, l1, lambda x: -l1, lambda x:l1)
I12, err12 = scipy.integrate.dblquad(E_1, -l2, l2, lambda x: -l2, lambda x:l2)
O = np.zeros(np.size(theta), complex)
O1 = O.copy()
O2 = O.copy()
for i in range(0,np.size(theta)):
    E_2 = lambda x,y:   A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))**2
    E_3_real = lambda x,y:   A(x,y,0)*A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))*np.cos(Phi(x,y,0)-Phi(x*np.cos(theta[i]),y,x*np.sin(theta[i])))
    E_3_img = lambda x,y:   A(x,y,0)*A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))*np.sin(Phi(x,y,0)-Phi(x*np.cos(theta[i]),y,x*np.sin(theta[i])))
    # QPD积分
    I2, err2 = scipy.integrate.dblquad(E_2, 0, r, 0, lambda x: (r**2-x**2)**0.5)
    I3_real, err3_real = scipy.integrate.dblquad(E_3_real, 0, r, 0, lambda x: (r**2-x**2)**0.5)
    I3_img, err3_img =  scipy.integrate.dblquad(E_3_img, 0, r, 0, lambda x: (r**2-x**2)**0.5)
    # CCD像元积分
    I21, err21 = scipy.integrate.dblquad(E_2, -l1, l1, lambda x: -l1, lambda x:l1)
    I22, err22 = scipy.integrate.dblquad(E_2, -l2, l2, lambda x: -l2, lambda x:l2)
    I31_real, err31_real    = scipy.integrate.dblquad(E_3_real, -l1, l1, lambda x: -l1, lambda x:l1)
    I31_img, err31_img      = scipy.integrate.dblquad(E_3_img,  -l1, l1, lambda x: -l1, lambda x:l1)
    I32_real, err32_real    = scipy.integrate.dblquad(E_3_real, -l2, l2, lambda x: -l2, lambda x:l2)
    I32_img, err32_img      = scipy.integrate.dblquad(E_3_img,  -l2, l2, lambda x: -l2, lambda x:l2)
    
    I31 = I31_real-1j*I31_img
    I32 = I32_real-1j*I32_img
    I3 = I3_real-1j*I3_img
    
    O[i] = I3/(I2*I1)**0.5
    O1[i] = I31/(I21*I11)**0.5
    O2[i] = I32/(I22*I12)**0.5
#%%
plt.figure(figsize = (8,6), dpi = 1200)
plt.plot(theta*1e3,O.real**2, linewidth = 3)
plt.grid()
plt.xticks([-10,-2,0,2,10],size = 15)
plt.yticks(size = 15)
plt.title("外差效率随偏转角度的变化", size = 20)
plt.xlabel("偏转角度/mrad",  size = 20)
plt.ylabel("外差效率",  size = 20)
tk = plt.gca()
tk.spines["top"].set_linewidth(0)
tk.spines["right"].set_linewidth(0)
tk.spines["left"].set_linewidth(3)
tk.spines["bottom"].set_linewidth(3)
plt.savefig("干涉效率随偏转角度的变化")
#plt.legend(["QDP:r=0.5mm","矩形:l=10um","矩形:l=100um"],loc = 'right')
#plt.savefig("不同探测面的外差效率随偏转角度的变化")