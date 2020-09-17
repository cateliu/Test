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
## 高斯光参数定义区
w_0 = 0.3125e-3
lambda1 = 1064e-9
pi = math.pi
k = 2*pi/lambda1
z_r = pi*w_0**2/lambda1
## QPD尺寸参数
r = np.array([5e-4,3e-5,5e-6])#np.arange(10e-6, 1e-4, 1e-6)# 0.5e-3
# 角度变化量
theta = np.arange(-10e-3,10e-3,1e-4)#np.array([1e-3,5e-3,10e-3])
# CCD像元尺寸

# 函数定义区
w =     lambda z:       w_0*(1+(z/z_r)**2)**0.5
R_z =   lambda z:       z/2/(z**2 + z_r**2)
phi_z = lambda z:       np.arctan(z/z_r)
A =     lambda x,y,z:   1/w(z)*np.exp(-(x**2+y**2)/w(z)**2)
Phi =   lambda x,y,z:   k*(-(x**2+y**2)*R_z(z)+z)-phi_z(z)
E_1 =   lambda x,y:         A(x,y,0)**2

def HeterInQPD(r,theta):      # 外差效率-(偏转角度,半径)
    O = np.zeros([r.size, theta.size], complex)
    for i in range(0,np.size(theta)):
        E_2 = lambda x,y:   A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))**2
        E_3_real = lambda x,y:   A(x,y,0)*A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))*np.cos(Phi(x,y,0)-Phi(x*np.cos(theta[i]),y,x*np.sin(theta[i])))
        E_3_img  = lambda x,y:   A(x,y,0)*A(x*np.cos(theta[i]),y,x*np.sin(theta[i]))*np.sin(Phi(x,y,0)-Phi(x*np.cos(theta[i]),y,x*np.sin(theta[i])))
        for m in range(0,np.size(r)):
            I1, err1 = scipy.integrate.dblquad(E_1, 0, r[m], 0, lambda x: (r[m]**2-x**2)**0.5)
            I2, err2 = scipy.integrate.dblquad(E_2, 0, r[m], 0, lambda x: (r[m]**2-x**2)**0.5)
            I3_real, err3_real = scipy.integrate.dblquad(E_3_real, 0, r[m], 0, lambda x: (r[m]**2-x**2)**0.5)
            I3_img, err3_img =  scipy.integrate.dblquad(E_3_img, 0, r[m], 0, lambda x: (r[m]**2-x**2)**0.5)
            
            I3 = I3_real - 1j*I3_img
            O[m,i] = I3/(I2*I1)**0.5
        print("完成度：",(i+1)/np.size(theta))
        Eta = np.abs(O)**2
    return Eta

Eta = HeterInQPD(r,theta)  
#%%
plt.figure(figsize = (8,6), dpi = 1200)
# plt.plot(theta*1e3,Eta[0,:], linewidth = 3)
# plt.plot(theta*1e3,Eta[0,:],theta*1e3,Eta[1,:],theta*1e3,Eta[2,:], linewidth = 3)
plt.plot(r*1e6,Eta[:,0],r*1e6,Eta[:,1],r*1e6,Eta[:,2],linewidth = 3)
plt.grid()
plt.xticks([10,30,100],size = 15)
plt.yticks(size = 15)
plt.title("外差效率随偏转角度的变化", size = 20)
plt.xlabel("偏转角度/mrad",  size = 20)
plt.ylabel("外差效率",  size = 20)
tk = plt.gca()
tk.spines["top"].set_linewidth(0)
tk.spines["right"].set_linewidth(0)
tk.spines["left"].set_linewidth(3)
tk.spines["bottom"].set_linewidth(3)
# tk.set_xscale('log')
#plt.savefig("干涉效率随偏转角度的变化")
# plt.legend(["QPD:r=0.5mm","QPD:r=0.05mm","QPD:r=0.005mm"],loc = 'right')
plt.legend(["偏转角度：1mrad","偏转角度：5mrad","偏转角度：10mrad"],loc = 'upper right')
#plt.savefig("不同QPD的像元尺寸的外差效率随偏转角度的变化")