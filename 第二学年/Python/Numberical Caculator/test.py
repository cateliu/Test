import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

w_0 = 0.3125e-3
lambda1 = 1064e-9
pi = math.pi
k = 2*pi/lambda1
z_r = pi*w_0**2/lambda1

z = np.arange(-1,1,1e-2)
theta = (lambda1/pi/z_r)**0.5
z1 = theta*z

def w(z):
    return w_0*(1+(z/z_r)**2)**0.5

def R_z(z):
    return z/2/(z**2 + z_r**2)

def phi_z(z):
    return np.arctan(z/z_r)

def A(x,y,z):
    return 1/w(z)*np.exp(-(x**2+y**2)/w(z)**2)
def Phi(x,y,z):
    return k*(-(x**2+y**2)*R_z(z)+z)-phi_z(z)

x = np.linspace(-1e-2,1e-2,640)
wx = A(x,0*x,0*x)*np.exp(-1j*Phi(x,0*x,0*x))

# 验证高斯光模型的正确性
fig = plt.figure()
f1 = fig.add_subplot(121)
f2 = fig.add_subplot(122)
f1.plot(z,w(z),'k',z,-w(z),'k',z,theta*z,'k--',z,-theta*z,'k--')
f1.set_title("光束半径在z轴的变化")
f1.set_xticks([0])
f1.set_yticks([0])
f2.plot(x, wx.real, 'k')
f2.set_title("光强在x轴的变化")
f2.set_xticks([0])
f2.set_yticks([0])


# 求解相位的平面分布

y = np.linspace(-0.8e-2,0.8e-2,512)
Y, X = np.meshgrid(x,y)
phi_l = Phi(X,Y,5e-2*np.ones(X.shape)) # z = 5cm 处的相位分布

fig2 = plt.figure()
ax2 = fig2.add_subplot(121, projection = '3d')
ax2.plot_surface(X,Y,phi_l-2.949e5,  cmap = 'rainbow')
ax2.set_xticks([-0.008,0,0.008])
ax2.set_yticks([-0.008,0,0.008])
ax2.set_title("光束未偏转的相位分布")

# 光束旋转后的相位分布
alpha = 10e-3   # 绕x轴旋转alpha角度
beta = 10e-3    # 绕y轴旋转beta角度

E_L = A(X,Y,5e-2*np.ones(X.shape))*np.exp(-1j*Phi(X,Y,5e-2*np.ones(X.shape))) # 本地光的高斯模型
Ralpha = np.array([(math.cos(alpha), 0, math.sin(alpha))
    ,(0, 1, 0)
    ,(-math.sin(alpha), 0, math.cos(alpha))])
Rbeta = np.array([(1,0,0)
    ,(0,math.cos(beta),math.sin(beta))
    ,(0,-math.sin(beta),math.cos(beta))])
X1 = np.zeros(X.shape)
Y1 = np.zeros(X.shape)
Z1 = np.zeros(X.shape)      # 旋转后的坐标与旋转前的坐标之间的关系
for m in range(X.shape[0]):
    for j in range(X.shape[1]):
        p = np.array([X[m,j], Y[m,j], 5e-2])@Ralpha@Rbeta
        X1[m,j], Y1[m,j],Z1[m,j] = p

phi_r = Phi(X1, Y1, Z1)     # 旋转后的相位分布
E_Phi = phi_l - phi_r       # 干涉信号的相位分布
E_R= A(X1, Y1, Z1)*np.exp(-1j*phi_r)

E_S = (E_L + E_R)*np.conj(E_L + E_R)

ax3 = fig2.add_subplot(122, projection = '3d')
ax3.plot_surface(X,Y,E_Phi, cmap = 'rainbow')
ax3.set_xticks([-0.8e-2,0, 0.8e-2])
ax3.set_yticks([-1e-2,0, 1e-2])

fig3 = plt.figure()
ax4 = Axes3D(fig3)
ax4.plot_surface(X,Y, E_S.real, rstride = 1, cmap = 'rainbow')

plt.show()