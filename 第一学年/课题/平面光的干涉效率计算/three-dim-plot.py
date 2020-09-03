from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.io import scio

m = scio.loadmat('Plane_Beams_Data.mat')
m.keys()