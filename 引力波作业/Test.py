# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz

# the ipython magic below must be commented out in the .py file, since it doesn't work.
# %matplotlib inline
# %config InlineBackend.figure_format = 'retina'
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import h5py

# LIGO-specific readligo.py 
import readligo as rl
#----------------------------------------------------------------
NRtime, NR_H1 = np.genfromtxt('observed-H.txt').transpose()
# NRtime = NRtime - 6.9e-3
# NR_H1 = -NR_H1
fs = 4096
# pick a shorter FTT time interval, like 1/8 of a second:
NFFT = fs/8
# and with a lot of overlap, to resolve short-time features:
NOVL = NFFT*15/16
plt.figure(2)
plt.plot(NRtime,NR_H1)
# spec_H1, freqs, bins, im = plt.specgram(NR_H1, NFFT=NFFT, Fs=fs,xextent=[0.25,0.47])
plt.xlabel('time (s) since ')
plt.ylabel('Frequency (Hz)')
# plt.colorbar()
# plt.axis([0.25, 0.5, 32, 512])
plt.title('aLIGO H1 strain data near GW150914')
plt.show()
plt.savefig('GW150914_H1_spectrogram.png')
