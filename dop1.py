import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#plt.rcParams["figure.figsize"] = (12, 12) # (w, h)
F = 10e9
c0 = 3e8

lamb = c0/F
PRT = 1e-3
hits_scan = 256

delta_v=lamb/(2*hits_scan*PRT)

v_amb = lamb/(4 * PRT)

vel_axis = np.linspace(-v_amb,v_amb,hits_scan)

dR = 1

R_max = 15e3

Range_bins =  int(np.round(R_max/dR))

range_axis = np.linspace(0, Range_bins - 1, Range_bins)*dR

target_pos = 5e3

tau = 1e3
vt = 3

data = []
target_range_index = []

target_pos = 5e3
target_width = 1e3

target_vel = 3

for i in range(hits_scan):
    target_pos = 5e3 + vt * PRT * (i - 1)
    target_range = target_pos + np.arange(0, tau-1, dR)
    ph = 2*np.pi*2*target_range/lamb
    target_range_index = np.round(target_range/dR)
    target_range_index = target_range_index.astype(int)
    signal = np.zeros(Range_bins)
    signal = signal.astype(complex)
    signal[target_range_index] = np.exp(1j * ph)
    data.append(signal)

data = np.reshape(data, (hits_scan, Range_bins))
plt.imshow(np.abs(data), extent=(1, Range_bins, 1, hits_scan))
plt.show();