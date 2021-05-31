import numpy as np
from scipy.fftpack import fft, ifft, fftshift

def DS_simulatorV2(SNR, m0, mu, sigma, n, v_amb, N):
    '''This is the code for the time domain (slow time) form of a Gaussian shaped velocity spectrum for weather radars.
    This function produces a high definition (HD) frequency domain signal (with n points) internally and based on the
    number of sweeps available physically (N) for a weather radar, it down-samples the HD spectrum using 
    the area under the curve principles so that the power of the spectrum is conserved. It returns the time domain 
    variant of the N point frequency domain Gaussian spectrum. It is usually referred to as the slow time domain
    signal for Doppler processing.


    Input: 
    SNR in linear scale (SNR), 
    zeroth moment of the Gaussian spectrum (m0), 
    Mean velocity (mu),
    Doppler Spectrum Width (sigma), 
    Number of points on velocity axis (n, usually high resolution like 512 or 1024),
    Maximum unambiguous velocity (v_amb), 
    Number of points needed in the output of the velocity spectrum (N). 

    Output: 
    Time domain data with N points (N should be less than n) (data), 
    N point frequency domain Gaussian velocity amplitude spectrum (data_f), 
    n point HD amplitude Gaussian spectrum (data_f_full), 
    time domian variant of the HD Gaussian spectrum with n points (data_full). 

    '''

    axis = np.linspace(-n / 2, n / 2 - 1, n) / n
    vel_axis = 2 * v_amb * axis

    dv = vel_axis[1] - vel_axis[0]

    X_full = np.random.uniform(0, 1, size=(1, n))
    Theta_full = 2 * np.pi * np.random.uniform(0, 1, size=(1, n))

    S_ = m0 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(-(vel_axis - mu) ** 2 / (2 * sigma ** 2))

    Noise_full = np.sum(S_) / (n * SNR)

    P_full = -(S_ + Noise_full) * np.log(X_full)

    data_f_full = np.sqrt(P_full) * np.exp(1j * Theta_full)

    data_full = ifft(fftshift(np.sqrt(n) * data_f_full))

    axis_permitted = np.linspace(-N / 2, N / 2 - 1, N) / N
    vel_axis_permitted = 2 * v_amb * axis_permitted

    idx = np.array(np.zeros(N))

    for i in range(N):
        idx[i] = (np.abs(vel_axis - vel_axis_permitted[i])).argmin()

    idx_for_integral = np.round(np.mean([idx[0:-1], idx[1:]], axis=0))
    idx_for_integral = np.append(idx_for_integral, idx[-1])
    idx_for_integral = np.insert(idx_for_integral, 0, idx[0])

    idx_for_integral = idx_for_integral.astype(int)

    # print(np.size(idx_for_integral))
    S = np.array(np.zeros(N))

    for k in range(N):
        num = idx_for_integral[k + 1] - idx_for_integral[k] + 1
        S[k] = np.sum(S_[idx_for_integral[k]:idx_for_integral[k + 1]] * dv) / (num * dv)

    idxn = np.linspace(0, n - 1, n)
    idxN = np.linspace(np.min(idxn), np.max(idxn), N)

    print(np.size(idxn))
    print(np.size(idxN))
    print(np.size(X_full))
    X = np.array(np.zeros(N))
    X = np.interp(idxN, idxn, X_full)
    Theta = np.array(np.zeros(N))
    Theta = np.interp(idxN, idxn, Theta_full)
    Noise = np.sum(S) / (N * SNR)

    P = -(S + Noise) * np.log(X)
    data_f = np.sqrt(P) * np.exp(1j * Theta)
    data = ifft(fftshift(np.sqrt(N) * data_f))

    return (data, data_f, data_f_full, data_full, X_full, Theta_full)