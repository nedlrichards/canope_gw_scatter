import numpy as np
import scipy.signal as sig
import scipy.io as load_mat
from math import pi

import matplotlib.pyplot as plt

from src import xponder

plt.ion()

xp = xponder()

#for day in range(250, 260):
for day in [253]:
    arr_11 = []
    arr_115 = []
    arr_12 = []


    for hr in np.arange(23):

        load_file = 'nav_' + f'{day}' + f'{hr:02}' + '5458.nc'
        try:
            p_raw, p_raw_ft = xp.load_raw(load_file)
        except FileNotFoundError:
            continue

        p_filt_11 = xp.filter_raw(0, p_raw_ft)
        arr_11.append(xp.window_sb(p_filt_11))

        p_filt_115 = xp.filter_raw(1, p_raw_ft)
        arr_115.append(xp.window_sb(p_filt_115))

        p_filt_12 = xp.filter_raw(2, p_raw_ft)
        arr_12.append(xp.window_sb(p_filt_12))

    arr_11 = np.array(arr_11)
    arr_115 = np.array(arr_115)
    arr_12 = np.array(arr_12)

    t_0_i  = np.argmin(np.abs(xp.sb_t_a))
    norm_11 = np.mean(np.abs(arr_11[:, t_0_i]))
    norm_115 = np.mean(np.abs(arr_115[:, t_0_i]))
    norm_12 = np.mean(np.abs(arr_12[:, t_0_i]))

    fig, ax = plt.subplots()
    ax.plot(xp.sb_t_a, 20 * np.log10(np.mean(np.abs(arr_11), axis=0)).T - 20 * np.log10(norm_11))
    ax.plot(xp.sb_t_a, 20 * np.log10(np.mean(np.abs(arr_115), axis=0)).T - 20 * np.log10(norm_115))
    ax.plot(xp.sb_t_a, 20 * np.log10(np.mean(np.abs(arr_12), axis=0)).T - 20 * np.log10(norm_12))
    ax.grid()


1/0
fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6.5, 6))
ax[0].plot(xp.sb_t_a, 20 * np.log10(np.abs(arr_11)).T - 20 * np.log10(norm_11), 'C0')
ax[1].plot(xp.sb_t_a, 20 * np.log10(np.abs(arr_115)).T - 20 * np.log10(norm_115), 'C1')
ax[2].plot(xp.sb_t_a, 20 * np.log10(np.abs(arr_12)).T - 20 * np.log10(norm_12), 'C2')

ax[0].set_ylim(-30, 5)


