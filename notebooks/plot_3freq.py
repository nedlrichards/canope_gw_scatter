import numpy as np
import scipy.signal as sig
import scipy.io as load_mat
from math import pi

import matplotlib.pyplot as plt

from src import xponder

#plt.ion()

xp = xponder()

for hr in range(24):
    load_file = 'nav_253' + f'{hr:02}' + '5458.nc'
    try:
        p_raw, p_raw_ft = xp.load_raw(load_file)
    except:
        continue

    p_filt_11 = xp.filter_raw(0, p_raw_ft)
    p_win_11 = xp.window_sb(p_filt_11)
    p_filt_115 = xp.filter_raw(1, p_raw_ft)
    p_win_115 = xp.window_sb(p_filt_115)
    p_filt_12 = xp.filter_raw(2, p_raw_ft)
    p_win_12 = xp.window_sb(p_filt_12)

    fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6.5, 6))
    ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)).T, 'C0')
    ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)).T - 24, '0.4')
    ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)).T - 24, '0.4')
    ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)).T - 24, '0.4')
    ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)).T, 'C1')
    ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)).T - 24, '0.4')
    ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)).T - 24, '0.4')
    ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)).T - 24, '0.4')
    ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)).T, 'C2')
    ax[0].set_ylim(110, 160)
    ax[0].set_xlim(7.5, 9.0)

    fig.savefig('notebooks/figures/' + load_file.split('.')[0] + '.png', dpi=300)
    plt.close(fig)



