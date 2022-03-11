import numpy as np
import scipy.signal as sig
import scipy.io as load_mat
from math import pi

import matplotlib.pyplot as plt

from src import xponder

#plt.ion()

xp = xponder()

hr = 5
min_overlap_front = 0.2
min_overlap_back = 0.02

def pick_sb(p_series):
    dB = 20 * np.log10(np.abs(p_series))
    peak_i, _ = sig.find_peaks(dB, height=125, distance=2800)
    s_i = np.sort(np.argsort(dB[peak_i])[-2:])
    peak_i = peak_i[s_i]

    return peak_i[1]

for day in [269]:
    for hr in range(24):

        load_file = 'nav_' + f'{day}' + f'{hr:02}' + '5458.nc'
        p_raw, p_raw_ft = xp.load_raw(load_file)

        p_filt_11 = xp.filter_raw(0, p_raw_ft)
        peaks_11 = pick_sb(p_filt_11)

        p_filt_115 = xp.filter_raw(1, p_raw_ft)
        peaks_115 = pick_sb(p_filt_115)

        p_filt_12 = xp.filter_raw(2, p_raw_ft)
        peaks_12 = pick_sb(p_filt_12)

        # choose non-interference frequencies
        sb_arr_i = np.array([peaks_11, peaks_115, peaks_12])
        arrival_order = np.argsort(sb_arr_i)

        mof_i = int(min_overlap_front * xp.fs)
        mob_i = int(min_overlap_back * xp.fs)

        sb_arrs = []
        if (sb_arr_i[arrival_order[2]] - sb_arr_i[arrival_order[1]]) > mob_i:
            sb_arrs.append(arrival_order[2])

        if ((sb_arr_i[arrival_order[1]] - sb_arr_i[arrival_order[0]]) > mob_i) \
            & ((sb_arr_i[arrival_order[2]] - sb_arr_i[arrival_order[1]]) > mof_i):
            sb_arrs.append(arrival_order[1])

        if (sb_arr_i[arrival_order[1]] - sb_arr_i[arrival_order[0]]) > mof_i:
            sb_arrs.append(arrival_order[0])

        fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6.5, 6))
        ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)), 'C0')
        ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)) - 24, '0.4')
        ax[0].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)) - 24, '0.4')

        ax[0].plot(xp.t_a_filt[peaks_11], 20 * np.log10(np.abs(p_filt_11))[peaks_11], 'o', color='C0')
        if 0 in sb_arrs:
            ax[0].plot(xp.t_a_filt[peaks_11], 20 * np.log10(np.abs(p_filt_11))[peaks_11], '*', color='C0', markersize=12)

        ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)) - 24, '0.4')
        ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)), 'C1')

        ax[1].plot(xp.t_a_filt[peaks_115], 20 * np.log10(np.abs(p_filt_115))[peaks_115], 'o', color='C1')
        if 1 in sb_arrs:
            ax[1].plot(xp.t_a_filt[peaks_115], 20 * np.log10(np.abs(p_filt_115))[peaks_115], '*', color='C1', markersize=12)

        ax[1].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)) - 24, '0.4')
        ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_11)) - 24, '0.4')
        ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_115)) - 24, '0.4')
        ax[2].plot(xp.t_a_filt, 20 * np.log10(np.abs(p_filt_12)), 'C2')

        ax[2].plot(xp.t_a_filt[peaks_12], 20 * np.log10(np.abs(p_filt_12))[peaks_12], 'o', color='C2')
        if 2 in sb_arrs:
            ax[2].plot(xp.t_a_filt[peaks_12], 20 * np.log10(np.abs(p_filt_12))[peaks_12], '*', color='C2', markersize=12)
        ax[0].set_ylim(110, 160)
        ax[0].set_xlim(7.5, 9.0)

        fig.savefig('notebooks/figures/' + load_file.split('.')[0] + '.png', dpi=300)
        plt.close(fig)


