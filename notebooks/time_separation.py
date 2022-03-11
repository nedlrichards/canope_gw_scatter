import numpy as np
import scipy.signal as sig
import scipy.io as load_mat
from math import pi

import matplotlib.pyplot as plt

from src import xponder

plt.ion()

xp = xponder()

hr = 5
min_overlap_front = 0.2
min_overlap_back = 0.02

def pick_sb(p_series):
    dB = 20 * np.log10(np.abs(p_series))
    peak_i, _ = sig.find_peaks(dB, height=125, distance=2800)
    s_i = np.sort(np.argsort(dB[peak_i])[-2:])
    peak_i = peak_i[s_i]

    return peak_i[:2]

mof_i = int(min_overlap_front * xp.fs)
mob_i = int(min_overlap_back * xp.fs)
ping_axis = np.arange(-mob_i, mof_i) / xp.fs

pings = []
dir_arrs = []

#for day in [268, 269, 270, 271]:
for day in range(250, 280):

    dir_arrs = []

    for hr in range(24):

        load_file = 'nav_' + f'{day}' + f'{hr:02}' + '5458.nc'
        try:
            p_raw, p_raw_ft = xp.load_raw(load_file)
        except:
            continue

        filt_ts = np.empty((3, xp.t_a_filt.size), dtype=np.complex128)
        sb_arr_i = []
        dir_arr = []

        for i, ts in enumerate(filt_ts):

            filt_ts[i] = xp.filter_raw(i, p_raw_ft)
            arrs = pick_sb(ts)
            sb_arr_i.append(arrs[1])
            dir_arr.append(np.abs(p_raw_ft[arrs[0]]))

        dir_arrs.append(np.array(dir_arr))

        sb_arr_i = np.array(sb_arr_i)
        arrival_order = np.argsort(sb_arr_i)

        # choose non-interference frequencies
        sb_arrs = []
        if (sb_arr_i[arrival_order[2]] - sb_arr_i[arrival_order[1]]) > mob_i:
            sb_arrs.append(arrival_order[2])

        if ((sb_arr_i[arrival_order[1]] - sb_arr_i[arrival_order[0]]) > mob_i) \
            & ((sb_arr_i[arrival_order[2]] - sb_arr_i[arrival_order[1]]) > mof_i):
            sb_arrs.append(arrival_order[1])

        if (sb_arr_i[arrival_order[1]] - sb_arr_i[arrival_order[0]]) > mof_i:
            sb_arrs.append(arrival_order[0])

        # create ping library
        e_ping = np.full((3, ping_axis.size), np.nan, dtype=np.complex128)

        for sb in sb_arrs:
            e_ping[sb] = filt_ts[sb, sb_arr_i[sb] - mob_i: sb_arr_i[sb] + mof_i]

        pings.append(e_ping)

dir_arrs = np.array(dir_arrs)
dir_amp = 10 * np.log10(np.mean(dir_arrs ** 2, axis=0))

pings = np.array(pings)

# total mean
mean_int = np.nanmean(np.abs(pings) ** 2, axis=0)
z_i = np.argmin(np.abs(ping_axis))

fig, ax = plt.subplots()
#ax.plot(ping_axis, 10 * np.log10(mean_int).T - dir_amp.T)
ax.plot(ping_axis, 10 * (np.log10(mean_int[0, :]) - np.log10(mean_int[0, z_i])))
ax.plot(ping_axis, 10 * (np.log10(mean_int[1, :]) - np.log10(mean_int[1, z_i])))
ax.plot(ping_axis, 10 * (np.log10(mean_int[2, :]) - np.log10(mean_int[2, z_i])))

# daily mean
num_days = int(pings.shape[0] / (4 * 24))
day_pings = np.array_split(pings, num_days, axis=0)

plot_index = 0

fig, ax = plt.subplots()
for day in day_pings:
    num_pings = np.sum(~np.isnan(day[:, plot_index, 0]), axis=0)
    print(num_pings)
    if num_pings < 24:
        continue
    day_mean_int = np.nanmean(np.abs(day[:, plot_index, :]) ** 2, axis=0)

    norm = 20 * np.log10(mean_int[plot_index, z_i])

    ax.plot(ping_axis, 10 * np.log10(day_mean_int) - 10 * np.log10(mean_int[plot_index, z_i]))

ax.set_ylim(-40, 3)
ax.grid()

plot_index = 1

fig, ax = plt.subplots()
for day in day_pings[:-2]:
    num_pings = np.sum(~np.isnan(day[:, plot_index, 0]), axis=0)
    print(num_pings)
    if num_pings < 24:
        continue
    day_mean_int = np.nanmean(np.abs(day[:, plot_index, :]) ** 2, axis=0)

    norm = 20 * np.log10(mean_int[plot_index, z_i])

    ax.plot(ping_axis, 10 * np.log10(day_mean_int) - 10 * np.log10(mean_int[plot_index, z_i]))

ax.set_ylim(-40, 3)
ax.grid()



1/0
plot_index = 2

fig, ax = plt.subplots()
for day in day_pings[3:]:
    num_pings = np.sum(~np.isnan(day[:, plot_index, 0]), axis=0)
    print(num_pings)
    if num_pings < 24:
        continue
    day_mean_int = np.nanmean(np.abs(day[:, plot_index, :]) ** 2, axis=0)

    norm = 20 * np.log10(mean_int[plot_index, z_i])

    ax.plot(ping_axis, 10 * np.log10(day_mean_int) - 10 * np.log10(mean_int[plot_index, z_i]))

ax.set_ylim(-40, 3)
ax.grid()

