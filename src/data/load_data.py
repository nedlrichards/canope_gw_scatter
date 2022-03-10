import numpy as np
import scipy.signal as sig
import scipy.io as load_mat
import netCDF4
from math import pi

from os import path

data_path = 'data/raw'

class xponder:

    def __init__(self):
        """setup common parameters"""
        # file samples
        self.fs = 1e7 / 256
        self.t_load = (5, 12)
        samples = np.array(self.t_load) * self.fs
        self.samples = samples.astype(np.int_)
        self.t_a = np.arange(self.samples[0], self.samples[1]) / self.fs

        # frequency domain specifications
        self.nfft = 2 ** (int(np.log2(self.t_a.size)) + 1)
        self.f_a = np.arange(self.nfft) * self.fs / self.nfft

        # hydrophone sensitivities
        self.bits2volts = 2.5 / 2 ** 23
        self.ampgain = 10 ** (12 / 20)
        self.hysens = 10 ** (-168 / 20)

        # bandpass filter specifications
        self.bp_numtaps = 2 ** 9
        self.bp_bw = 0.5e3
        self.bp_trans_width = 0.2e3
        self.ping_fc = [11e3, 11.5e3, 12e3]

        # replica pulse specifications
        dt = 1 / self.fs
        self.pulse_T = 0.009
        pulse_N = self.pulse_T // dt
        self.pulse_t_a = dt * np.arange(pulse_N)

        # filter and pulse replica banks
        filter_bank = []
        pulse_bank = []

        for f in self.ping_fc:
            bp_edges = [0, f - self.bp_bw / 2 - self.bp_trans_width,
                        f - self.bp_bw / 2, f + self.bp_bw / 2,
                        f + self.bp_bw / 2 + self.bp_trans_width, self.fs / 2]
            filter_bank.append(sig.remez(self.bp_numtaps, bp_edges,
                                         [0, 1, 0], Hz=self.fs))
            pulse_bank.append(np.sin(2 * pi * f * self.pulse_t_a))

        self.filter_bank = np.array(filter_bank)
        self.pulse_bank = np.array(pulse_bank)

        # Fourier transform used in filtering
        self.filter_bank_ft = np.fft.fft(self.filter_bank, n=self.nfft)
        self.pulse_bank_ft = np.fft.fft(self.pulse_bank, n=self.nfft)

        # clip filter result
        self.num_edge = int(np.maximum(self.bp_numtaps, pulse_N) - 1)
        self.t_a_filt = self.t_a[self.num_edge: ]

        # surface bounce bounds
        self.sb_tbounds = (-0.1, 0.5)
        num_sb = np.ceil((self.sb_tbounds[1] - self.sb_tbounds[0]) * self.fs)
        self.sb_t_a = np.arange(num_sb) / self.fs + self.sb_tbounds[0]

    def load_raw(self, load_file):
        """Load ping data with basic processing"""
        nav_path = path.join(data_path, 'nav_ping')
        nc_file = netCDF4.Dataset(path.join(nav_path, load_file), 'r')

        p_raw = nc_file.variables['samples'][self.samples[0]: self.samples[1]]
        p_raw = np.asarray(p_raw, dtype=np.float64)

        p_raw *= self.bits2volts    #  Convert bits to volts
        p_raw /= (self.ampgain * self.hysens)   #  Convert volts to input pressure (uPa)
        p_raw_ft = np.fft.fft(p_raw, n=self.nfft)
        return p_raw, p_raw_ft

    def filter_raw(self, fc_i, p_raw_ft):
        """match filter, bandpass and demodulate"""

        filt_ft = self.filter_bank_ft[fc_i]
        pulse_ft = self.pulse_bank_ft[fc_i]

        p_filt = p_raw_ft * filt_ft * pulse_ft.conj()

        # hilbert transform
        p_filt[self.f_a >= self.fs / 2] = 0 + 0j
        p_filt = np.fft.ifft(p_filt)[:self.t_a.size]
        # demodulate by carrier
        p_filt = p_filt * np.exp(-2j * pi * self.ping_fc[fc_i] * self.t_a)

        p_filt = p_filt[self.num_edge: ]
        return p_filt

    def pick_first_arrival(self, p_filt):
        """just pick second arrival for now"""
        abv = np.abs(p_filt)
        i_sep = int(0.1 * self.fs)
        peaks = sig.find_peaks(abv, distance=i_sep, height=5e6)[0]
        dir_peak = self.t_a_filt[peaks[0]]
        return dir_peak


    def window_sb(self, p_filt):
        """just pick second arrival for now"""
        abv = np.abs(p_filt)
        i_sep = int(0.1 * self.fs)
        peaks = sig.find_peaks(abv, distance=i_sep, height=5e6)[0]

        sb_peak = peaks[1]

        temp_t_a = self.t_a_filt - self.t_a_filt[sb_peak]
        time_inds = (temp_t_a >= self.sb_tbounds[0]) \
                  & (temp_t_a <= self.sb_tbounds[1])

        return p_filt[time_inds]

