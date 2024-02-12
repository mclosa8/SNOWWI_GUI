import numpy as np
from os import stat
from os.path import getmtime
from glob import glob, iglob
from struct import unpack
from time import gmtime, time, sleep, strftime
from sys import argv, exit

from platform import system
if system() == 'Darwin':   # Mac OS
    from matplotlib import use
    use('MacOSX')  # Workaround for animation problem in default backend

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
from matplotlib.widgets import Button

import os


def second_to_last(data_dir):
    pat = data_dir + '/' + '*.dat'
    try:
        lof = sorted(iglob(pat), key=os.path.getctime)
        return lof[-2]
    except:
        return ''


def parse_params(params_fname):
    inpFile = open(params_fname, 'r')
    while 1:
        l = inpFile.readline()

        if len(l) == 0:
            break

        lparts = l.rpartition(' ')

        if (len(lparts) > 0):
            var = lparts[0].strip(' \r\n')
            val = lparts[2].strip(' \r\n')
            if var == 'prf':
                prf = float(val)
            if var == 'tp':
                tp = float(val)
            if var == '--N0':
                Nrx = int(val)*4
            if var == 'fsrx':
                fsrx = float(val)
            if var == 'numGates':  # 2014-5 version
                Nrx = int(val)
                tp = 0.001
            if var == 'fs':   # 2014-5 version
                fsrx = float(val)
                tp = 0.001
            if var == 'decimate_by':
                decimate_test = int(val)
                if decimate_test > 1:
                    print(
                        "!!!!!!!!!!!!!!!!!!!!!!WARNING: YOU ARE DECIMATING DATA!!!!!!!!!!!!!!!!!")

    inpFile.close()
    if Nrx == 0:
        Nrx = int(fsrx * tp) + 1

    prf = 1e3
    tp = 1/prf
    fsrx = 491.52e6
    return prf, tp, Nrx, fsrx


def prep(fsrx, tp, Nrx, rmax, ver, Navg, header):
    # Set up constants prior to reading the data.
    c = 3e8
    B = 80e6
    r0 = 0  # Range correction - half of the chirp length TODO dinamically
    N = Nrx - header  # - 1  # True frame size without the extra metadata
    n = np.arange(0, N)
    nfft = np.uint32(2**(np.ceil(np.log2(N))))   # next higher power of 2
    nfft *= 2   # double fft length => interpolates in range
    Nf = nfft/2 + 1
    Nmax = int(np.ceil(2*rmax / c * fsrx))
    Nmax = Nrx - header
    # print(type(Nmax), Nmax)
    if Nmax > Nf:
        Nmax = int(Nf)
        rmax = Nmax * c/(2.0*B) * (N*1.0) / nfft - r0
    x = np.arange(0, nfft//2+1)
    f = x * fsrx / nfft
    # print(f.shape)
    r = x * c/(2.0*B) * (N*1.0) / nfft - r0
    w = .5 * (1-np.cos(2*np.pi * np.linspace(0, 1., N, endpoint=False)))
    w = w / N / 32767  # correct amplitude and normalize
    if ver == 0:
        w[0] = 0    # Force the metadata sample to be excluded
    W = np.tile(w, (Navg, 1))
    num_rows = 100
    M1 = np.zeros((num_rows, Nmax))
    M2 = np.zeros((num_rows, Nmax))

    # print(M1.shape, M2.shape)

    return c, B, r0, N, n, nfft, Nf, Nmax, rmax, x, f, r, w, W, M1, M2


def setupFig(N, r0, rmax, M1, M2, fs):

    fig = plt.figure('SNOWWI DAQ', figsize=(26, 14), dpi=80)

    ax1 = plt.subplot(411)
    Iraw, Qraw = plt.plot([], [], 'C0--', [], [], 'C1--')
    qPlot, iPlot = plt.plot([], [], 'r-', [], [], 'k-')
    plt.xlabel('Sample', fontsize=9)
    plt.ylabel('Bit value', fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.axis([-N*.005, N*1.005, -2**11 - 20, 2**11 + 20])
    title1 = plt.title('Radar data', fontsize=9)
    plt.grid('on')

    ax2 = plt.subplot(412)
    r2Plot, = plt.plot([], [], 'r-')
    r1Plot, = plt.plot([], [], 'k-')
    plt.axis([0, fs/2, 0, 80])
    plt.xlabel('Frequency [Hz]', fontsize=9)
    plt.ylabel('Spectrum [dBm]', fontsize=9)
    # plt.xlabel('Range [m]',fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.grid('on')

    min_dB = -25   # M color values
    max_dB = +25

    ax3 = plt.subplot(425)
    im1Plot = plt.imshow(M1, origin='upper',
                         interpolation='nearest', cmap='jet', aspect='auto')
    im1Plot.set_clim((min_dB, max_dB))
    ax3.xaxis.set_visible(0)
    plt.yticks(fontsize=9)
    plt.ylabel('Channel 0', fontsize=9)

    ax4 = plt.subplot(427)
    im2Plot = plt.imshow(M2, origin='upper',
                         interpolation='nearest', cmap='jet', aspect='auto')
    im2Plot.set_clim((min_dB, max_dB))
    ax4.xaxis.set_visible(0)
    plt.yticks(fontsize=9)
    plt.ylabel('Channel 1', fontsize=9)

    # Plots for histogram
    ax5 = plt.subplot(426)
    hist0 = plt.plot([], [])
    plt.xlabel('Digital value')
    plt.ylabel('Pdf')

    ax6 = plt.subplot(428)
    hist1 = plt.plot([], [])
    plt.xlabel('Digital value')
    plt.ylabel('Pdf')

    fig.tight_layout()

    return fig, iPlot, qPlot, r1Plot, r2Plot, im1Plot, im2Plot, title1, hist0, hist1, Iraw, Qraw


def readFrame(ver, file1, N):

    bytes = (N+1)*2*2+16  # 1 channel, real+imag, 16bits
    d = file1.read(bytes)
    if len(d) < bytes:
        exit('[python] Incomplete frame; reached end of data')
    fmt = ('%sh' % (N*2))   # e.g. '500h'  for 500 int16
    d2 = d[16:16+N*4]     # actually this should depend on ver
    vc0 = np.array(unpack(fmt, d2))
    if ver > 0:
        v1 = vc0[0::2]
        v2 = vc0[1::2]
    else:
        v1 = vc0[0:N]
        v2 = vc0[N:(2*N)]
    ts = np.array(unpack('1q', d[0:8]))[0].astype(
        np.longlong)   # unix time, seconds
    tf = np.array(unpack('1d', d[8:16]))[0].astype(
        np.double)  # fractional second
    return v1, v2, ts, tf


def processChunk(ver, ch0, ch1, N, Nmax, Navg, nfft, w, W, M1, M2):

    # print(f'Nmax: {Nmax}')

    w = .5 * (1-np.cos(2*np.pi * np.linspace(0, 1., len(ch0), endpoint=False)))
    w = 1
    fft0 = np.fft.rfft(w*ch0, nfft)
    fft1 = np.fft.rfft(w*ch1, nfft)
    fft1 = np.ones_like(fft0)

    # Vfm1a = Vfm1[:, 0:Nmax]
    # Vfm2a = Vfm2[:, 0:Nmax]
    # Pavg1 = np.sum(np.abs(Vfm1a)**2, a/is=0) / Navg
    # Pavg2 = np.sum(np.abs(Vfm2a)**2, axis=0) / Navg

    # Correction for 0 dBm at Ettus input (Single sinusoid)  (APPROX)
    dBm0 = 0
    P1 = 10*np.log10(abs(fft0)) + dBm0
    P2 = 10*np.log10(abs(fft1)) + dBm0

    # print(Nmax)

    # Store in matrix
    M1[1:, :] = M1[:-1, :]    # Shift down
    M2[1:, :] = M2[:-1, :]
    M1[0, :] = ch0[0:Nmax]
    M2[0, :] = ch1[0:Nmax]

    """
    try:
        tm = gmtime(ts)
    except:
        exit('[python] Invalid data')
    ts = strftime('%Y-%m-%d %H:%M:', tm)
    ts += format(tm.tm_sec + tf, '06.3f') + ' UTC'
    """
    return ch0, ch1, P1, P2, M1, M2


def peakSearch(P1, P2, r0, c, B, N, nfft):
    min_r = 25
    max_r = 500
    n1 = np.floor((min_r+r0)/(c/(2.0*B) * N*1.0/nfft))
    n2 = np.floor((max_r+r0)/(c/(2.0*B) * N*1.0/nfft))
    m1 = np.max(P1[n1:n2])
    m2 = np.max(P2[n1:n2])
    print('Peaks [dBm]: ', format(m1, '6.2f'), ' , ', format(m2, '6.2f'))
    return m1, m2


def read_and_reshape(fileName, header, skip_samples=0):

    prf, tp, Nrx, fsrx = parse_params('options.txt')

    # print(fileName)
    fir = np.fromfile(fileName, np.int16, offset=0)
    # non_fir = np.fromfile(fileName, np.int16, offset=0)[1::2]
    # window_meta = np.fromfile(fileName, np.uint16, offset=0)[1::2]
    # window_labels = window_meta % 2

    #start_idx = (np.diff(window_labels) == 1).nonzero()[
    #    0]  # Use rising edges of window_labels

    n_samp = int(Nrx)
    # print(n_samp)

    # Check if reshaping is needed
    if np.size(fir) % n_samp != 0:
        # Calculate the number of zeros to pad
        zeros_to_add = n_samp - (np.size(fir) % n_samp)
        # Pad the array with zeros
        fir = np.concatenate((fir, np.zeros(zeros_to_add)))

    reshaped_fir = fir.reshape((-1, n_samp))
    # reshaped_non_fir = non_fir[start_idx[0]:start_idx[-1]].reshape((-1, n_samp))

    # print(reshaped_fir.shape)

    reshaped_fir = np.delete(reshaped_fir, np.s_[:header], 1)
    # reshaped_non_fir = np.delete(reshaped_non_fir, np.s_[:header], 1)

    # start_idx = next(x for x, val in enumerate(reshaped_fir[0])if val > 160)
    reshaped_fir = np.delete(reshaped_fir, np.s_[:skip_samples], 1)
    # reshaped_non_fir = np.delete(reshaped_non_fir, np.s_[:skip_samples], 1)

    offset_fir = np.mean(np.mean(reshaped_fir, axis=0))
    # offset_non_fir = np.mean(np.mean(reshaped_non_fir, axis=0))

    n_samp = reshaped_fir.shape[1]
    # print(n_samp)

    # , (reshaped_non_fir - offset_non_fir)
    return n_samp, (reshaped_fir - offset_fir)/16


def grouped_avg(myArray, N=2):
    cum = np.cumsum(myArray, 0)
    result = cum[N-1::N]/float(N)
    result[1:] = result[1:] - result[:-1]

    remainder = myArray.shape[0] % N
    if remainder != 0:
        if remainder < myArray.shape[0]:
            lastAvg = (cum[-1]-cum[-1-remainder])/float(remainder)
        else:
            lastAvg = cum[-1]/float(remainder)
        result = np.vstack([result, lastAvg])

    return result
