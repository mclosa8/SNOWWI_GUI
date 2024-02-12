#!/usr/bin/env python3

"""
Real-time acquisition program for SNOWWI.
Runs the C++ data acquisition program while rendering plots & images in real time.

Usage:
    python3 snowwi-daq.py

Note:
    The C++ daq program can also be run stand-alone as:
        TODO: modify this

Marc Closa Tarres       August 2023
(mclosatarres@umass.edu)
"""

from snowwi_lib import *
from signal import SIGINT, SIGKILL, SIGTERM
from traceback import format_tb

import argparse
import atexit
import datetime
import os
import psutil
import shutil
import subprocess
import sys
import time

import datetime

output_date = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
print(output_date, end='\n')

# Add and parse args
parser = argparse.ArgumentParser()

# parser.add_argument('script')
# parser.add_argument('-N', '--num_samps', default=10000)
# parser.add_argument('-t', '--run_time', required=True)
parser.add_argument('-o', '--out_dir', required=False, default=output_date)
parser.add_argument('-r', '--max_range', default=None)
parser.add_argument('-nh', '--header', required=False,
                    default=52, type=int,help='Number of header samples')
parser.add_argument('-c', '--channel', required=False, type=int, default=0)
parser.add_argument('-nvme', '--nvme', required=False, default=2, type=int)

args = parser.parse_args()

# assert args.nvme == 1 or args.nvme == 0, "-nvme arg should be 0 or 1"

# Data directory
# TODO: Replace with proper path
root_path = '/home/frosty/Desktop/replay_trigger_Colorado_FINAL/build'
# data_path = os.path.join(root_path, f'save_data_nvme{args.nvme}', args.out_dir)
data_path = os.path.join(root_path, args.out_dir)
channel_path = os.path.join(data_path, f'chan{args.channel}')

save_path = os.path.join(root_path, f'save_data_nvme{args.nvme}')

# Create data directory
if (os.path.isdir(data_path)):
    print(f'Data directory {data_path} already exists.\n')
else:
    print(f'Creating data directory: {data_path}\n')
    os.mkdir(data_path)
    # Create directiories for channels
    os.mkdir(os.path.join(data_path, 'chan0'))
    os.mkdir(os.path.join(data_path, 'chan1'))
    os.mkdir(os.path.join(data_path, 'chan2'))
    os.mkdir(os.path.join(data_path, 'chan3'))

# Print paths
print('Working path:')
print(f'    {root_path}\n')

print('Data path:')
print(f'    {data_path}\n')

print('Channel path:')
print(f'    {channel_path}\n')

# Parse params from setup txt file
prf, tp, Nrx, fsrx = parse_params('options.txt')  # TODO: adapt to SNOWWI

# Script name
script_name = f"""
    ./bash_replay_trigger.sh options.txt {args.out_dir}
"""

print(f'Running script: {script_name}...')

# Pulses to integrate
N_avg = 25
print(f'Averaging {N_avg} chirps...')

# Max range to display
if (args.max_range is not None):
    r_max = args.max_range
    max_samp = 2 * r_max * fsrx / 3e8
    print(f'Maximum range displayed {r_max} m.')
else:
    max_samp = Nrx
    r_max = max_samp / fsrx * 3e8 / 2
    print(f'Maximum range displayed {r_max} m.')

# Define aux functions


def stop_daq():
    """
    Stops C++ acq program interrupting the subprocess.

    See: https://docs.python.org/3/library/subprocess.html
    and  https://docs.python.org/3/library/signal.html
    """
    try:
        os.system('pkill -SIGINT -f replay_trigger')
        print('C++ daq program interrupted by Python program.')
        time.sleep(1)
        print('Moving files to NVME...')
        shutil.move(data_path, save_path)
        shutil.copyfile(os.path.join(root_path, 'options.txt'),
                        os.path.join(save_path, args.out_dir, 'options.txt'))
        shutil.move(os.path.join(root_path, 'buffer_errors.txt'),
                    os.path.join(save_path, args.out_dir))
        shutil.move(os.path.join(root_path, 'terminal.txt'),
                    os.path.join(save_path, args.out_dir))
        shutil.copyfile(os.path.join(root_path, 'chirp_params.txt'),
                        os.path.join(save_path, args.out_dir, 'chirp_params.txt'))
        shutil.copyfile(os.path.join(root_path, 'chirp.csv'),
                        os.path.join(save_path, args.out_dir, 'chirp.csv'))
        shutil.copyfile(os.path.join(root_path, 'send.dat'),
                        os.path.join(save_path, args.out_dir, 'send.dat'))
        
        print(f'{data_path} moved to {save_path}')
    except:
        print('Error stopping daq process')


# Clean up any process related to daq
subprocess.call(['pkill', script_name.split()[0][2:]])

print('Killing alive processes...')
time.sleep(3)

# Configure handlers for exit - Function to be executed at exit
atexit.register(stop_daq)

def cleanup(exectype, value, tb):
    print('Exception \n', format_tb(tb), exectype, '; ', value)
    stop_daq()

sys.excepthook = cleanup

# Start the C++ program and wait for data
fname0 = second_to_last(channel_path)
fname = fname0  # Here it will be None
print('FILENAME', channel_path, fname, end='\n')

acq_process = subprocess.Popen(script_name.split())
while (fname == fname0):  # Loops until the fname updates
    time.sleep(0.1)
    fname = second_to_last(channel_path)
print('FILENAME', channel_path, fname, end='\n')


ver = 1


# Hardcoded values needed from acq params
tp = 1e-3
prf = 1/tp
fsrx = 491.52e6*2

"""
    TODO:
        - Number of samples from file
        - Sampling frequency
        - Prf / prt
"""

# Sets the constants prior to data reading
c, B, r0, N, n, nfft, Nf, Nmax, rmax, x, f, r, w, W, M1, M2 = prep(
    fsrx, tp, Nrx, r_max, ver, N_avg, args.header)  # TODO: adapt to SNOWWI

# Set up figure environment
fig, i_plot, q_plot, r0_plot, r1_plot, im0_plot, im1_plot, title1, hist0, hist1, Iraw, Qraw = setupFig(
    N, r0, r_max, M1, M2, fsrx)

# Define button and button event


def button_func(event):
    exit(0)


ax = plt.axes([0.9, 0.02, 0.06, 0.03])
b = Button(ax, 'Stop')
b.on_clicked(button_func)

npk = -1  # Not sure why is this for


while plt.fignum_exists(fig.number):

    # Gets second to last file
    fname = second_to_last(channel_path)

    n_samp, fir = read_and_reshape(fname, args.header)
    n_rows = fir.shape[0]
    n_samps = fir.shape[1]

    # print(f'Fir: {fir.shape}')

    # N rows after average
    n_steps = np.ceil(fir.shape[0]/N_avg)

    # Average N_avg rows
    avgd0 = grouped_avg(fir, N_avg)
    avgd1 = np.ones_like(avgd0)

    first_pulses0 = fir[3::N_avg]
    first_pulses1 = np.ones_like(first_pulses0)
    # print(avgd0.shape, avgd1.shape)

    for i, row in enumerate(avgd0):
        ch0 = row
        ch1 = avgd1[i]

        # Here you process the chunk
        # TODO: change process chunk
        v1, v2, P1, P2, M1, M2 = processChunk(
            ver, ch0, ch1, N, Nmax, N_avg, nfft, w, W, M1, M2)

        Iraw.set_data(n, first_pulses0[i])  # Raw data for Ch0
        Qraw.set_data(n, first_pulses1[i])

        i_plot.set_data(n, v1)  # Bit value for Ch0
        q_plot.set_data(n, v2)  # Bit value for Ch1

        r0_plot.set_data(f, P1)  # Range plot for Ch0
        r1_plot.set_data(f, P2)  # Range plot for Ch1

        im0_plot.set_data(M1)  # Heatmap plot for Ch0
        im1_plot.set_data(M2)  # Heatmap plot for Ch1

        data_for_hist0 = np.zeros_like(v1)
        data_for_hist1 = np.zeros_like(v2)

        # hist0, bin_edges0 = np.histogram(data_for_hist0, bins='auto')
        # hist1, bin_edges1 = np.histogram(data_for_hist1, bins='auto')

        plt.pause(0.001)
        time.sleep(tp*N_avg)

        # title1.set_text(ts)

# Stop acquisition if the figure is closed from the close button (X)
stop_daq()
