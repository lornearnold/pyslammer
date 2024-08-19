from pathlib import Path
import os
import tkinter.filedialog as tkf
import tkinter.simpledialog as dlg
import numpy as np
import datetime as dtm
import csv

G_EARTH = 9.80665

__all__ = ['csv_time_hist', 'sample_ground_motions']


def select_csv():
    """
    Open a file dialog to select a time history CSV file.

    Returns:
    file: The selected file.
    """
    file = tkf.askopenfilenames(title='Select a Time History CSV File')
    return file[0]

def sample_ground_motions():
    sgms = {}

    # Get the path to the sample_ground_motions folder
    folder_path = Path(__file__).resolve().parent / "sample_ground_motions"

    # Iterate over all files in the folder
    for file_path in folder_path.glob("*.csv"):
        # Add the file name to the list
        sgms[file_path.name[:-4]] = csv_time_hist(file_path)

    return sgms


def csv_time_hist(filename: str):
    """
    Read a CSV file containing time history data and return it as a numpy array.

    Returns:
        numpy.ndarray: A 2D numpy array containing time history data.
            The first row represents time values (seconds), and the second row represents acceleration values (multiples of g).
    """
    file = open(filename, 'r')
    if file is None:
        return None
    else:
        pass
    reader = csv.reader(file)
    time = []
    accel = []
    for row in reader:
        if '#' in row[0]:
            continue
        else:
            pass
        if len(row) == 2:
            time.append(float((row[0])))
            accel.append(float((row[1])))
        else:
            accel.append(float((row[0])))
    if len(time) == 0:
        dt = dlg.askfloat('Enter dt', 'Enter a time interval (s): ')
        time = np.arange(0, len(accel)*dt, dt)
    else:
        pass
    time = np.array(time)
    accel = np.array(accel)
    time_history = np.vstack((time, accel))
    return time_history


# def test_time_hist(freq_hz: float, write=False, output_dir=''):
#     """
#     Generate a 30 second time history of acceleration for a given frequency. 
#     Optionally write the time history to a CSV file.

#     Parameters:
#     freq_hz (float): The frequency in hertz.
#     write (bool, optional): Flag to indicate whether to write the time history to a file. Defaults to False.
#     output_dir (str, optional): The output directory path. Defaults to ''.

#     Returns:
#     numpy.ndarray: The time history array.

#     """     
#     time = np.arange(0, 30, 0.001)
#     accel = np.sin(freq_hz*2*np.pi*time) * G_EARTH
#     time_history = np.vstack((time, accel))
#     if write == True:
#         if output_dir == '':
#             output_dir = tkf.askdirectory(title='Select an Output Directory')
#         else:
#             pass
#         output_name = '/' + str(freq_hz) + '_Hz_harmonic_' + dtm.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.csv'
#         output = output_dir + output_name
#         with open(output, 'w', newline='') as file:
#             writer = csv.writer(file)
#             writer.writerows(time_history.T)
#     return time_history


# def write_output(block_data: np.ndarray, output_name='', output_dir='', input_file=''):
#     """
#     Write the block data to a CSV file.

#     Args:
#         block_data (numpy.ndarray): The block data to be written.
#         output_name (str, optional): The name of the output file. If not provided, a timestamp will be used.
#         output_dir (str, optional): The directory where the output file will be saved. If not provided, a directory will be selected using a file dialog.
#         input_file (str, optional): The input file path. If provided, the output directory and name will be derived from the input file.

#     Returns:
#         None
#     """
#     if input_file != '':
#         output_dir = Path(input_file).parent
#         output_name = Path(input_file).stem + '_output.csv'
#     else:
#         if output_dir == '':
#             output_dir = tkf.askdirectory(title='Select an Output Directory')
#             if output_dir == '':
#                 return
#         else:
#             pass
#         if output_name == '':
#             output_name = 'pySLAMMER_' + dtm.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '_output.csv'
#         else:
#             output_name = output_name + '_' + dtm.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.csv'
#     output = output_dir / output_name
#     with open(output, 'w', newline='') as file:
#         writer = csv.writer(file)
#         writer.writerows(block_data.T)