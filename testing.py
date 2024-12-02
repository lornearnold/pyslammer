import csv
import numpy as np
import datetime as dtm
import matplotlib as mpl
import tkinter.filedialog as tkf
from pathlib import Path
from matplotlib import pyplot as plt
from pyslammer.rigid_block import RigidBlock


G_EARTH = 9.80665

class SlammerData:

    def __init__(self):
        self.name = None
        self.station = None
        self.rigid_normal = None
        self.rigid_inverse = None
        self.rigid_average = None
        self.k_y = None

    def __str__(self):
        return self.name
    

class PySlammerData:

    def __init__(self):
        self.name = None
        self.station = None
        self.rigid_normal = None
        self.rigid_inverse = None
        self.rigid_average = None
        self.k_y = None

    def __str__(self):
        return self.name
    
    def average(self):
        self.rigid_average = (self.rigid_normal + self.rigid_inverse) / 2


def plot_classic(data):
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 8), sharex=True)
    fig.suptitle('Rigid Block Analysis for Loma Prieta LGP-000')
    ax[0].plot(data.time, data.gnd_acc/G_EARTH, label='Ground Acceleration')
    ax[0].set_ylabel('Ground \nAcceleration \n(g)')
    ax[1].plot(data.time, data.block_vel, label='Block Velocity')
    ax[1].set_ylabel('Block \nVelocity \n(m/s)')
    ax[2].plot(data.time, data.block_disp, label='Block Displacement')
    ax[2].set_ylabel('Block \nDisplacement \n(m)')
    ax[-1].set_xlabel('Time (s)')
    #vel = np.array(data.block_vel)
    #ax[0].fill_between(data.time, min(data.gnd_acc/G_EARTH), max(data.gnd_acc/G_EARTH), where=vel > 0, facecolor='gray', alpha=0.5)
    #ax[1].fill_between(data.time, 0, max(data.block_vel), where=vel > 0, facecolor='gray', alpha=0.5)
    #ax[2].fill_between(data.time, 0, max(data.block_disp), where=vel > 0, facecolor='gray', alpha=0.5)
    plt.subplots_adjust(left=0.15)
    plt.show()


def csv_time_hist(filename: str):
    """
    Read a CSV file containing time history acceleration data and return a 1D numpy array and a timestep

    Returns:
        data: A 2D numpy array containing time and acceleration data
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
    data = np.vstack((time, accel))
    return data


def analytical_test(src_dir: str, write: bool=False, manual: bool=False, method: str=None):
    # Open the harmonic solutions file and read the data.
    file = open(src_dir)
    if file is None:
        exit()
    reader = csv.reader(file)
    freq = []
    k_y = []
    analytical_displacement = []
    # ky and frequency values that have been analyzed.
    ky_analyzed = []
    freq_analyzed = []
    # SLAMMER displacements for each frequency and ky value.
    slammer_disp = []
    for row in reader:
        if '#' in row[0]:
            # Ignore headers
            continue
        freq.append(float(row[0]))
        if freq[-1] not in freq_analyzed:
            freq_analyzed.append(freq[-1])
        k_y.append(float(row[1]))
        if k_y[-1] not in ky_analyzed:
            ky_analyzed.append(k_y[-1])
        analytical_displacement.append(float(row[2]))
        if manual:
            slammer_disp.append(float(row[8]))
    file.close()

    # Calculate the displacement for each frequency and ky value.
    error = []
    if manual:
        for i in range(len(slammer_disp)):
            error.append(((slammer_disp[i] - analytical_displacement[i])/analytical_displacement[i]) * 100)
    else:
        disp = []
        if write:
            output_dir = tkf.askdirectory(title='Select an Output Directory')
        for i in range(len(freq)):
            # Create time and acceleration arrays using each frequency and ky value
            # for one half of the period.
            time = np.arange(0, 1/(freq[i]), 0.001)
            accel = np.sin(freq[i]*2*np.pi*time)
            time_history = np.vstack((time, accel))

            # Write the time history to a CSV file if the write flag is set.
            if write:
                output_name = '/' + str(freq[i]) + '_Hz_harmonic_' + dtm.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.csv'
                output = output_dir + output_name
                with open(output, 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerows(time_history.T)

            # Calculate the displacement using the desired method and append the total
            # displacement to the disp list. User to modify the next three lines as
            # needed to test the desired method.
            data = RigidBlock(time_history)
            if method == 'jibson':
                data.downslope_jibson(k_y[i])
            elif method == 'dgr':
                data.downslope_dgr(k_y[i])
            disp.append(data.block_disp[-1])

        # Calculate the percent error between the analytical solution and the calculated.
        error.append(((data.block_disp[-1] - analytical_displacement[i])/analytical_displacement[i]) * 100)

    # Create a data array to hold frequency and error data (dimension 1) sorted by frequency 
    # tested (dimension 1) and ky value (dimension 2).
    data = np.ndarray(shape=(2, len(freq_analyzed), len(ky_analyzed)))
    ky_idx = [0] * len(ky_analyzed)

    i = 0
    for i in range(len(freq)):
        # Find the index of the ky value in the ky_analyzed list to save frequency and error 
        # data in the correct location in the data array based on ky.
        idx = ky_analyzed.index(float(k_y[i]))
        data[0][ky_idx[idx]][idx] = float(freq[i])
        data[1][ky_idx[idx]][idx] = float(error[i])
        ky_idx[idx] += 1
        i += 1

    # Plot the error vs frequency data for each ky value.
    markers = ['o', 's', 'D', '^', 'v']
    fig, ax = plt.subplots()
    for j in range(len(ky_analyzed)):
        ax.plot(data[0, :, j], data[1, :, j], markers[j], label='$\mathregular{k_y}$ = ' + str(ky_analyzed[j]))
    if manual:
        ax.set_title('SLAMMER Analytical Solution Test')
    if method == 'jibson':
        ax.set_title('Relative Acceleration Method Analytical Solution Test')
    elif method == 'dgr':
        ax.set_title('Relative Velocity Method Analytical Solution Test')
    else:
        pass
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Error (%)')
    ax.set_ylim(-0.7, 0.005)
    ax.set_xscale('log')
    ax.legend()
    plt.show()


def ground_motion_comp(slam_dir: str, data_dir: str):
    # Open the pySLAMMER results file and read the data.
    file = open(slam_dir)
    if file is None:
        exit()
    reader = csv.reader(file)
    slammer_results = []
    names = []
    for row in reader:
        if '#' in row[0]:
            continue
        if row[0] not in names:
            data = SlammerData()
            name = row[0]
            names.append(name)
            # Remove commas, periods, and spaces from the name to create a valid file name.
            name = name.translate(str.maketrans({ord(','): None, ord('.'): None, ord(' '): '_'}))           
            data.station = row[1]
            data.name = name + '_' + data.station
            data.rigid_normal = float(row[2])/100
            data.rigid_inverse = float(row[3])/100
            data.rigid_average = float(row[4])/100
            data.k_y = float(row[15])
            slammer_results.append(data)
        else:
            continue
    file.close()
    pyslammer_results = [None] * len(slammer_results)
    data_dir = Path(data_dir)
    
    for file in data_dir.iterdir():
        # Iterate over each record in the data directory and create a RigidBlock object for each.
        filename = str(file.stem)
        # Find the index of the slammer results that corresponds to the record.
        idx = [i for i, x in enumerate(slammer_results) if x.name == filename]
        # Calculate the downslope displacement using the ky value from the pySLAMMER results.
        data = PySlammerData()
        sba = RigidBlock(csv_time_hist(file), name=filename)
        data.name = filename
        data.station = slammer_results[idx[0]].station
        # Calculate normal dispalcement.
        sba.downslope_jibson(slammer_results[idx[0]].k_y)
        data.rigid_normal = sba.total_disp
        # Calculate inverse displacement.
        sba.invert()
        sba.downslope_jibson(slammer_results[idx[0]].k_y)
        data.rigid_inverse = sba.total_disp
        # Calculate average displacement.
        data.average()
        pyslammer_results[idx[0]] = data

    ref_data = []
    pyslam_data_normal = []
    pyslam_data_inverse = []
    pyslam_data_average = []
    i = 0
    for i in range(len(slammer_results)):
        ref_data.append(slammer_results[i].rigid_normal/slammer_results[i].rigid_normal)
        pyslam_data_normal.append(pyslammer_results[i].rigid_normal/slammer_results[i].rigid_normal)
        pyslam_data_inverse.append(pyslammer_results[i].rigid_inverse/slammer_results[i].rigid_inverse)
        pyslam_data_average.append(pyslammer_results[i].rigid_average/slammer_results[i].rigid_average)
        i += 1
        
    
    plt.plot(names, ref_data, label='SLAMMER (All)')
    plt.plot(names, pyslam_data_normal, 'o', label='pySLAMMER Normal')
    plt.plot(names, pyslam_data_inverse, 's', label='pySLAMMER Inverse')
    plt.plot(names, pyslam_data_average, 'D', label='pySLAMMER Average')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.title('Rigid Block Displacement Comparison Results')
    plt.ylabel('Displacement Ratio (d/$\mathregular{d_{SLAMMER}}$)')
    plt.xlabel('Earthquake Record')
    plt.legend()
    plt.show()

    pass


harmonic_solutions = "tests/Trimmed Data/harmonic_errors.csv"
loma_prieta = csv_time_hist("tests/Trimmed Data/Loma_Prieta_1989_HSP-000.csv")

mpl.rcParams.update({'font.family': 'times new roman', 'font.size': 12})
analytical_test(harmonic_solutions, manual=True)

mpl.rcParams.update({'font.family': 'times new roman', 'font.size': 20})
data = RigidBlock(loma_prieta)
data.downslope_jibson(0.25)
#plot_classic(data)

slammer_results = 'test/Unmodified Data/pySLAMMER_motion_sample_results.csv'
data_location = 'tests/Ground Motion Test Set'
# ground_motion_comp(slammer_results, data_location)