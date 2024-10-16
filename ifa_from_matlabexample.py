# -*- coding: utf-8 -*-
"""
IFA Antenna Simulation using openEMS and Python

Converted from MATLAB code

Tested with
 - python 3.x
 - openEMS v0.0.34+

"""

### Import Libraries
import os
root = r"C:\opt\openEMS"
os.add_dll_directory(root)

import numpy as np
import tempfile

from pylab import *
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

#############################################################################
#                substrate_width
#  _______________________________________________    __ substrate_thickness
# | A                        ifa_l                |\  __
# | |ifa_e         __________________________     | |
# | |             |    ___  _________________| w2 | |
# | |       ifa_h |   |   ||                      | |
# |_V_____________|___|___||______________________| |
# |                _w1   _wf\                     | |
# |                   |_fp|  \                    | |
# |                       |    feed point         | |
# |                       |                       | | substrate_length
# |<- substrate_width/2 ->|                       | |
# |                                               | |
# |_______________________________________________| |
#  \_______________________________________________\|
#
# Note: It's not checked whether your settings make sense, so check
#       graphical output carefully.
#
##############################################################################

c0 = C0
### Simulation parameters
unit = 1e-3  # all lengths in mm

# Substrate parameters
substrate_width = 80
substrate_length = 80
substrate_thickness = 1.5
substrate_cells = 4

# IFA parameters
ifa_h = 8       # height of short circuit stub
ifa_l = 22.5    # length of radiating element
ifa_w1 = 4      # width of short circuit stub
ifa_w2 = 2.5    # width of radiating element
ifa_wf = 1      # width of feed element
ifa_fp = 4      # position of feed element relative to short circuit stub
ifa_e = 10      # distance to edge

# Substrate material properties
substrate_epsR = 4.5
substrate_kappa = 1e-3 * 2 * pi * 2.45e9 * EPS0 * substrate_epsR

# Feeding setup
feed_R = 50

#minmax frequency
min_freq = 2.4e9
center_freq = 2.45e9
max_freq = 2.5e9

# Show structure
showCad = True

# Simulation box size
SimBox = np.array([substrate_width * 2, substrate_length * 2, 150])

# FDTD parameters
f0 = 2.5e9  # center frequency
fc = 0.5e9    # 20 dB corner frequency

FDTD = openEMS(NrTS=60000)#, EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

# Initialize CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

# Initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
mesh.AddLine('z', [-SimBox[2]/2, SimBox[2]/2])

# Create substrate
substrate = CSX.AddMaterial('substrate', epsilon=substrate_epsR, kappa=substrate_kappa)
start = [-substrate_width/2, -substrate_length/2, 0]
stop = [substrate_width/2, substrate_length/2, substrate_thickness]
substrate.AddBox(start=start, stop=stop, priority=1)

# Add extra cells to discretize the substrate thickness
mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells+1))

# Create ground plane
gnd = CSX.AddMetal('groundplane')
start = [-substrate_width/2, -substrate_length/2, substrate_thickness]
stop = [substrate_width/2, substrate_length/2 - ifa_e, substrate_thickness]
gnd.AddBox(start=start, stop=stop, priority=10)

FDTD.AddEdges2Grid(dirs='xy', properties=gnd)

# Create IFA
ifa = CSX.AddMetal('ifa')
tl = np.array([0, substrate_length / 2 - ifa_e, substrate_thickness])  # translation vector

# Feed element
start = np.array([0, 0.5, 0]) + tl
stop = start + np.array([ifa_wf, ifa_h - 0.5, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

# Short circuit stub
start = np.array([-ifa_fp, 0, 0]) + tl
stop = start + np.array([-ifa_w1, ifa_h, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

# Radiating element
start = np.array([-ifa_fp - ifa_w1, ifa_h, 0]) + tl
stop = start + np.array([ifa_l, -ifa_w2, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

# Add IFA edges to grid
FDTD.AddEdges2Grid(dirs='xy', properties=ifa, metal_edge_res=0.5)

# Apply the excitation & resistor as a current source
start = np.array([0, 0, 0]) + tl
stop = start + np.array([ifa_wf, 0.5, 0])
#port = FDTD.AddLumpedPort(5, feed_R, start, stop, 'y', 1.0, priority=5)
#port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'z', 1.0, priority=5, edges2grid='xy')
port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'y', 1.0, priority=5, edges2grid='xz')
# Finalize the mesh
# Generate a smooth mesh with max. cell size: lambda_min / 20
mesh_res = c0 / (f0 + fc) / unit / 20
mesh.SmoothMeshLines('all', mesh_res, 1.4)

# Add the nf2ff recording box
nf2ff = FDTD.CreateNF2FFBox()

# Prepare simulation folder
Sim_Path = os.path.abspath('tmp_MATLAB_IFA')
Sim_CSX = 'IFA.xml'
post_proc_only = True  # Set to True to skip simulation

if os.path.exists(Sim_Path):
    import shutil
    shutil.rmtree(Sim_Path)
os.mkdir(Sim_Path)
CSX_file = os.path.join(Sim_Path, Sim_CSX)
CSX.Write2XML(CSX_file)

# Show the structure
if showCad:
    from CSXCAD import AppCSXCAD_BIN
    os.system(root+'\\'+AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

# Run openEMS
if not post_proc_only:
    try:
        #FDTD.Run(Sim_Path, verbose=3, cleanup=True)
        FDTD.Run(Sim_Path, verbose=3, cleanup=True, debug=True)
    except Exception as e:
        print("An error occurred during simulation:", e)
        import traceback
        traceback.print_exc()

# Post-processing & plotting
freq = np.linspace(max(1e9, f0 - fc), f0 + fc, 501)
port.CalcPort(Sim_Path, freq)

Zin = port.uf_tot / port.if_tot
s11 = port.uf_ref / port.uf_inc
s11_dB = 20.0*np.log10(np.abs(s11))
P_in = np.real(0.5 * port.uf_tot * np.conj(port.if_tot))

# Create a figure with subplots, 2 rows and 2 columns
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot feed point impedance in the top-left subplot (axs[0, 0])
axs[0, 0].plot(freq / 1e6, np.real(Zin), 'k-', linewidth=2, label='Re{Zin}')
axs[0, 0].plot(freq / 1e6, np.imag(Zin), 'r--', linewidth=2, label='Im{Zin}')
# Add a horizontal line at 50 Ohms
axs[0, 0].axhline(feed_R, color='blue', linestyle='--', linewidth=1, label=f'{feed_R} Ohm')
axs[0, 0].axvline(min_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{min_freq/1e6} MHz')
axs[0, 0].axvline(center_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{center_freq/1e6} MHz')
axs[0, 0].axvline(max_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{max_freq/1e6} MHz')


# Customize the grid, title, labels, and legend
axs[0, 0].grid()
axs[0, 0].set_title('Feed point impedance')
axs[0, 0].set_xlabel('Frequency f / MHz')
axs[0, 0].set_ylabel('Impedance Zin / Ohm')
axs[0, 0].legend()


# Plot reflection coefficient S11 in the top-right subplot (axs[0, 1])
axs[0, 1].plot(freq / 1e6, 20 * np.log10(np.abs(s11)), 'k-', linewidth=2,label='S11(db)')

cutoffDB = -6
axs[0, 1].axhline(cutoffDB, color='blue', linestyle='--', linewidth=1, label=f'{cutoffDB} dB Cutoff')
axs[0, 1].axvline(min_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{min_freq/1e6} MHz')
axs[0, 1].axvline(center_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{center_freq/1e6} MHz')
axs[0, 1].axvline(max_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{max_freq/1e6} MHz')

axs[0, 1].grid()
axs[0, 1].set_title('Reflection coefficient S11')
axs[0, 1].set_xlabel('Frequency f / MHz')
axs[0, 1].set_ylabel('Reflection coefficient |S11| (dB)')
axs[0, 1].legend()
# For polar subplots, calculate NF2FF if the resonance frequency is found
thetaRange = np.arange(0, 182, 2)
phiRange = np.arange(0, 362, 2) - 180

idx = np.where((s11_dB < -10) & (s11_dB == np.min(s11_dB)))[0]
if len(idx) != 1:
    print('No resonance frequency found for far-field calculation')
else:
    # Find resonance frequency from s11
    f_res_ind = np.argmin(np.abs(s11))
    f_res = freq[f_res_ind]

    theta = np.arange(-180.0, 180.0, 2.0)
    print("Calculate NF2FF")
    nf2ff_res_phi0 = nf2ff.CalcNF2FF(Sim_Path, f_res, theta, 0, verbose=1, outfile='3D_Pattern.h5')

    # Create the bottom-left polar subplot (axs[1, 0]) for the xz-plane
    ax_polar1 = plt.subplot(223, polar=True)
    E_norm = 20.0 * np.log10(nf2ff_res_phi0.E_norm / np.max(nf2ff_res_phi0.E_norm)) + nf2ff_res_phi0.Dmax
    ax_polar1.plot(np.deg2rad(theta), 10 ** (np.squeeze(E_norm) / 20), linewidth=2, label='xz-plane')
    ax_polar1.grid(True)
    ax_polar1.set_xlabel('theta (deg)')
    ax_polar1.set_theta_zero_location('N')
    ax_polar1.set_theta_direction(-1)
    ax_polar1.legend(loc=3)

    # Create the bottom-right polar subplot (axs[1, 1]) for the xy-plane
    phi = theta
    nf2ff_res_theta90 = nf2ff.CalcNF2FF(Sim_Path, f_res, 90, phi, center=np.array([0, 0, 0]) * unit, read_cached=True, outfile='nf2ff_xy.h5')

    ax_polar2 = plt.subplot(224, polar=True)
    E_norm = 20.0 * np.log10(nf2ff_res_theta90.E_norm / np.max(nf2ff_res_theta90.E_norm)) + nf2ff_res_theta90.Dmax
    ax_polar2.plot(np.deg2rad(phi), 10 ** (np.squeeze(E_norm) / 20), linewidth=2, label='xy-plane')
    ax_polar2.grid(True)
    ax_polar2.set_xlabel('phi (deg)')
    plt.suptitle('Bent Patch Antenna Pattern\nFrequency: {:.2f} GHz'.format(f_res / 1e9), fontsize=14)
    ax_polar2.legend(loc=3)

    print('Radiated power: Prad = {:.2e} Watt'.format(nf2ff_res_theta90.Prad[0]))
    print('Directivity:    Dmax = {:.1f} ({:.1f} dBi)'.format(nf2ff_res_theta90.Dmax[0], 10 * np.log10(nf2ff_res_theta90.Dmax[0])))
    print('Efficiency:   nu_rad = {:.1f} %'.format(100 * nf2ff_res_theta90.Prad[0] / np.real(P_in[idx[0]])))

# Adjust the layout to ensure everything fits
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout and add space for suptitle

# Show all plots
plt.show()