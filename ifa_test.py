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

from datetime import datetime as dt
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
# | <---ifa_e---->| w1   _wf\                     | |
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

# Get the current timestamp
timestamp = dt.now().strftime('%Y%m%d_%H%M%S')

# Add the timestamp to the file name
Sim_Path = os.path.abspath(f'runs\\tmp_IFA_{timestamp}')
Sim_CSX = 'IFA.xml'

# Show structure
showCad = False #set true to show cad
post_proc_only = False  # Set to True to skip simulation

c0 = C0
### Simulation parameters
unit = 1e-3  # all lengths in mm

# Substrate parameters
substrate_width = 21
substrate_length = 83.15 #not including antenna element
substrate_thickness = 1.5
gndplane_position=0#referenced from top layer
substrate_cells = 4

# IFA parameters
ifa_h = 5.586       # height of short circuit stub
ifa_l = 16.6      # length of radiating element
ifa_w1 =0.613      # width of short circuit stub
ifa_w2 =0.472      # width of radiating element
ifa_wf =0.425      # width of feed element
ifa_fp =1.108      # position of feed element relative to short circuit stub
ifa_e = 0.5      # distance to edge

ifa_h+ifa_e
# Substrate material properties
substrate_epsR = 4.5
substrate_kappa = 1e-3 * 2 * pi * 2.45e9 * EPS0 * substrate_epsR

# Feeding setup
feed_R = 50

#minmax frequency
min_freq = 2.4e9
center_freq = 2.45e9
max_freq = 2.5e9



# Simulation box size
SimBox = np.array([substrate_width * 3, substrate_length * 2, 150])

# FDTD parameters
f0 = 2.45e9  # center frequency
fc = 1.0e9    # 20 dB corner frequency

FDTD = openEMS(NrTS=600000)#, EndCriteria=1e-5)
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
stop = [substrate_width/2, substrate_length/2+ifa_h+ifa_e, substrate_thickness]
substrate.AddBox(start=start, stop=stop, priority=1)

# Add extra cells to discretize the substrate thickness
mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells+1))

# Create ground plane
gnd = CSX.AddMetal('groundplane')
start = [-substrate_width/2+ifa_e, -substrate_length/2+ifa_e, substrate_thickness+gndplane_position]
stop = [substrate_width/2-ifa_e, substrate_length/2 - ifa_e, substrate_thickness+gndplane_position]
gnd.AddBox(start=start, stop=stop, priority=10)

FDTD.AddEdges2Grid(dirs='xy', properties=gnd)

# Create IFA
ifa = CSX.AddMetal('ifa')
tl = np.array([-substrate_width/2+ifa_e+ifa_fp+ifa_w1, substrate_length / 2 - ifa_e, substrate_thickness])  # translation vector


def extend_line(start, stop):
    # Calculate the total range for each dimension
    total_range = [stop[i] - start[i] for i in range(3)]

    # Calculate 1/6 of the range for each dimension
    one_sixth = [r / 6.0 for r in total_range]

    # Calculate new start and stop points for inside lines (2/3 inside the figure)
    inside_start = [start[i] + one_sixth[i] for i in range(3)]
    inside_end = [stop[i] - one_sixth[i] for i in range(3)]

    # Calculate the outside extensions (1/6 beyond start and stop points)
    outside_start = [start[i] - one_sixth[i] for i in range(3)]
    outside_end = [stop[i] + one_sixth[i] for i in range(3)]

    # Construct the coordinate arrays
    x = [outside_start[0], inside_start[0], inside_end[0], outside_end[0]]
    y = [outside_start[1], inside_start[1], inside_end[1], outside_end[1]]
    z = [outside_start[2], inside_start[2], inside_end[2], outside_end[2]]

    return x, y, z
# Feed element
start = np.array([0, ifa_w2, 0]) + tl
stop = start + np.array([ifa_wf, ifa_h - ifa_w2, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

meshlines = extend_line(start, stop)
mesh.AddLine('x',meshlines[0] )
mesh.AddLine('y', meshlines[1])

# Short circuit stub
start = np.array([-ifa_fp, 0, 0]) + tl
stop = start + np.array([-ifa_w1, ifa_h, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

meshlines = extend_line(start, stop)
mesh.AddLine('x',meshlines[0] )
mesh.AddLine('y', meshlines[1])

# Radiating element
start = np.array([-ifa_fp - ifa_w1, ifa_h, 0]) + tl
stop = start + np.array([ifa_l, -ifa_w2, 0])
ifa.AddBox(start=start, stop=stop, priority=10)

meshlines = extend_line(start, stop)
mesh.AddLine('x',meshlines[0] )
mesh.AddLine('y', meshlines[1])
# Add IFA edges to grid
#FDTD.AddEdges2Grid(dirs='xy', properties=ifa, metal_edge_res=0.5)

# Apply the excitation & resistor as a current source
start = np.array([0, 0, 0]) + tl
stop = start + np.array([ifa_wf, ifa_w2, 0])
port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'y', 1.0, priority=5, edges2grid='xz')

meshlines = extend_line(start, stop)
mesh.AddLine('x',meshlines[0] )
mesh.AddLine('y', meshlines[1])

# Finalize the mesh
# Generate a smooth mesh with max. cell size: lambda_min / 20
mesh_res = c0 / (f0 + fc) / unit / 20
print(mesh_res)
mesh.SmoothMeshLines('all', mesh_res, 1.5)

# Add the nf2ff recording box
nf2ff = FDTD.CreateNF2FFBox()

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
    plt.suptitle('IFA Antenna Pattern\nFrequency: {:.2f} GHz'.format(f_res / 1e9), fontsize=14)
    ax_polar2.legend(loc=3)

    print('Radiated power: Prad = {:.2e} Watt'.format(nf2ff_res_theta90.Prad[0]))
    print('Directivity:    Dmax = {:.1f} ({:.1f} dBi)'.format(nf2ff_res_theta90.Dmax[0], 10 * np.log10(nf2ff_res_theta90.Dmax[0])))
    print('Efficiency:   nu_rad = {:.1f} %'.format(100 * nf2ff_res_theta90.Prad[0] / np.real(P_in[idx[0]])))

# Prepare the parameters text with table-like formatting
parameters_text = (
    f"{'IFA Parameters:':<30}\n"
    f"{'ifa_h:':<20}{ifa_h:>10.3f} mm "
    f"{'ifa_l:':<20}{ifa_l:>10.3f} mm\n"
    f"{'ifa_w1:':<20}{ifa_w1:>10.3f} mm "
    f"{'ifa_w2:':<20}{ifa_w2:>10.3f} mm\n"
    f"{'ifa_wf:':<20}{ifa_wf:>10.3f} mm "
    f"{'ifa_fp:':<20}{ifa_fp:>10.3f} mm\n"
    f"{'ifa_e:':<20}{ifa_e:>10.3f} mm "
    f"{'ifa_h + ifa_e:':<20}{ifa_h+ifa_e:>10.3f} mm\n\n"
    f"{'Substrate Properties:':<30}\n"
    f"{'substrate_epsR:':<20}{substrate_epsR:>10.3f} "
    f"{'substrate_kappa:':<20}{substrate_kappa:>10.3e} S/m\n\n"
    f"{'Feeding Setup:':<30}\n"
    f"{'feed_R:':<20}{feed_R:>10} Ω"
)
# Add the parameters text to the figure at the bottom
fig.text(0.5, 0.01, parameters_text, ha='center', va='bottom', fontsize=12, family='monospace')

# Adjust layout to prevent overlap with the text
plt.tight_layout(rect=[0, 0.2, 1, 0.95])  # Adjust rect to leave space at bottom and top

# Show all plots
plt.figure()
plt.show()