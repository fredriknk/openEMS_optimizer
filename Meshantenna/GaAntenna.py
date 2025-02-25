# -*- coding: utf-8 -*-
"""
IFA/MIFA Antenna Simulation using openEMS and Python

Converted from MATLAB code

Tested with
 - python 3.1x
 - openEMS v0.0.34+

"""

### Import Libraries
from dotenv import load_dotenv
import os
import gc
import psutil

load_dotenv('.env')
root = os.getenv('rootdir')
csxcad_location = os.getenv('csxcad_location')
os.add_dll_directory(root)

base_path = os.path.abspath(f'runs')

from datetime import datetime as dt
from random import randint
import numpy as np

import hashlib
import atexit
import shutil

from pylab import *
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *


def generate_mesh_lines(startpoint, center, stoppoint, axis, min_size, max_size, length, meshlines):
    """
    Generates mesh lines for a given axis, starting from a start point to a stop point
    through a center point. Mesh lines expand geometrically from the start to the center
    and are mirrored from the center to the stop point.

    Parameters:
    - startpoint (float): The starting point of the mesh.
    - center (float): The center point of the mesh.
    - stoppoint (float): The stopping point of the mesh.
    - axis (str or int): The axis along which to generate the mesh lines.
    - min_size (float): The minimum size of mesh increments.
    - max_size (float): The maximum size of mesh increments.
    - length (float): The length of the interval from start to stop.
    - meshlines (dict): A dictionary containing lists of mesh lines for each axis.

    Returns:
    - None: The function modifies the 'meshlines' dictionary in place.
    """

    if startpoint != stoppoint:
        # Normalize the direction vector
        direction = int(length / abs(length))

        # Add the first point offset by the direction
        meshlines[axis].append(startpoint + direction * min_size * 2 / 3)
        # Initialize the current point and size
        meshlines[axis].append(startpoint - direction * min_size * 1 / 3)
        # Generate mirrored points from center to stop
        inverted = (np.array(meshlines[axis][::-1]) - startpoint) * -1 + stoppoint
        for point in inverted:
            meshlines[axis].append(point)

        # Ensure the distance between two middle points is not larger than max_size
        mid_index = len(meshlines[axis]) // 2
        # if abs(meshlines[axis][mid_index] - meshlines[axis][mid_index - 1]) > max_size:
        #    meshlines[axis].insert(mid_index, center)
        if abs(meshlines[axis][mid_index] - meshlines[axis][mid_index - 1]) > length / 2:
            meshlines[axis].insert(mid_index, center)


def extend_line(start, stop, min_size=0.2, max_size=4., min_cells=3, max_cells=10):
    # Calculate the total range for each dimension
    total_range = [stop[i] - start[i] for i in range(3)]

    lines = {"x": [], "y": [], "z": []}
    meshlines = [[], [], []]
    axis = 0
    min_size_bcp = 0.1
    for axis in range(3):
        startpoint = start[axis]
        stoppoint = stop[axis]
        center = (startpoint + stoppoint) / 2
        length = startpoint - stoppoint
        retries = 0
        generate_mesh_lines(startpoint, center, stoppoint, axis, min_size, max_size, length, meshlines)

    return meshlines


def initialize_simulation(f0, fc, max_timesteps):
    """
    Initialize the FDTD and CSX structures.

    Returns:
    - FDTD object
    - CSX object
    """
    # Initialize openEMS
    FDTD = openEMS(NrTS=max_timesteps)  # , EndCriteria=1e-5)
    FDTD.SetGaussExcite(f0, fc)
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

    # Initialize CSXCAD geometry & mesh
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    return FDTD, CSX


def setup_mesh(mesh, SimBox, unit):
    """
    Initialize the mesh with the air-box dimensions.
    """
    mesh.SetDeltaUnit(unit)
    # Initialize the mesh with the "air-box" dimensions
    mesh.AddLine('x', [-SimBox[0] / 2, SimBox[0] / 2])
    mesh.AddLine('y', [-SimBox[1] / 2, SimBox[1] / 2])
    mesh.AddLine('z', [-SimBox[2] / 2, SimBox[2] / 2])


def create_substrate(CSX, parameters, mesh, FDTD):
    """
    Create the substrate in the CSX structure using parameters from the dictionary.
    """
    substrate_width = parameters['substrate_width']
    substrate_length = parameters['substrate_length']
    substrate_thickness = parameters['substrate_thickness']
    substrate_epsR = parameters['substrate_epsR']
    substrate_kappa = parameters.get('substrate_kappa', 1e-3 * 2 * np.pi * parameters['center_freq'] * EPS0 * substrate_epsR)
    substrate_cells = parameters['substrate_cells']
    unit = parameters['unit']
    ant_h = parameters['ant_h']
    ant_e = parameters['ant_e']
    
    # Create substrate material
    substrate = CSX.AddMaterial('substrate', epsilon=substrate_epsR, kappa=substrate_kappa)
    start = [-substrate_width / 2, -substrate_length / 2, 0]
    stop = [substrate_width / 2, substrate_length / 2 + ant_h + ant_e, substrate_thickness]
    substrate.AddBox(start=start, stop=stop, priority=1)

    # Add edges and additional mesh lines
    FDTD.AddEdges2Grid(dirs='xy', properties=substrate)
    mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells + 1))


def create_ground_plane(CSX, parameters, mesh, FDTD):
    """
    Create the ground plane in the CSX structure using parameters from the dictionary.
    """
    substrate_width = parameters['substrate_width']
    substrate_length = parameters['substrate_length']
    substrate_thickness = parameters['substrate_thickness']
    ant_e = parameters['ant_e']
    
    gndplane_position = parameters['gndplane_position']

    overlap=parameters['overlap']
    ant_h = parameters['ant_h']  # antenna height (y-direction)
    ant_l = parameters['ant_l']  # antenna length (x-direction)
    antenna_grid = parameters['antenna_grid']  # 2D grid of the antenna pattern
    num_cells_x = antenna_grid.shape[1]
    num_cells_y = antenna_grid.shape[0]
    # Define the physical size of each cell
    cell_size_x = ant_l / num_cells_x
    cell_size_y = ant_h / num_cells_y

    # Create ground plane material
    gnd = CSX.AddMetal('groundplane')
    start = [-substrate_width / 2 + ant_e, -substrate_length / 2 + ant_e, substrate_thickness + gndplane_position]
    stop = [substrate_width / 2 - ant_e, substrate_length / 2 - ant_e, substrate_thickness + gndplane_position]
    #start = [-substrate_width / 2 + ant_e-cell_size_x*overlap, -substrate_length / 2 + ant_e-ant_e-cell_size_y*overlap, substrate_thickness + gndplane_position]
    #stop = [substrate_width / 2 - ant_e+cell_size_x*overlap, substrate_length / 2 - ant_e+cell_size_y*overlap, substrate_thickness + gndplane_position]
    gnd.AddBox(start=start, stop=stop, priority=10)

    # Add edges and optional mesh line
    FDTD.AddEdges2Grid(dirs='xy', properties=gnd)

    start = [-substrate_width / 2 + ant_e-cell_size_x*overlap, -substrate_length / 2 + ant_e-ant_e, substrate_thickness + gndplane_position]
    stop = [substrate_width / 2 - ant_e+cell_size_x*overlap, substrate_length / 2 - ant_e+cell_size_y*overlap, substrate_thickness + gndplane_position]
    gnd.AddBox(start=start, stop=stop, priority=10)


    # Optional: add mesh line along 'x' if needed
    # mesh.AddLine("x", start)

def find_contiguous_blocks(grid):
    """
    Find contiguous blocks of filled cells in the grid.
    Returns a list of rectangles defined by (x_start, y_start, x_end, y_end)
    """
    blocks = []
    visited = np.zeros_like(grid, dtype=bool)

    for y in range(grid.shape[0]):
        for x in range(grid.shape[1]):
            if grid[y, x] == 1 and not visited[y, x]:
                # Start a new block
                x_start, y_start = x, y
                x_end, y_end = x, y

                # Expand in x-direction
                while x_end + 1 < grid.shape[1] and grid[y, x_end + 1] == 1 and not visited[y, x_end + 1]:
                    x_end += 1

                # Expand in y-direction
                expand_in_y = True
                while expand_in_y and y_end + 1 < grid.shape[0]:
                    for xi in range(x_start, x_end + 1):
                        if grid[y_end + 1, xi] == 0 or visited[y_end + 1, xi]:
                            expand_in_y = False
                            break
                    if expand_in_y:
                        y_end += 1

                # Mark cells as visited
                for yi in range(y_start, y_end + 1):
                    for xi in range(x_start, x_end + 1):
                        visited[yi, xi] = True

                # Add the block to the list
                blocks.append((x_start, y_start, x_end, y_end))

    return blocks

def makearray(num_cells_x,num_cells_y,antenna_grid = None, makesafe=False):
    if antenna_grid is None:
        antenna_grid = np.zeros((num_cells_y, num_cells_x), dtype=int)
    import random
    #random.seed(78)
    for x in range(0, num_cells_x):
        for y in range(0, num_cells_y):
            if random.random() < 0.4:
                antenna_grid[y, x] = 1
    
    return antenna_grid

def create_ga(FDTD, CSX, mesh, parameters):
    """
    Create the GA-based antenna in the CSX structure.
    """
    # Extract parameters needed
    unit = parameters['unit']
    substrate_width = parameters['substrate_width']
    substrate_length = parameters['substrate_length']
    substrate_thickness = parameters['substrate_thickness']
    ant_h = parameters['ant_h']  # antenna height (y-direction)
    ant_l = parameters['ant_l']  # antenna length (x-direction)
    ant_fp = parameters['ant_fp']  # feedpoint position along x
    ant_e = parameters['ant_e']  # edge distance
    overlap = parameters['overlap']  # overlap distance
    gndplane_position = parameters['gndplane_position']  # depth of the ground plane in z
    antenna_grid = parameters['antenna_grid']  # 2D grid of the antenna pattern

    # Create IFA material
    ifa_material = CSX.AddMetal('ifa')

    # Define the top-left coordinate (origin for the antenna grid)
    tl = np.array([-substrate_width / 2 + ant_e, substrate_length / 2 +ant_h- ant_e,
                   substrate_thickness])  # translation vector
    
    
    num_cells_x = antenna_grid.shape[1]
    num_cells_y = antenna_grid.shape[0]
    # Define the physical size of each cell
    cell_size_x = ant_l / num_cells_x
    cell_size_y = ant_h / num_cells_y
    

    # Store cell sizes in parameters for use in add_feed
    parameters['cell_size_x'] = cell_size_x
    parameters['cell_size_y'] = cell_size_y

    # Define the feedpoint coordinates
    feed_cell_x = int(ant_fp / cell_size_x)
    feed_cell_y = num_cells_y - 1  # Assuming feed is at the bottom row
    
    
    for i in range(0,num_cells_x-1):
        antenna_grid[feed_cell_y, i] = 0  # never cover the feed position
        #antenna_grid[feed_cell_y, feed_cell_x-i] = 0  # never cover the feed position
    
    antenna_grid[-5:-1, feed_cell_x]=1
    
    feed_x = tl[0] + feed_cell_x * cell_size_x + cell_size_x / 2
    feed_y = tl[1] - feed_cell_y * cell_size_y - cell_size_y / 2
    feed_z = tl[2]

    # Store feedpoint coordinates in parameters for later use
    parameters['feed_point'] = np.array([feed_x, feed_y, feed_z])

    # Create a test grid (this should be replaced with your GA-generated grid)
    

    # Example pattern (modify this as needed)
    # Example pattern: a diagonal line

    # Find contiguous blocks in the grid
    blocks = find_contiguous_blocks(antenna_grid)

    # For each block, create a box in CSXCAD
    for block in blocks:
        x_start, y_start, x_end, y_end = block

        # Calculate the physical coordinates
        start_x = tl[0] + x_start * cell_size_x-cell_size_x*overlap
        stop_x = tl[0] + (x_end + 1) * cell_size_x+cell_size_x*overlap # +1 because cells are zero-indexed

        start_y = tl[1] - y_start * cell_size_y+cell_size_y*overlap
        stop_y = tl[1] - (y_end + 1) * cell_size_y-cell_size_y*overlap

        start = np.array([start_x, stop_y, tl[2]])  # Note that stop_y < start_y because y decreases
        stop = np.array([stop_x, start_y, tl[2]])   # z-coordinate is constant (substrate_thickness)

        # Add the box to the CSX structure
        ifa_material.AddBox(start, stop, priority=10)

    x_meshlines = [tl[0] + i * cell_size_x/2 for i in range(num_cells_x*2 + 1)]
    y_meshlines = [tl[1] - i * cell_size_y/2 for i in range(num_cells_y*2 + 1)]

    # Add meshlines along X and Y
    mesh.AddLine('x', x_meshlines)
    mesh.AddLine('y', y_meshlines)



def add_feed(FDTD, CSX, mesh, parameters):
    """
    Apply the excitation & resistor as a current source using the feed point from parameters.
    Create the feed as a box the same size as the metal pads, aligned to the grid.
    """
    # Extract necessary parameters
    feed_point = parameters['feed_point']
    gndplane_position = parameters['gndplane_position']
    substrate_thickness = parameters['substrate_thickness']
    substrate_width = parameters['substrate_width']
    overlap = parameters['overlap']
    ant_e = parameters['ant_e']
    feed_R = parameters['feed_R']
    cell_size_x = parameters['cell_size_x']
    cell_size_y = parameters['cell_size_y']
    ifa_material = CSX.AddMetal('ifa')  # Using the same metal material as the antenna
    print(feed_point)
    # Determine the feed box dimensions
    start_x = feed_point[0]-cell_size_x/2-cell_size_x*overlap
    stop_x = feed_point[0] + cell_size_x/2+cell_size_x*overlap
    start_y = feed_point[1]-cell_size_y/2-cell_size_y*overlap
    stop_y = feed_point[1]+cell_size_y/2+cell_size_y*overlap
    start_z = feed_point[2]

    if gndplane_position != 0:
        # Feed connects vertically (z-direction) to the ground plane at z = substrate_thickness + gndplane_position
        feed_direction = 'z'
        ground_plane_z = substrate_thickness + gndplane_position
        stop_z = ground_plane_z
        start_coord = np.array([start_x, start_y, start_z])
        stop_coord = np.array([stop_x, stop_y, stop_z])
    else:
        # Feed connects horizontally (x-direction) to the ground plane at x = -substrate_width / 2 + ant_e
        feed_direction = 'y'
        ground_plane_x = -substrate_width / 2 + ant_e
        #
        start_coord = np.array([start_x, start_y, start_z])
        stop_coord = np.array([stop_x, stop_y, start_z])  # z remains the same

    # Add the feed box to the CSX structure
    #ifa_material.AddBox(start_coord, stop_coord, priority=20)  # Higher priority for meshing

    # Add meshlines along the feed box edges
    #mesh.AddLine('x', [start_coord[0], stop_coord[0]])
    #mesh.AddLine('y', [start_coord[1], stop_coord[1]])
    mesh.AddLine('z', [start_coord[2], stop_coord[2]])

    # Now, add the lumped port between the feed box and the ground plane
    if gndplane_position != 0:
        # Lumped port is at the bottom face of the feed box in z-direction
        port_start = np.array([start_x, start_y, start_z])
        port_stop = np.array([stop_x, stop_y, stop_z])
    else:
        # Lumped port is at the side face of the feed box in y-direction
        port_start = np.array([start_x, start_y, start_z])
        port_stop = np.array([stop_x, stop_y, start_z])

    # Add the lumped port to the FDTD simulation
    port = FDTD.AddLumpedPort(1, feed_R, port_start, port_stop, feed_direction, 1.0, priority=5)

    return port

def prepare_simulation_directory(Sim_Path, Sim_CSX, CSX, showCad, csxcad_location):
    """
    Prepare the simulation directory and write the CSX file.
    """
    CSX_file = os.path.join(Sim_Path, Sim_CSX)
    CSX.Write2XML(CSX_file)

    # Show the structure
    if showCad:
        print("showing cad")
        print(f"csxcad location{csxcad_location}")
        print(f"csxfile location {CSX_file}")
        os.system(csxcad_location + ' "{}"'.format(CSX_file))


def run_simulation(FDTD, Sim_Path, sim_file, temp_file,parameters):
    """
    Run the simulation.
    """
    if not os.path.exists(sim_file):
        try:
            with open(temp_file, 'w') as f:
                f.write('Running')
            FDTD.Run(Sim_Path, verbose=0, cleanup=False,numThreads=parameters.get('numThreads',0))
            os.remove(temp_file)
            with open(sim_file, 'w') as f:
                f.write('Completed')
        except Exception as e:
            print("An error occurred during simulation:", e)
            import traceback
            traceback.print_exc()


def post_process_results(Sim_Path, port, freq, delete_simulation_files, plot, center_freq, nf2ff, parameters):
    """
    Post-process the simulation results and plot.
    """
    port.CalcPort(Sim_Path, freq)

    Zin = port.uf_tot / port.if_tot
    s11 = port.uf_ref / port.uf_inc
    s11_dB = 20.0 * np.log10(np.abs(s11))
    P_in = np.real(0.5 * port.uf_tot * np.conj(port.if_tot))

    thetaRange = np.arange(0, 182, 2)
    phiRange = np.arange(0, 362, 2) - 180
    idx = np.where((s11_dB < -10) & (s11_dB == np.min(s11_dB)))[0]

    if len(idx) != 1:
        print('No resonance frequency found for far-field calculation')
    else:
        f_res_ind = np.argmin(np.abs(s11))
        f_res = freq[f_res_ind]

        theta = np.arange(-180.0, 180.0, 2.0)
        # Create the bottom-right polar subplot (axs[1, 1]) for the xy-plane
        phi = theta
        nf2ff_res_theta90 = nf2ff.CalcNF2FF(Sim_Path, f_res, 90, phi, center=np.array([0, 0, 0]), read_cached=True,
                                            outfile='nf2ff_xy.h5')

        print('Radiated power: Prad = {:.2e} Watt'.format(nf2ff_res_theta90.Prad[0]))
        print('Directivity:    Dmax = {:.1f} ({:.1f} dBi)'.format(nf2ff_res_theta90.Dmax[0],
                                                                  10 * np.log10(nf2ff_res_theta90.Dmax[0])))
        print('Efficiency:   nu_rad = {:.1f} %'.format(100 * nf2ff_res_theta90.Prad[0] / np.real(P_in[idx[0]])))

        print(f"Resonance frequency: {f_res / 1e9} GHz")
        # s11 at closste to center frequency
        print(f"S11 at resonance frequency: {s11_dB[f_res_ind]} dB")

    print(f"S11 at center frequency{s11_dB[freq == center_freq]} dB")

    if plot:
        # Create a figure with subplots, 2 rows and 2 columns
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        # Plot feed point impedance in the top-left subplot (axs[0, 0])
        axs[0, 0].plot(freq / 1e6, np.real(Zin), 'k-', linewidth=2, label='Re{Zin}')
        axs[0, 0].plot(freq / 1e6, np.imag(Zin), 'r--', linewidth=2, label='Im{Zin}')
        # Add a horizontal line at 50 Ohms
        axs[0, 0].axhline(parameters["feed_R"], color='blue', linestyle='--', linewidth=1,
                          label=f'{parameters["feed_R"]} Ohm')
        axs[0, 0].axvline(parameters["min_freq"] / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{parameters["min_freq"] / 1e6} MHz')
        axs[0, 0].axvline(parameters["center_freq"] / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{parameters["center_freq"] / 1e6} MHz')
        axs[0, 0].axvline(parameters["max_freq"] / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{parameters["max_freq"] / 1e6} MHz')

        # Customize the grid, title, labels, and legend
        axs[0, 0].grid()
        axs[0, 0].set_title('Feed point impedance')
        axs[0, 0].set_xlabel('Frequency f / MHz')
        axs[0, 0].set_ylabel('Impedance Zin / Ohm')
        axs[0, 0].legend()

        # Plot reflection coefficient S11 in the top-right subplot (axs[0, 1])
        axs[0, 1].plot(freq / 1e6, 20 * np.log10(np.abs(s11)), 'k-', linewidth=2, label='S11(db)')

        cutoffDB = -6
        axs[0, 1].axhline(cutoffDB, color='blue', linestyle='--', linewidth=1, label=f'{cutoffDB} dB Cutoff')
        axs[0, 1].axvline(parameters["min_freq"] / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{parameters["min_freq"] / 1e6} MHz')
        axs[0, 1].axvline(center_freq / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{center_freq / 1e6} MHz')
        axs[0, 1].axvline(parameters["max_freq"] / 1e6, color='green', linestyle='--', linewidth=1,
                          label=f'{parameters["max_freq"] / 1e6} MHz')

        axs[0, 1].grid()
        axs[0, 1].set_title('Reflection coefficient S11')
        axs[0, 1].set_xlabel('Frequency f / MHz')
        axs[0, 1].set_ylabel('Reflection coefficient |S11| (dB)')
        axs[0, 1].legend()
        # For polar subplots, calculate NF2FF if the resonance frequency is found

        if len(idx) != 1:
            print('No resonance frequency found for far-field calculation')
        else:
            # Find resonance frequency from s11
            print("Calculate NF2FF")
            nf2ff_res_phi0 = nf2ff.CalcNF2FF(Sim_Path, f_res, theta, 0, verbose=0, outfile='3D_Pattern.h5')

            # Create the bottom-left polar subplot (axs[1, 0]) for the xz-plane
            ax_polar1 = plt.subplot(223, polar=True)
            E_norm = 20.0 * np.log10(nf2ff_res_phi0.E_norm / np.max(nf2ff_res_phi0.E_norm)) + nf2ff_res_phi0.Dmax
            ax_polar1.plot(np.deg2rad(theta), 10 ** (np.squeeze(E_norm) / 20), linewidth=2, label='xz-plane')
            ax_polar1.grid(True)
            ax_polar1.set_xlabel('theta (deg)')
            ax_polar1.set_theta_zero_location('N')
            ax_polar1.set_theta_direction(-1)
            ax_polar1.legend(loc=3)

            ax_polar2 = plt.subplot(224, polar=True)
            E_norm = 20.0 * np.log10(
                nf2ff_res_theta90.E_norm / np.max(nf2ff_res_theta90.E_norm)) + nf2ff_res_theta90.Dmax
            ax_polar2.plot(np.deg2rad(phi), 10 ** (np.squeeze(E_norm) / 20), linewidth=2, label='xy-plane')
            ax_polar2.grid(True)
            ax_polar2.set_xlabel('phi (deg)')
            plt.suptitle('IFA Antenna Pattern\nFrequency: {:.2f} GHz'.format(f_res / 1e9), fontsize=14)
            ax_polar2.legend(loc=3)

        # Prepare the parameters text with table-like formatting
        parameters_text = (
            f"{'IFA Parameters:':<30}\n"
            f"{'ant_h:':<20}{parameters['ant_h']:>10.3f} mm "
            f"{'ant_l:':<20}{parameters['ant_l']:>10.3f} mm\n"
            f"{'ant_fp:':<20}{parameters['ant_fp']:>10.3f} mm\n"
            f"{'ant_e:':<20}{parameters['ant_e']:>10.3f} mm "
            f"{'ant_h + ant_e:':<20}{parameters['ant_h'] + parameters['ant_e']:>10.3f} mm\n\n"
            f"{'Substrate Properties:':<30}\n"
            f"{'substrate_epsR:':<20}{parameters.get('substrate_epsR', 'N/A'):>10} "
            f"{'substrate_kappa:':<20}{parameters.get('substrate_kappa', 'N/A'):>10} S/m\n\n"
            f"{'Feeding Setup:':<30}\n"
            f"{'feed_R:':<20}{parameters['feed_R']:>10} Ω"
        )
        # Add the parameters text to the figure at the bottom
        fig.text(0.5, 0.01, parameters_text, ha='center', va='bottom', fontsize=12, family='monospace')

        # Adjust layout to prevent overlap with the text
        plt.tight_layout(rect=[0, 0.2, 1, 0.95])  # Adjust rect to leave space at bottom and top

        # Show all plots
        plt.figure()
        plt.show()

    return freq, s11_dB, Zin, P_in

def mesh_divide(mesh, axes=['x', 'y', 'z'], num_parts=2):
    for axis in axes:
        mesh_lines = mesh.GetLines(axis)
        new_lines = []

        for i in range(len(mesh_lines) - 1):
            start_line = mesh_lines[i]
            end_line = mesh_lines[i + 1]

            # Add the starting line
            new_lines.append(start_line)

            # Calculate and add 'divisor' lines between start_line and end_line
            for j in range(1, num_parts):
                interpolated_line = start_line + (end_line - start_line) * (j / (num_parts))
                new_lines.append(interpolated_line)

        # Add the last line
        new_lines.append(mesh_lines[-1])

        # Set the new lines for the axis
        mesh.SetLines(axis, new_lines)
        
def ga_simulation(    parameters = {
        'Sim_CSX' : 'IFA.xml',
        'unit': 1e-3,
        'substrate_width': 21,
        'substrate_length': 20,
        'substrate_thickness': 1.5,
        'substrate_epsR': 4.5,
        'gndplane_position': 0,
        'substrate_cells': 4,
        'ant_h': 14,
        'ant_l': 20,
        'ant_fp': 5,
        'ant_e': 0.5,
        'feed_R': 50,
        'min_freq': 2.4e9,
        'center_freq': 2.45e9,
        'max_freq': 2.5e9,
        'overlap': 0.1,
        'fc': 1.0e9,
        'max_timesteps': 1000,
        'override_min_global_grid': None,
        'plot': True,  # Set to True to plot results
        'showCad': True,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': makearray(20, 20),
        'numthreads': 4,
        #'lambdamultiplier': 2,
    }):
    #############################################################################
    #                substrate_width
    #  _______________________________________________    __ substrate_thickness
    # | A        X                                    |\  __
    # | |        XXXX      X   X                      | |
    # | |ant_h      XXXX  XX                          | |
    # | |         X    XXXX   X                       | | ______
    # | |                   xxxx                      | ||mant_edgedistance (minimum value)
    # |_V______________________x______________________| ||______
    # | <-ant_e->|            |x                      | |
    # |                       |                       | |
    # |                       |                       | |
    # |                       |                       | | substrate_length
    # |<- ant_fp----------- ->|       x=metal         | |
    # |                                               | |
    # |_______________________________________________| |
    #  \_______________________________________________\|
    #
    # Note: It's not checked whether your settings make sense, so check
    #       graphical output carefully.
    #
    ##############################################################################
    
    # Create IFA


    params_tuple = tuple(parameters.items())

    params_str = '_'.join(map(str, params_tuple))
    # Generate a SHA256 hash of the parameters string
    params_hash = hashlib.sha256(params_str.encode('utf-8')).hexdigest()
    # Use the first 8 characters of the hash for brevity
    hash_prefix = params_hash[:12]
    # Create the simulation path using the hash
    Sim_Path = os.path.join(base_path, f'tmp_IFA_{hash_prefix}')

    substrate_width=parameters['substrate_width']
    substrate_length=parameters['substrate_length']
    ant_h=parameters['ant_h']
    max_timesteps=parameters['max_timesteps']
    unit=parameters['unit']
    xmultiplier = parameters.get('xmultiplier', 2)
    ymultiplier = parameters.get('ymultiplier', 2)
    zmultiplier = parameters.get('zmultiplier', 2)
    # Simulation box size
    
    
    if parameters.get('f0', None) is None:
        f0 = parameters["center_freq"]
    else:
        f0 = parameters["f0"]
    
    fc = parameters["fc"]
    
    if parameters.get('lambdamultiplier', None) is not None:
        quarter_wavelength = C0 / (f0 - fc) / unit / 4
        lambdamultiplier = parameters['lambdamultiplier']
        SimBox = np.array([lambdamultiplier*quarter_wavelength + substrate_width, lambdamultiplier*quarter_wavelength + (substrate_length + ant_h),
                        lambdamultiplier*quarter_wavelength + (substrate_length + ant_h)])
    else:
        SimBox = np.array([substrate_width * xmultiplier, (substrate_length + 2 * ant_h) * ymultiplier, (substrate_length + 2 * ant_h) * zmultiplier])
    
    
    
    # Initialize simulation
    FDTD, CSX = initialize_simulation(f0, fc, max_timesteps)
    mesh = CSX.GetGrid()
    setup_mesh(mesh, SimBox, unit)
    
    
    # Create substrate
    create_substrate(CSX, parameters, mesh, FDTD)

    # Create ground plane
    create_ground_plane(CSX, parameters, mesh, FDTD)

    

    create_ga(FDTD, CSX, mesh, parameters)

    port = add_feed(FDTD, CSX, mesh, parameters)

    # Finalize the mesh

    Sim_CSX = parameters["Sim_CSX"]
    
    
    mesh_res = C0 / (f0 + fc) / unit / 20
    if parameters.get('override_min_global_grid', None) is not None:
        mesh_res = parameters['override_min_global_grid']
        
    mesh.SmoothMeshLines('all', mesh_res, 1.5)
    if parameters.get('mesh_divide', None) is not None:
        mesh_divide(mesh, axes=['x', 'y', 'z'], num_parts=parameters.get('mesh_divide'))
    nf2ff = FDTD.CreateNF2FFBox()
    #finalize_mesh(mesh, min_size, f0, fc, unit, override_min_global_grid,
    #              {"x": [], "y": [], "z": [0, substrate_thickness + gndplane_position, substrate_thickness]})

    temp_file = os.path.join(Sim_Path, 'incomplete_run.flag')

    if os.path.exists(Sim_Path):
        if os.path.exists(temp_file) or parameters.get('allways_rerun', True):
            import shutil
            print(f"Cleaning up incomplete run: {Sim_Path}")
            shutil.rmtree(Sim_Path)

    if not os.path.exists(Sim_Path):
        print(f"Creating directory {Sim_Path}")
        try:
            os.mkdir(Sim_Path)
        except OSError:
            print(f"Creation of the directory {Sim_Path} failed")
            #make the missing folder
            missingpath = os.path.dirname(Sim_Path)
            os.makedirs(missingpath)
            os.mkdir(Sim_Path)
            

    CSX_file = os.path.join(Sim_Path, Sim_CSX)
    CSX.Write2XML(CSX_file)

    showCad = parameters.get('showCad', False)
    # Show the structure
    if showCad:
        print("showing cad")
        print(f"csxcad location{csxcad_location}")
        print(f"csxfile location {CSX_file}")
        os.system(csxcad_location + ' "{}"'.format(CSX_file))

    sim_file = os.path.join(Sim_Path, 'complete_run.flag')
    
    post_proc_only=parameters.get('post_proc_only', False)
    if not post_proc_only and not os.path.exists(sim_file):
        run_simulation(FDTD, Sim_Path, sim_file, temp_file,parameters)
    delete_simulation_files = parameters.get('delete_simulation_files', True)
    plot = parameters.get('plot', False)
    if os.path.exists(sim_file):
        # Post-processing & plotting
        freq = np.linspace(max(0, f0 - fc), f0 + fc, 501)
        freq, s11_dB, Zin, P_in = post_process_results(Sim_Path, port, freq, delete_simulation_files, plot, f0,
                                                    nf2ff, parameters)
        return freq, s11_dB, Zin, P_in, hash_prefix

if __name__ == "__main__":
    ga_simulation()
