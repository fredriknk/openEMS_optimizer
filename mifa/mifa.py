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
load_dotenv('.env')
root = os.getenv('rootdir')
csxcad_location= os.getenv('csxcad_location')
os.add_dll_directory(root)

base_path=os.path.abspath(f'runs')

from datetime import datetime as dt
from random import randint
import numpy as np
import tempfile
import hashlib

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
        meshlines[axis].append( startpoint - direction * min_size * 1 / 3)
        # Generate mirrored points from center to stop
        inverted = (np.array(meshlines[axis][::-1]) - startpoint) * -1 + stoppoint
        for point in inverted:
            meshlines[axis].append(point)

        # Ensure the distance between two middle points is not larger than max_size
        mid_index = len(meshlines[axis]) // 2
        #if abs(meshlines[axis][mid_index] - meshlines[axis][mid_index - 1]) > max_size:
        #    meshlines[axis].insert(mid_index, center)
        if abs(meshlines[axis][mid_index] - meshlines[axis][mid_index - 1]) > length/2:
            meshlines[axis].insert(mid_index, center)

def extend_line(start, stop,min_size=0.2,max_size =4.,min_cells=3,max_cells=10):
    # Calculate the total range for each dimension
    total_range = [stop[i] - start[i] for i in range(3)]

    lines = {"x":[],"y":[],"z":[]}
    meshlines=[[],[],[]]
    axis = 0
    min_size_bcp = 0.1
    for axis in range(3):
        startpoint= start[axis]
        stoppoint = stop[axis]
        center = (startpoint+stoppoint)/2
        length = startpoint-stoppoint
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
    FDTD = openEMS(NrTS=max_timesteps)#, EndCriteria=1e-5)
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
    mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
    mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
    mesh.AddLine('z', [-SimBox[2]/2, SimBox[2]/2])

def create_substrate(CSX, substrate_width, substrate_length, substrate_thickness, ifa_h, ifa_e, substrate_epsR, substrate_kappa, substrate_cells, unit, mesh):
    """
    Create the substrate in the CSX structure.
    """
    # Create substrate
    substrate = CSX.AddMaterial('substrate', epsilon=substrate_epsR, kappa=substrate_kappa)

    start = [-substrate_width/2, -substrate_length/2, 0]
    stop = [substrate_width/2, substrate_length/2+ifa_h+ifa_e, substrate_thickness]
    substrate.AddBox(start=start, stop=stop, priority=1)

    # Add extra cells to discretize the substrate thickness
    mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells+1))

def create_ground_plane(CSX, substrate_width, substrate_length, substrate_thickness, ifa_e, gndplane_position, mesh, min_size, max_size):
    """
    Create the ground plane in the CSX structure.
    """
    gnd = CSX.AddMetal('groundplane')
    start = [-substrate_width/2+ifa_e, -substrate_length/2+ifa_e, substrate_thickness+gndplane_position]
    stop = [substrate_width/2-ifa_e, substrate_length/2 - ifa_e, substrate_thickness+gndplane_position]
    gnd.AddBox(start=start, stop=stop, priority=10)

    meshlines = extend_line(start, stop, min_size, max_size)
    mesh.AddLine('x', meshlines[0])
    mesh.AddLine('y', meshlines[1])

def create_ifa(CSX, mesh, parameters):
    """
    Create the IFA in the CSX structure.
    """
    # Extract parameters needed
    unit = parameters['unit']
    substrate_width = parameters['substrate_width']
    substrate_length = parameters['substrate_length']
    substrate_thickness = parameters['substrate_thickness']
    ifa_h = parameters['ifa_h']
    ifa_l = parameters['ifa_l']
    ifa_w1 = parameters['ifa_w1']
    ifa_w2 = parameters['ifa_w2']
    ifa_wf = parameters['ifa_wf']
    ifa_fp = parameters['ifa_fp']
    ifa_e = parameters['ifa_e']
    mifa_meander = parameters['mifa_meander']
    mifa_tipdistance = parameters['mifa_tipdistance']
    mifa_meander_edge_distance = parameters['mifa_meander_edge_distance']
    gndplane_position = parameters['gndplane_position']
    min_size = parameters['min_size']
    max_size = parameters['max_size']

    # Create IFA
    ifa = CSX.AddMetal('ifa')
    tl = np.array([-substrate_width/2+ifa_e+ifa_fp, substrate_length / 2 - ifa_e, substrate_thickness])  # translation vector
    via_diameter=0.3
    if gndplane_position != 0:
        print(f"offset groundplane, add vias")
        feedpoint = -via_diameter #set feedpoint to groundplane edge
        short_circuit_stub = via_diameter
    else:
        feedpoint = ifa_wf
        short_circuit_stub = 0
    # Feed element
    start = np.array([0, feedpoint, 0]) + tl
    stop = start + np.array([ifa_wf, ifa_h - feedpoint, 0])
    ifa.AddBox(start=start, stop=stop, priority=10)
    meshlines = extend_line(start, stop, min_size, max_size)
    #mesh.AddLine('x',meshlines[0] )
    #mesh.AddLine('y', meshlines[1])

    # Short circuit stub
    start = np.array([-ifa_fp, -short_circuit_stub, 0]) + tl
    stop = start + np.array([ifa_w1, ifa_h+short_circuit_stub, 0])
    ifa.AddBox(start=start, stop=stop, priority=10)

    meshlines = extend_line(start, stop,min_size,max_size)
    #mesh.AddLine('x',meshlines[0] )
    #mesh.AddLine('y', meshlines[1])

    if gndplane_position != 0:
        #Add Via to groundplane
        start =start+np.array([+ifa_w1/2,+via_diameter/2,0])
        stop = start+np.array([0,0,gndplane_position])
        ifa.AddCylinder(start,stop,via_diameter/2,priority=10)

        start = stop+np.array([-via_diameter/2,+via_diameter/2,0])
        stop = start + np.array([via_diameter, via_diameter,-gndplane_position])
        meshlines = extend_line(start, stop,min_size,max_size)
        mesh.AddLine('x', meshlines[0])
        mesh.AddLine('y', meshlines[1])
        mesh.AddLine('z', meshlines[2])

    # Radiating element
    if ifa_l < substrate_width-2*ifa_e:
        print('Ifa: ifa_l smaller than substrate width just do normal radiating element')
        start = np.array([-ifa_fp , ifa_h, 0]) + tl
        stop = start + np.array([ifa_l, -ifa_w2, 0])
        ifa.AddBox(start=start, stop=stop, priority=10)
    else:
                #check if we can add a tip element
        length_diff = (ifa_l - (substrate_width-2*ifa_e))
        print('ifa_l larger than substrate width starting meandering')

        radiating_element_start = np.array([-ifa_fp , ifa_h, 0]) + tl
        radiating_element_stop = radiating_element_start + np.array([substrate_width-2*ifa_e, -ifa_w2, 0])


        #distance to gndedge = ifa_h
        max_length_mifa = ifa_h-mifa_meander_edge_distance
        max_edgelength_tip = ifa_h-mifa_tipdistance

        if length_diff < max_edgelength_tip:
            print(f"only the tip {length_diff} < {max_edgelength_tip}")
            start= radiating_element_stop + np.array([0,+ifa_w2,0])
            stop = start + np.array([-ifa_w2, -length_diff-ifa_w2, 0])

            ifa.AddBox(start=start, stop=stop, priority=10)
            meshlines = extend_line(start, stop,min_size,max_size)
            mesh.AddLine('x',meshlines[0] )
            mesh.AddLine('y', meshlines[1])

            length_diff -= max_edgelength_tip
            ifa.AddBox(start=radiating_element_start, stop=radiating_element_stop, priority=10)
            meshlines = extend_line(radiating_element_start, radiating_element_stop,min_size,max_size)
            mesh.AddLine('x',meshlines[0] )
            mesh.AddLine('y', meshlines[1])
            length_diff -= max_edgelength_tip

        else:
            #maximize tip element
            start= radiating_element_stop + np.array([0,+ifa_w2,0])
            stop = start + np.array([-ifa_w2, -max_edgelength_tip-ifa_w2, 0])
            ifa.AddBox(start=start, stop=stop, priority=10)
            meshlines = extend_line(start, stop,min_size,max_size)
            #mesh.AddLine('x',meshlines[0] )
            #mesh.AddLine('y', meshlines[1])
            length_diff -= max_edgelength_tip
            print("Adding tip element{max_edgelength_tip}")

            #add meanders for remainder of distance
            ldfiff_ratio = length_diff /(max_length_mifa*2)
            print(f"ldfiff_ratio: {ldfiff_ratio}")
            if ldfiff_ratio > 0:
                stop = start+np.array([-ifa_w2,0,0])
                while ldfiff_ratio > 0:
                    print(f"current ldfiff_ratio: {ldfiff_ratio}")
                    if ldfiff_ratio < 1:
                        current_meander = ldfiff_ratio
                        ldfiff_ratio = 0
                    else:
                        current_meander = 1
                        ldfiff_ratio -= 1
                    print(f"Adding meanders ratio: {current_meander}")
                    if current_meander < 0.005:
                        print(f"current_meander too small: {current_meander}")
                        break
                    #first topline
                    start = stop+np.array([ifa_w2,0,0])
                    stop = start + np.array([-mifa_meander-ifa_w2,-ifa_w2,0])
                    ifa.AddBox(start=start, stop=stop, priority=10)
                    meshlines = extend_line(start, stop, min_size, max_size)
                    mesh.AddLine('x', meshlines[0])
                    mesh.AddLine('y', meshlines[1])

                    #meander down
                    start = stop + np.array([0,ifa_w2,0])
                    stop = start + np.array([ifa_w2,-current_meander*max_length_mifa-ifa_w2,0])
                    ifa.AddBox(start=start, stop=stop, priority=10)
                    meshlines = extend_line(start, stop, min_size, max_size)
                    mesh.AddLine('x', meshlines[0])
                    mesh.AddLine('y', meshlines[1])

                    #lower meander sideways
                    start = stop+ np.array([0,0,0])
                    stop = start + np.array([-mifa_meander-ifa_w2,ifa_w2,0])
                    ifa.AddBox(start=start, stop=stop, priority=10)
                    meshlines = extend_line(start, stop, min_size, max_size)
                    mesh.AddLine('x', meshlines[0])
                    mesh.AddLine('y', meshlines[1])

                    #meander up
                    start = stop + np.array([+ifa_w2,-ifa_w2,0])
                    stop = start + np.array([-ifa_w2,+current_meander*max_length_mifa+ifa_w2,0])
                    ifa.AddBox(start=start, stop=stop, priority=10)
                    meshlines = extend_line(start, stop, min_size, max_size)
                    mesh.AddLine('x', meshlines[0])
                    mesh.AddLine('y', meshlines[1])

                #connect to shorting stub
                start = radiating_element_start
                stop = stop+ np.array([ifa_w2,-ifa_w2,0])
                ifa.AddBox(start=start, stop=stop, priority=10)

    meshlines = extend_line(start, stop,min_size,max_size)
    mesh.AddLine('x',meshlines[0] )
    mesh.AddLine('y', meshlines[1])

    if gndplane_position != 0:
        feed_direction = 'z'
    else:
        feed_direction = 'x'
    #Apply the excitation & resistor as a current source
    start = np.array([0, 0, 0]) + tl
    stop = start + np.array([ifa_wf, feedpoint, gndplane_position])

    # (Continue with the rest of the code for creating IFA)

def add_feed(FDTD, CSX, mesh, tl, ifa_wf, feedpoint, gndplane_position, feed_R, min_size, max_size):
    """
    Apply the excitation & resistor as a current source.
    """
    if gndplane_position != 0:
        feed_direction = 'z'
    else:
        feed_direction = 'x'

    start = np.array([0, 0, 0]) + tl
    stop = start + np.array([ifa_wf, feedpoint, gndplane_position])

    port = FDTD.AddLumpedPort(1, feed_R, start, stop, feed_direction, 1.0, priority=5)
    meshlines = extend_line(start, stop, min_size, max_size)
    if gndplane_position != 0:
        mesh.AddLine('x',meshlines[0] )
        mesh.AddLine('y', meshlines[1])
        mesh.AddLine('z', stop[2])
    else:
        mesh.AddLine('x',meshlines[0] )
        mesh.AddLine('y', meshlines[1])

    return port

def finalize_mesh(mesh, min_size, f0, fc, unit, override_min_global_grid, immunity_array):
    """
    Finalize the mesh.
    """
    # Generate a smooth mesh with max. cell size: lambda_min / 20
    mesh_res = C0 / (f0 + fc) / unit / 20
    if override_min_global_grid is not None:
        mesh_res = override_min_global_grid

    TOL = 1e-9  # Tolerance for floating point comparison

    # Function to check if a line is immune
    def is_immune(line, immune_lines):
        for immune_line in immune_lines:
            if abs(line - immune_line) < TOL:
                return True
        return False

    for axis in ['x', 'y', 'z']:
        lines = mesh.GetLines(axis)
        print(f"meshlines {axis} before filtering: {lines}")

        # Sort lines to process them in order
        lines.sort()

        # Get immune lines for this axis
        immune_lines = immunity_array.get(axis, [])

        # Initialize filtered lines with the first line
        filtered_lines = [lines[0]]

        # Iterate over the remaining lines
        for l in lines[1:]:
            last_line = filtered_lines[-1]
            distance = abs(l - last_line)

            immune_l = is_immune(l, immune_lines)
            immune_last_line = is_immune(last_line, immune_lines)

            if distance >= min_size / 4:
                # Lines are sufficiently apart, keep the current line
                filtered_lines.append(l)
            else:
                # Lines are too close
                if immune_l and immune_last_line:
                    # Both lines are immune, keep both
                    filtered_lines.append(l)
                elif immune_l:
                    # Current line is immune, remove last_line if not immune
                    if not immune_last_line:
                        print(f"Removing line {last_line} on {axis}-axis; too close to immune line {l} (distance {distance})")
                        filtered_lines[-1] = l  # Replace last_line with l
                    else:
                        # Both lines are immune, keep both
                        filtered_lines.append(l)
                elif immune_last_line:
                    # Last line is immune, remove l
                    print(f"Removing line {l} on {axis}-axis; too close to immune line {last_line} (distance {distance})")
                    # Do not add l to filtered_lines
                    pass
                else:
                    # Neither line is immune, remove current line l
                    print(f"Removing line {l} on {axis}-axis; too close to {last_line} (distance {distance})")
                    # Do not add l to filtered_lines
                    pass

        mesh.SetLines(axis, filtered_lines)
        print(f"meshlines {axis} after filtering: {mesh.GetLines(axis)}")

    mesh.SmoothMeshLines('all', mesh_res, 1.5)
    # Add the nf2ff recording box
    return mesh

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

def run_simulation(FDTD, Sim_Path, sim_file, temp_file):
    """
    Run the simulation.
    """
    if not os.path.exists(sim_file):
        try:
            with open(temp_file, 'w') as f:
                f.write('Running')
            FDTD.Run(Sim_Path, verbose=0, cleanup=False)
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
    s11_dB = 20.0*np.log10(np.abs(s11))
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
        nf2ff_res_theta90 = nf2ff.CalcNF2FF(Sim_Path, f_res, 90, phi, center=np.array([0, 0, 0]) , read_cached=True, outfile='nf2ff_xy.h5')

        print('Radiated power: Prad = {:.2e} Watt'.format(nf2ff_res_theta90.Prad[0]))
        print('Directivity:    Dmax = {:.1f} ({:.1f} dBi)'.format(nf2ff_res_theta90.Dmax[0],
                                                                10 * np.log10(nf2ff_res_theta90.Dmax[0])))
        print('Efficiency:   nu_rad = {:.1f} %'.format(100 * nf2ff_res_theta90.Prad[0] / np.real(P_in[idx[0]])))

        print(f"Resonance frequency: {f_res/1e9} GHz")
        #s11 at closste to center frequency
        print(f"S11 at resonance frequency: {s11_dB[f_res_ind]} dB")

    print(f"S11 at center frequency{ s11_dB[freq == center_freq]} dB")

    if plot:
        # Create a figure with subplots, 2 rows and 2 columns
            fig, axs = plt.subplots(2, 2, figsize=(12, 10))

            # Plot feed point impedance in the top-left subplot (axs[0, 0])
            axs[0, 0].plot(freq / 1e6, np.real(Zin), 'k-', linewidth=2, label='Re{Zin}')
            axs[0, 0].plot(freq / 1e6, np.imag(Zin), 'r--', linewidth=2, label='Im{Zin}')
            # Add a horizontal line at 50 Ohms
            axs[0, 0].axhline(parameters["feed_R"], color='blue', linestyle='--', linewidth=1, label=f'{parameters["feed_R"]} Ohm')
            axs[0, 0].axvline(parameters["min_freq"]/1e6, color='green', linestyle='--', linewidth=1, label=f'{parameters["min_freq"]/1e6} MHz')
            axs[0, 0].axvline(parameters["center_freq"]/1e6, color='green', linestyle='--', linewidth=1, label=f'{parameters["center_freq"]/1e6} MHz')
            axs[0, 0].axvline(parameters["max_freq"]/1e6, color='green', linestyle='--', linewidth=1, label=f'{parameters["max_freq"]/1e6} MHz')


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
            axs[0, 1].axvline(parameters["min_freq"]/1e6, color='green', linestyle='--', linewidth=1, label=f'{parameters["min_freq"]/1e6} MHz')
            axs[0, 1].axvline(center_freq/1e6, color='green', linestyle='--', linewidth=1, label=f'{center_freq/1e6} MHz')
            axs[0, 1].axvline(parameters["max_freq"]/1e6, color='green', linestyle='--', linewidth=1, label=f'{parameters["max_freq"]/1e6} MHz')

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
                E_norm = 20.0 * np.log10(nf2ff_res_theta90.E_norm / np.max(nf2ff_res_theta90.E_norm)) + nf2ff_res_theta90.Dmax
                ax_polar2.plot(np.deg2rad(phi), 10 ** (np.squeeze(E_norm) / 20), linewidth=2, label='xy-plane')
                ax_polar2.grid(True)
                ax_polar2.set_xlabel('phi (deg)')
                plt.suptitle('IFA Antenna Pattern\nFrequency: {:.2f} GHz'.format(f_res / 1e9), fontsize=14)
                ax_polar2.legend(loc=3)

            # Prepare the parameters text with table-like formatting
            parameters_text = (
                f"{'IFA Parameters:':<30}\n"
                f"{'ifa_h:':<20}{parameters['ifa_h']:>10.3f} mm "
                f"{'ifa_l:':<20}{parameters['ifa_l']:>10.3f} mm\n"
                f"{'ifa_w1:':<20}{parameters['ifa_w1']:>10.3f} mm "
                f"{'ifa_w2:':<20}{parameters['ifa_w2']:>10.3f} mm\n"
                f"{'ifa_wf:':<20}{parameters['ifa_wf']:>10.3f} mm "
                f"{'ifa_fp:':<20}{parameters['ifa_fp']:>10.3f} mm\n"
                f"{'ifa_e:':<20}{parameters['ifa_e']:>10.3f} mm "
                f"{'ifa_h + ifa_e:':<20}{parameters['ifa_h'] + parameters['ifa_e']:>10.3f} mm\n\n"
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

def ifa_simulation(Sim_CSX='IFA.xml',
                   showCad=True,
                   post_proc_only=False,
                   unit=1e-3,
                   substrate_width=21,
                   substrate_length=83.15,
                   substrate_thickness=1.5,
                   gndplane_position=0,
                   substrate_cells=4,
                   ifa_h=14.054,
                   ifa_l=19.500,
                   ifa_w1=0.608,
                   ifa_w2=0.400,
                   ifa_wf=0.762,
                   ifa_fp=5.368,
                   ifa_e=0.5,
                   mifa_meander=3.0,
                   mifa_tipdistance=2.0,
                   mifa_meander_edge_distance=3.0,
                   substrate_epsR=4.5,
                   feed_R=50,
                   min_freq=2.4e9,
                   center_freq=2.45e9,
                   max_freq=2.5e9,
                   min_size=0.2,#minimum automesh size
                   max_size=4.0,#maximum automesh size
                   fc=1.0e9,  # 20 dB corner frequency
                   max_timesteps=600000,
                   override_min_global_grid=None,
                   plot=True,
                   delete_simulation_files=True):
    #############################################################################
    #                substrate_width
    #  _______________________________________________    __ substrate_thickness
    # | A                        ifa_l(total length)  |\  __
    # | |ifa_e         _____________________     w2   | |
    # | |             |    ___  _________  |____|  |  | |
    # | |       ifa_h |   |   ||         |_________|  | | ______
    # | |             |   |   ||         mifa_meander | ||mifa_edgedistance (minimum value)
    # |_V_____________|___|___||______________________| ||______
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
    f0=center_freq  # center frequency
    unit = 1e-3  # all lengths in mm
    # Derived parameter
    substrate_kappa = 1e-3 * 2 * pi * 2.45e9 * EPS0 * substrate_epsR

    params_tuple = (
        Sim_CSX,
        unit,
        substrate_width,
        substrate_length,
        substrate_thickness,
        gndplane_position,
        substrate_cells,
        ifa_h,
        ifa_l,
        ifa_w1,
        ifa_w2,
        ifa_wf,
        ifa_fp,
        ifa_e,
        mifa_meander,
        mifa_tipdistance,
        mifa_meander_edge_distance,
        substrate_epsR,
        feed_R,
        min_size,
        max_size,
        max_timesteps,
        override_min_global_grid,
        f0,
        fc,
    )

    params_str = '_'.join(map(str, params_tuple))
    # Generate a SHA256 hash of the parameters string
    params_hash = hashlib.sha256(params_str.encode('utf-8')).hexdigest()
    # Use the first 8 characters of the hash for brevity
    hash_prefix = params_hash[:12]
    # Create the simulation path using the hash
    Sim_Path = os.path.join(base_path, f'tmp_IFA_{hash_prefix}')

    # Simulation box size
    SimBox = np.array([substrate_width * 2, (substrate_length+2*ifa_h) * 2, 150])

    # Initialize simulation
    FDTD, CSX = initialize_simulation(f0, fc, max_timesteps)
    mesh = CSX.GetGrid()
    setup_mesh(mesh, SimBox, unit)

    # Create substrate
    create_substrate(CSX, substrate_width, substrate_length, substrate_thickness, ifa_h, ifa_e, substrate_epsR, substrate_kappa, substrate_cells, unit, mesh)

    # Create ground plane
    create_ground_plane(CSX, substrate_width, substrate_length, substrate_thickness, ifa_e, gndplane_position, mesh, min_size, max_size)

    # Create IFA
    parameters = {
        'unit': unit,
        'substrate_width': substrate_width,
        'substrate_length': substrate_length,
        'substrate_thickness': substrate_thickness,
        "substrate_epsR": substrate_epsR,
        "substrate_kappa": substrate_kappa,
        'ifa_h': ifa_h,
        'ifa_l': ifa_l,
        'ifa_w1': ifa_w1,
        'ifa_w2': ifa_w2,
        'ifa_wf': ifa_wf,
        'ifa_fp': ifa_fp,
        'ifa_e': ifa_e,
        'mifa_meander': mifa_meander,
        'mifa_tipdistance': mifa_tipdistance,
        'mifa_meander_edge_distance': mifa_meander_edge_distance,
        'gndplane_position': gndplane_position,
        'min_size': min_size,
        'max_size': max_size,
        "feed_R": feed_R,
        "min_freq": min_freq,
        "max_freq": max_freq,
        "center_freq": center_freq,
        "fc": fc, 
    }
    f0 = center_freq
    fc = parameters["fc"]
    create_ifa(CSX, mesh, parameters)

    # Apply the excitation & resistor as a current source
    tl = np.array([-substrate_width/2+ifa_e+ifa_fp, substrate_length / 2 - ifa_e, substrate_thickness])  # translation vector
    via_diameter=0.3
    if gndplane_position != 0:
        feedpoint = -via_diameter #set feedpoint to groundplane edge
    else:
        feedpoint = ifa_wf

    port = add_feed(FDTD, CSX, mesh, tl, ifa_wf, feedpoint, gndplane_position, feed_R, min_size, max_size)

    # Finalize the mesh
    nf2ff = FDTD.CreateNF2FFBox()
    finalize_mesh(mesh, min_size, f0, fc, unit, override_min_global_grid, {"x": [], "y": [], "z": [0, substrate_thickness+gndplane_position,substrate_thickness]})

    temp_file = os.path.join(Sim_Path, 'incomplete_run.flag')

    if os.path.exists(Sim_Path):
        if os.path.exists(temp_file):
            import shutil
            print(f"Cleaning up incomplete run: {Sim_Path}")
            shutil.rmtree(Sim_Path)

    if not os.path.exists(Sim_Path):
        print(f"Creating directory {Sim_Path}")
        os.mkdir(Sim_Path)

    CSX_file = os.path.join(Sim_Path, Sim_CSX)
    CSX.Write2XML(CSX_file)

    # Show the structure
    if showCad:
        print("showing cad")
        print(f"csxcad location{csxcad_location}")
        print(f"csxfile location {CSX_file}")
        os.system( csxcad_location + ' "{}"'.format(CSX_file))

    sim_file = os.path.join(Sim_Path, 'complete_run.flag')

    if not post_proc_only and not os.path.exists(sim_file):
        run_simulation(FDTD, Sim_Path, sim_file, temp_file)
    
    if os.path.exists(sim_file):
        # Post-processing & plotting
        freq = np.linspace(max(0, f0 - fc), f0 + fc, 501)
        freq, s11_dB, Zin, P_in = post_process_results(Sim_Path, port, freq, delete_simulation_files, plot, center_freq, nf2ff,parameters)

        return freq, s11_dB, Zin, P_in , hash_prefix
    return None, None, None, None, None

if __name__ == "__main__":
    ifa_simulation()
