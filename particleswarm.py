from ifa_test import ifa_simulation
import logging
import numpy as np
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx

# Function to evaluate the swarm
def evaluation_fun(swarm):
    # Extract parameter arrays from swarm
    ifa_h = swarm[:, 0]
    ifa_l = swarm[:, 1]
    ifa_fp = swarm[:, 2]
    ifa_w = swarm[:, 3]

    # Constants for the simulation
    Sim_CSX = 'IFA.xml'
    showCad = False
    post_proc_only = False
    unit = 1e-3
    substrate_width = 21
    substrate_length = 83.15
    substrate_thickness = 1.5
    gndplane_position = 0
    substrate_cells = 4
    ifa_e = 0.5
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.2  # Minimum automesh size
    plot = False

    # Prepare an array to store the fitness results
    fitness_results = np.zeros(swarm.shape[0])

    # Loop over each particle in the swarm
    for i in range(swarm.shape[0]):
        # Call the simulation for each particle
        try:
            freq, s11_dB, Zin, P_in = ifa_simulation(
                Sim_CSX=Sim_CSX,
                showCad=showCad,
                post_proc_only=post_proc_only,
                unit=unit,
                substrate_width=substrate_width,
                substrate_length=substrate_length,
                substrate_thickness=substrate_thickness,
                gndplane_position=gndplane_position,
                substrate_cells=substrate_cells,
                ifa_h=ifa_h[i],
                ifa_l=ifa_l[i],
                ifa_fp=ifa_fp[i],
                ifa_w1=ifa_w[i],
                ifa_w2=ifa_w[i],
                ifa_wf=ifa_w[i],
                ifa_e=ifa_e,
                substrate_epsR=substrate_epsR,
                feed_R=feed_R,
                min_freq=min_freq,
                center_freq=center_freq,
                max_freq=max_freq,
                plot=plot,
                min_size=min_size
            )
            
            # Assume that the function returns arrays and we need the index where frequency equals the center frequency
            fitness_results[i] = s11_dB[freq == center_freq]  # Assign fitness based on s11 at the center frequency

        except exception as e:
            logging.warning(f"IndexError during port calculation: {e}. Continuing with the next simulation.")
            fitness_results[i] = 0
    
        # Log the parameters and the corresponding fitness value
        logging.info(f"Iteration: {i}, ifa_h: {ifa_h[i]}, ifa_l: {ifa_l[i]}, ifa_fp: {ifa_fp[i]}, ifa_w: {ifa_w[i]}, fitness: {fitness_results[i]}")


    return fitness_results

if __name__ == "__main__":
    # Setup file logging
    logging.basicConfig(filename='logs\\pso_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

    # Define bounds
    bounds = (np.array([1, 10, 1, 0.4]), np.array([14, 19.5, 10, 1]))
    dimensions = 4

    # Setup optimizer
    options = {'c1': 1.5, 'c2': 1.5, 'w': 0.9}
    optimizer = ps.single.GlobalBestPSO(n_particles=50, dimensions=4, options=options, bounds=bounds)

    # Optimize with a step callback to log progress
    cost, pos = optimizer.optimize(evaluation_fun, iters=100)
